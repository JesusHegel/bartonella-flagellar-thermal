#!/bin/bash
# 06. Panel de genero y ortologia con la cepa peruana USM-LMMB07.
#
# Busca por BLASTp los 31 genes flagelares (conjunto del script 01) en siete
# especies de Bartonella, y calcula ortologos reciprocos entre KC583 y
# USM-LMMB07. Descarga los proteomas del NCBI (necesita internet) y guarda la
# huella md5 de cada uno, para saber si el NCBI los cambia en el futuro.
set -uo pipefail
source scripts/rutas.sh
RAIZ=$(pwd)
GC=$RAIZ/$DIR_RES/genomica_comparada
AN=$RAIZ/$DIR_RES/anotacion

mkdir -p "$DIR_DESC/proteomas" && cd "$DIR_DESC/proteomas"

# KC583 se toma de la referencia incluida en el repositorio.
cp "$RAIZ/$REF/protein.faa" kc583.faa
echo "GCF_000015445.1  $(grep -c '>' kc583.faa) proteinas (referencia local)"

for A in GCF_000253015.1 GCF_039555305.1 GCF_001281405.1 \
         GCF_000196435.1 GCF_009936175.1 GCF_019930925.1 GCF_001624625.1; do
  [ -d $A ] || { datasets download genome accession $A --include protein --filename $A.zip \
                 && unzip -q -o $A.zip -d $A; }
  echo "$A  $(grep -c '>' "$(find $A -name protein.faa | head -1)") proteinas"
done
cp "$(find GCF_001624625.1 -name protein.faa | head -1)" usm.faa

cat > nombres.txt <<'EONOM'
kc583 bacilliformis
GCF_000253015.1 clarridgeiae
GCF_039555305.1 schoenbuchensis
GCF_001281405.1 ancashensis
GCF_000196435.1 tribocorum
GCF_009936175.1 quintana
GCF_019930925.1 henselae
EONOM
proteoma() { [ "$1" = kc583 ] && echo kc583.faa || find "$1" -name protein.faa | head -1; }

: > "$GC/md5_proteomas_descargados.txt"
while read A N; do echo "$(md5sum < "$(proteoma $A)" | cut -d' ' -f1)  $N" >> "$GC/md5_proteomas_descargados.txt"; done < nombres.txt
echo "$(md5sum < usm.faa | cut -d' ' -f1)  USM-LMMB07" >> "$GC/md5_proteomas_descargados.txt"

# Proteinas de los 31 genes flagelares
cut -f1 "$AN/conjunto_flagelar.tsv" | sort -u > lt31.txt
awk -F'\t' 'NR==FNR{k[$1];next} $1 in k {print $2}' lt31.txt "$AN/locus_a_proteina.tsv" | sort -u > wp31.txt
python3 - <<'PY' > flag31.faa
want={l.strip() for l in open('wp31.txt') if l.strip()}
keep=False
for l in open('kc583.faa'):
    if l.startswith('>'): keep = l[1:].split()[0] in want
    if keep: print(l, end='')
PY
echo "Secuencias en flag31.faa: $(grep -c '>' flag31.faa)"

mkdir -p db
while read A N; do
  makeblastdb -in "$(proteoma $A)" -dbtype prot -out db/$N >/dev/null 2>&1
done < nombres.txt

echo; echo "=== PANEL DE GENERO (BLASTp, e<1e-5) ==="
: > "$GC/panel_genero_blastp.tsv"
while read A N; do
  blastp -query flag31.faa -db db/$N -max_hsps 1 -max_target_seqs 1 \
    -evalue 1e-5 -num_threads 4 -outfmt '6 qseqid pident' 2>/dev/null \
    | sort -u -k1,1 | awk -v n=$N '{print $1"\t"n"\t"$2}' >> "$GC/panel_genero_blastp.tsv"
  printf "%-18s %s de 31\n" "$N" "$(awk -v n=$N '$2==n' "$GC/panel_genero_blastp.tsv" | wc -l)"
done < nombres.txt

echo; echo "=== ORTOLOGIA KC583 vs USM-LMMB07 (coincidencias reciprocas) ==="
makeblastdb -in kc583.faa -dbtype prot -out db/kcdb  >/dev/null 2>&1
makeblastdb -in usm.faa   -dbtype prot -out db/usmdb >/dev/null 2>&1
blastp -query usm.faa   -db db/kcdb  -evalue 1e-10 -max_hsps 1 -num_threads 4 \
       -outfmt '6 qseqid sseqid pident bitscore' > usm_vs_kc.tsv
blastp -query kc583.faa -db db/usmdb -evalue 1e-10 -max_hsps 1 -num_threads 4 \
       -outfmt '6 qseqid sseqid pident bitscore' > kc_vs_usm.tsv
sort -k1,1 -k4,4gr usm_vs_kc.tsv | awk '!v[$1]++' > best_usm_kc.tsv
sort -k1,1 -k4,4gr kc_vs_usm.tsv | awk '!v[$1]++' > best_kc_usm.tsv
awk 'NR==FNR{m[$1]=$2;next} m[$2]==$1' best_kc_usm.tsv best_usm_kc.tsv > "$GC/ortologos_reciprocos_KC583_USM-LMMB07.tsv"
echo "Pares reciprocos : $(wc -l < "$GC/ortologos_reciprocos_KC583_USM-LMMB07.tsv")"
echo "Flagelares con ortologo reciproco: $(awk 'NR==FNR{k[$1];next} $2 in k' wp31.txt "$GC/ortologos_reciprocos_KC583_USM-LMMB07.tsv" | wc -l) de 31"
