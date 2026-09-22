#!/bin/bash
# Panel de genero y ortologia con la cepa peruana USM-LMMB07.
# Requiere NCBI datasets y BLAST+. Se ejecuta desde ~/bbrna.
# Depende de s01_anotacion.sh (usa salida/flag_gff.tsv y salida/lt2prot.tsv).
set -u
cd ~/bbrna
DATASETS=$(ls $HOME/miniforge3/envs/*/bin/datasets $HOME/miniforge3/bin/datasets 2>/dev/null | head -1)
BLASTP=$(ls $HOME/miniforge3/envs/*/bin/blastp 2>/dev/null | head -1)
MAKEDB=$(ls $HOME/miniforge3/envs/*/bin/makeblastdb 2>/dev/null | head -1)

mkdir -p salida/faa && cd salida/faa
for A in GCF_000015445.1 GCF_000253015.1 GCF_039555305.1 GCF_001281405.1 \
         GCF_000196435.1 GCF_009936175.1 GCF_019930925.1 GCF_001624625.1; do
  [ -d $A ] || { $DATASETS download genome accession $A --include protein --filename $A.zip \
                 && unzip -q -o $A.zip -d $A; }
  echo "$A  $(grep -c '>' "$(find $A -name protein.faa | head -1)") proteinas"
done

cat > nombres.txt <<'EONOM'
GCF_000015445.1 bacilliformis
GCF_000253015.1 clarridgeiae
GCF_039555305.1 schoenbuchensis
GCF_001281405.1 ancashensis
GCF_000196435.1 tribocorum
GCF_009936175.1 quintana
GCF_019930925.1 henselae
EONOM

cp "$(find GCF_000015445.1 -name protein.faa | head -1)" kc583.faa
cp "$(find GCF_001624625.1 -name protein.faa | head -1)" usm.faa

# Proteinas de los 31 genes flagelares (conjunto definido a priori en s01)
cut -f1 ../flag_gff.tsv | sort -u > lt31.txt
awk -F'\t' 'NR==FNR{k[$1];next} $1 in k {print $2}' lt31.txt ../lt2prot.tsv | sort -u > wp31.txt
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
  $MAKEDB -in "$(find $A -name protein.faa | head -1)" -dbtype prot -out db/$N >/dev/null 2>&1
done < nombres.txt

echo; echo "=== PANEL DE GENERO (BLASTp, e<1e-5) ==="
: > ../panel.tsv
while read A N; do
  $BLASTP -query flag31.faa -db db/$N -max_hsps 1 -max_target_seqs 1 \
    -evalue 1e-5 -num_threads 4 -outfmt '6 qseqid pident' 2>/dev/null \
    | sort -u -k1,1 | awk -v n=$N '{print $1"\t"n"\t"$2}' >> ../panel.tsv
  printf "%-18s %s de 31\n" "$N" "$(awk -v n=$N '$2==n' ../panel.tsv | wc -l)"
done < nombres.txt

echo; echo "=== ORTOLOGIA KC583 vs USM-LMMB07 (coincidencias reciprocas) ==="
$MAKEDB -in kc583.faa -dbtype prot -out db/kcdb  >/dev/null 2>&1
$MAKEDB -in usm.faa   -dbtype prot -out db/usmdb >/dev/null 2>&1
$BLASTP -query usm.faa   -db db/kcdb  -evalue 1e-10 -max_hsps 1 -num_threads 4 \
        -outfmt '6 qseqid sseqid pident bitscore' > usm_vs_kc.tsv
$BLASTP -query kc583.faa -db db/usmdb -evalue 1e-10 -max_hsps 1 -num_threads 4 \
        -outfmt '6 qseqid sseqid pident bitscore' > kc_vs_usm.tsv
sort -k1,1 -k4,4gr usm_vs_kc.tsv | awk '!v[$1]++' > best_usm_kc.tsv
sort -k1,1 -k4,4gr kc_vs_usm.tsv | awk '!v[$1]++' > best_kc_usm.tsv
awk 'NR==FNR{m[$1]=$2;next} m[$2]==$1' best_kc_usm.tsv best_usm_kc.tsv > ../rbh.tsv
echo "Pares reciprocos : $(wc -l < ../rbh.tsv)"
echo "Flagelares con ortologo reciproco: $(awk 'NR==FNR{k[$1];next} $2 in k' wp31.txt ../rbh.tsv | wc -l) de 31"
