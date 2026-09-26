#!/bin/bash
# 06. Panel de genero y ortologia con la cepa peruana USM-LMMB07.
#
# Busca por BLASTp los 31 genes flagelares (conjunto del script 01) en siete
# especies de Bartonella, y calcula ortologos reciprocos entre KC583 y
# USM-LMMB07. En las especies donde no estan los 31, busca de vuelta en KC583
# cada proteina encontrada, para saber si es el gen flagelar o una paraloga. Descarga los proteomas del NCBI (necesita internet) y guarda la
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

echo; echo "=== ESPECIES SIN LOS 31: cada acierto se busca de vuelta en KC583 ==="
# En las especies donde no estan los 31 genes, cada proteina encontrada se busca
# de vuelta en KC583. Si su mejor coincidencia es el mismo gen flagelar, el
# acierto es un ortologo reciproco. Si es otra proteina de KC583, el acierto es
# una proteina emparentada (paraloga), no el gen flagelar.
# Se piden las 5 mejores coincidencias y se elige la de mayor puntaje (bitscore),
# porque -max_target_seqs 1 no garantiza devolver la mejor.
while read A N; do echo "$N $(proteoma $A)"; done < nombres.txt > rutas_proteomas.txt
: > ida.tsv
while read A N; do
  [ "$(awk -v n=$N '$2==n' "$GC/panel_genero_blastp.tsv" | wc -l)" -lt 31 ] || continue
  blastp -query flag31.faa -db db/$N -max_hsps 1 -max_target_seqs 5 -evalue 1e-5 -num_threads 4 \
    -outfmt '6 qseqid sseqid pident qcovs bitscore' 2>/dev/null \
    | sort -t$'\t' -k1,1 -k5,5gr -k2,2 | awk -v n=$N -F'\t' '!v[$1]++ {print n"\t"$0}' >> ida.tsv
done < nombres.txt

python3 - <<'PY'
import re
def leer(fa):   # identificador -> (producto, secuencia)
    d, cur = {}, None
    for l in open(fa):
        if l.startswith('>'):
            cur, _, desc = l[1:].rstrip('\n').partition(' ')
            desc = re.sub(r'\s*\[[^\]]*\]$', '', desc).replace('MULTISPECIES: ', '')
            d[cur] = [desc, []]
        else:
            d[cur][1].append(l.strip())
    return {k: (v[0], ''.join(v[1])) for k, v in d.items()}
rutas = dict(l.split() for l in open('rutas_proteomas.txt'))
prot = {}
with open('aciertos.faa', 'w') as f:
    for l in open('ida.tsv'):
        n, q, s = l.rstrip('\n').split('\t')[:3]
        prot.setdefault(n, leer(rutas[n]))
        f.write('>%s|%s\n%s\n' % (n, s, prot[n][s][1]))
PY
blastp -query aciertos.faa -db db/bacilliformis -max_hsps 1 -max_target_seqs 5 -evalue 1e-5 -num_threads 4 \
  -outfmt '6 qseqid sseqid bitscore' 2>/dev/null \
  | sort -t$'\t' -k1,1 -k3,3gr -k2,2 | awk -F'\t' '!v[$1]++' > vuelta.tsv

python3 - "$GC/panel_genero_busqueda_de_vuelta.tsv" <<'PY'
import re, sys
def productos(fa):
    d = {}
    for l in open(fa):
        if l.startswith('>'):
            i, _, desc = l[1:].rstrip('\n').partition(' ')
            d[i] = re.sub(r'\s*\[[^\]]*\]$', '', desc).replace('MULTISPECIES: ', '')
    return d
rutas = dict(l.split() for l in open('rutas_proteomas.txt'))
kc = productos(rutas['bacilliformis'])
vuelta = dict(l.rstrip('\n').split('\t')[:2] for l in open('vuelta.tsv'))
cab = ['especie', 'proteina_kc583', 'producto_kc583', 'acierto', 'producto_acierto',
       'identidad', 'cobertura', 'vuelta_en_kc583', 'producto_vuelta', 'reciproco']
filas, prods = [], {}
for l in open('ida.tsv'):
    n, q, s, pid, cov = l.rstrip('\n').split('\t')[:5]
    prods.setdefault(n, productos(rutas[n]))
    v = vuelta.get(n + '|' + s, '')
    filas.append([n, q, kc[q], s, prods[n][s], pid, cov, v or 'ninguna', kc.get(v, ''),
                  'si' if v == q else 'no'])
with open(sys.argv[1], 'w') as f:
    f.write('\t'.join(cab) + '\n')
    for r in filas: f.write('\t'.join(r) + '\n')
for r in filas:
    print('%-12s %s (%s) -> %s (%s), identidad %s %%, cobertura %s %%; de vuelta: %s (%s) -> reciproco: %s'
          % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9]))
PY
