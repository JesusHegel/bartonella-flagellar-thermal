#!/bin/bash
# 07. Busqueda del sistema de quimiotaxis Che y estructura del 5'UTR de la flagelina.
#
# Quimiotaxis: proteinas Che de cuatro alfaproteobacterias con sistema completo,
# buscadas por BLASTp en el proteoma de KC583. Descarga los proteomas del NCBI.
# UTR: region intergenica de 141 nt aguas arriba de la flagelina (RS05045,
# hebra negativa), plegada con RNAfold y comparada con 100 fragmentos
# aleatorios del mismo genoma (semilla 7).
set -uo pipefail
source scripts/rutas.sh
RAIZ=$(pwd)
GC=$RAIZ/$DIR_RES/genomica_comparada

mkdir -p "$DIR_DESC/quimiotaxis" && cd "$DIR_DESC/quimiotaxis"
for A in GCF_037023865.1 GCF_000022005.1 GCF_000092025.1 GCF_003324715.1; do
  [ -d $A ] || { datasets download genome accession $A --include protein --filename $A.zip \
                 && unzip -q -o $A.zip -d $A; }
done
for A in GCF_037023865.1 GCF_000022005.1 GCF_000092025.1 GCF_003324715.1; do
  echo "$(md5sum < "$(find $A -name protein.faa | head -1)" | cut -d' ' -f1)  $A"
done > "$GC/md5_proteomas_quimiotaxis.txt"

: > che_todos.faa
for D in GCF_000022005.1 GCF_000092025.1 GCF_003324715.1 GCF_037023865.1; do   # orden alfabetico
  F=$(find "$D" -name protein.faa | head -1)
  python3 - "$F" <<'PY' >> che_todos.faa
import sys,re
pat=re.compile(r'chemotaxis protein Che[ABDRWYZ]|protein-glutamate|CheR family|CheB family|CheW|CheA|methyl-accepting chemotaxis',re.I)
k=False
for l in open(sys.argv[1]):
    if l.startswith('>'): k=bool(pat.search(l))
    if k: sys.stdout.write(l)
PY
done
echo "Consultas Che: $(grep -c '>' che_todos.faa)"

makeblastdb -in "$RAIZ/$REF/protein.faa" -dbtype prot -out kc583db >/dev/null 2>&1
echo; echo "=== HOMOLOGOS Che EN KC583 (e<1e-3) ==="
blastp -query che_todos.faa -db kc583db -max_hsps 1 -max_target_seqs 1 \
  -evalue 1e-3 -num_threads 4 -outfmt '6 qseqid sseqid pident length qlen evalue' 2>/dev/null \
  | sort -k6,6g | tee "$GC/homologos_quimiotaxis_che.tsv"
echo "Aciertos totales: $(wc -l < "$GC/homologos_quimiotaxis_che.tsv")"
echo "Proteinas KC583 distintas alcanzadas: $(cut -f2 "$GC/homologos_quimiotaxis_che.tsv" | sort -u | wc -l)"

echo; echo "=== ESTRUCTURA DEL 5'UTR DE LA FLAGELINA ==="
cd "$RAIZ"
G=$REF/GCF_000015445.1_ASM1544v1_genomic.fna      # genoma completo, no las CDS
echo "Genoma: $G"
# Flagelina RS05045: 1076957-1078099 (-); gen vecino RS05055 empieza en 1078241.
# La region intergenica es 1078100-1078240 (141 nt), en complemento reverso.
python3 - "$G" <<'PY' > "$GC/utr_flagelina.fa"
import sys
s=''.join(l.strip() for l in open(sys.argv[1]) if not l.startswith('>'))
sub=s[1078099:1078240]
comp={'A':'T','T':'A','G':'C','C':'G','N':'N'}
print('>5UTR_flagellin_RS05045_141nt')
print(''.join(comp[c] for c in reversed(sub.upper())))
PY
RNAfold --noPS < "$GC/utr_flagelina.fa" | tee "$GC/utr_flagelina_mfe.txt"

python3 - "$G" <<'PY' > "$DIR_INT/utr_controles.fa"
import sys,random
random.seed(7)
s=''.join(l.strip() for l in open(sys.argv[1]) if not l.startswith('>')).upper()
for i in range(100):
    p=random.randint(0,len(s)-142)
    print(f'>ctrl{i}'); print(s[p:p+141])
PY
RNAfold --noPS < "$DIR_INT/utr_controles.fa" | grep -oP '\(\s*-?[0-9]+\.[0-9]+\)$' | tr -d '() ' > "$GC/utr_controles_mfe.txt"
python3 - "$GC/utr_controles_mfe.txt" "$GC/utr_flagelina_mfe.txt" <<'PY'
import sys, statistics as st
v=sorted(float(x) for x in open(sys.argv[1]))
u=float(open(sys.argv[2]).read().strip().split('(')[-1].rstrip(')').strip())
print('UTR flagelina MFE : %.1f' % u)
print('Control n=%d  media %.1f  mediana %.1f' % (len(v),st.mean(v),st.median(v)))
print('Controles mas estructurados que el UTR: %.0f %%' % (100*sum(1 for x in v if x<u)/len(v)))
PY
