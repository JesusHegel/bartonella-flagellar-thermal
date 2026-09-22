#!/bin/bash
# Busqueda por homologia del sistema de quimiotaxis y estructura del 5'UTR
# de la flagelina. Requiere BLAST+, NCBI datasets y ViennaRNA.
# Se ejecuta desde ~/bbrna. Depende de s07 (usa v2/faa/kc583.faa).
set -u
cd ~/bbrna
DATASETS=$(ls $HOME/miniforge3/envs/*/bin/datasets $HOME/miniforge3/bin/datasets 2>/dev/null | head -1)
BLASTP=$(ls $HOME/miniforge3/envs/*/bin/blastp 2>/dev/null | head -1)
MAKEDB=$(ls $HOME/miniforge3/envs/*/bin/makeblastdb 2>/dev/null | head -1)
RNAFOLD=$(ls $HOME/miniforge3/envs/*/bin/RNAfold 2>/dev/null | head -1)

mkdir -p v2/che && cd v2/che
for A in GCF_037023865.1 GCF_000022005.1 GCF_000092025.1 GCF_003324715.1; do
  [ -d $A ] || { $DATASETS download genome accession $A --include protein --filename $A.zip \
                 && unzip -q -o $A.zip -d $A; }
done

: > che_todos.faa
for D in GCF_*/; do
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

$MAKEDB -in ../faa/kc583.faa -dbtype prot -out kc583db >/dev/null 2>&1
echo; echo "=== HOMOLOGOS Che EN KC583 (e<1e-3) ==="
$BLASTP -query che_todos.faa -db kc583db -max_hsps 1 -max_target_seqs 1 \
  -evalue 1e-3 -num_threads 4 -outfmt '6 qseqid sseqid pident length qlen evalue' 2>/dev/null \
  | sort -k6,6g | tee ../chk_quimiotaxis.tsv
echo "Aciertos totales: $(wc -l < ../chk_quimiotaxis.tsv)"
echo "Proteinas KC583 distintas alcanzadas: $(cut -f2 ../chk_quimiotaxis.tsv | sort -u | wc -l)"

echo; echo "=== ESTRUCTURA DEL 5'UTR DE LA FLAGELINA ==="
cd ~/bbrna
# IMPORTANTE: el genoma completo, no cds_from_genomic.fna
G=$(ls ref/ncbi_dataset/data/GCF_000015445.1/*_genomic.fna | grep -v cds_from | head -1)
echo "Genoma: $G"
python3 - "$G" <<'PY' > v2/utr_flagelina.fa
import sys
s=''.join(l.strip() for l in open(sys.argv[1]) if not l.startswith('>'))
sub=s[1078099:1078240]
comp={'A':'T','T':'A','G':'C','C':'G','N':'N'}
print('>5UTR_flagellin_RS05045_141nt')
print(''.join(comp[c] for c in reversed(sub.upper())))
PY
$RNAFOLD --noPS < v2/utr_flagelina.fa | tee v2/utr_mfe.txt

python3 - "$G" <<'PY' > v2/control.fa
import sys,random
random.seed(7)
s=''.join(l.strip() for l in open(sys.argv[1]) if not l.startswith('>')).upper()
for i in range(100):
    p=random.randint(0,len(s)-142)
    print(f'>ctrl{i}'); print(s[p:p+141])
PY
$RNAFOLD --noPS < v2/control.fa | grep -oP '\(\s*-?[0-9]+\.[0-9]+\)$' | tr -d '() ' > v2/mfe_ctrl.txt
python3 -c "
import statistics as st
v=sorted(float(x) for x in open('v2/mfe_ctrl.txt'))
u=float(open('v2/utr_mfe.txt').read().strip().split('(')[-1].rstrip(')').strip())
print('UTR flagelina MFE : %.1f' % u)
print('Control n=%d  media %.1f  mediana %.1f' % (len(v),st.mean(v),st.median(v)))
print('Controles mas estructurados que el UTR: %.0f %%' % (100*sum(1 for x in v if x<u)/len(v)))"
