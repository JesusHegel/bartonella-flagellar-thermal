#!/bin/bash
# Busqueda por homologia de componentes del sistema de quimiotaxis.
# Consultas: 94 proteinas de cuatro alfaproteobacterias con sistema completo.

for A in GCF_037023865.1 GCF_000022005.1 GCF_000092025.1 GCF_003324715.1; do
  datasets download genome accession $A --include protein --filename ${A}.zip
  unzip -q -o ${A}.zip -d ${A}
done

# Extraccion de proteinas Che de cada organismo
for D in GCF_*/; do
  F=$(find $D -name "protein.faa" | head -1)
  python3 - "$F" <<'PY' >> che_todos.faa
import sys,re
pat=re.compile(r'chemotaxis protein Che[ABDRWYZ]|protein-glutamate|CheR family|CheB family|CheW|CheA|methyl-accepting chemotaxis',re.I)
k=False
for l in open(sys.argv[1]):
    if l.startswith('>'): k=bool(pat.search(l))
    if k: sys.stdout.write(l)
PY
done

blastp -query che_todos.faa -db kc583 -max_hsps 1 -max_target_seqs 1 \
  -evalue 1e-3 -outfmt '6 qseqid sseqid pident length qlen evalue' | sort -k6,6g
