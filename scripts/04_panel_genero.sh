#!/bin/bash
# Distribucion del regulon flagelar en siete especies del genero Bartonella.
# La identidad taxonomica de cada genoma se verifica por consulta a NCBI.

mkdir -p genero && cd genero

for A in GCF_000015445.1 GCF_000253015.1 GCF_039555305.1 \
         GCF_001281405.1 GCF_000196435.1 GCF_009936175.1 GCF_019930925.1; do
  datasets download genome accession $A --include protein --filename ${A}.zip
  unzip -q -o ${A}.zip -d ${A}
done

printf '%s\n' \
"GCF_000015445.1 bacilliformis" "GCF_000253015.1 clarridgeiae" \
"GCF_039555305.1 schoenbuchensis" "GCF_001281405.1 ancashensis" \
"GCF_000196435.1 tribocorum" "GCF_009936175.1 quintana" \
"GCF_019930925.1 henselae" > nombres.txt

# Bases de datos por especie
mkdir -p db
while read A N; do
  F=$(find $A -name "protein.faa" | head -1)
  makeblastdb -in $F -dbtype prot -out db/$N > /dev/null 2>&1
done < nombres.txt

# BLASTp de los 31 flagelares (conjunto definido a priori) contra cada proteoma
while read A N; do
  C=$(blastp -query flag31.faa -db db/$N -max_hsps 1 -max_target_seqs 1 \
      -evalue 1e-5 -outfmt '6 qseqid' | sort -u | wc -l)
  printf "%-16s %s de 31\n" "$N" "$C"
done < nombres.txt
