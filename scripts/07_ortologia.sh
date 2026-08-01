#!/bin/bash
# Ortologia entre B. bacilliformis KC583 y USM-LMMB07 por mejores
# coincidencias reciprocas (BLASTp), reteniendo el mejor bitscore por consulta.

makeblastdb -in kc583.faa -dbtype prot -out kcdb > /dev/null 2>&1
makeblastdb -in usm.faa   -dbtype prot -out usmdb > /dev/null 2>&1

blastp -query usm.faa -db kcdb -evalue 1e-10 -max_hsps 1 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen bitscore' > usm_vs_kc.tsv
blastp -query kc583.faa -db usmdb -evalue 1e-10 -max_hsps 1 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen bitscore' > kc_vs_usm.tsv

sort -k1,1 -k7,7gr usm_vs_kc.tsv | awk '!v[$1]++' > best_usm_kc.tsv
sort -k1,1 -k7,7gr kc_vs_usm.tsv | awk '!v[$1]++' > best_kc_usm.tsv
awk 'NR==FNR{m[$1]=$2;next} m[$2]==$1' best_kc_usm.tsv best_usm_kc.tsv > rbh.tsv

echo "Pares reciprocos: $(wc -l < rbh.tsv)"
echo "Consultas totales: $(grep -c '>' usm.faa)"
