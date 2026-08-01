#!/bin/bash
# Definicion a priori del conjunto flagelar desde la anotacion del genoma.
# IMPORTANTE: este script no consulta ningun dato de expresion.
# Es la evidencia de que el conjunto de 31 genes se definio antes del analisis.

GFF=ref/ncbi_dataset/data/GCF_000015445.1/genomic.gff

# Tabla de todas las CDS: locus_tag, inicio, fin, hebra, producto
awk -F'\t' '$3=="CDS"{
  lt=""; pr="";
  if (match($9,/locus_tag=[^;]+/)) lt=substr($9,RSTART+10,RLENGTH-10);
  if (match($9,/product=[^;]+/))   pr=substr($9,RSTART+8,RLENGTH-8);
  print lt"\t"$4"\t"$5"\t"$7"\t"pr
}' $GFF | sort -u > cds_todos.tsv

# Conjunto flagelar por patron sobre el campo product
grep -Ei 'flagell|chemotax|hook|motility|basal.?body' cds_todos.tsv \
  | cut -f1 | sort -u > flag_gff.tsv

echo "CDS totales:      $(wc -l < cds_todos.tsv)"
echo "Genes flagelares: $(wc -l < flag_gff.tsv)"
