#!/bin/bash
# 01. Anotacion del genoma KC583 y conjunto flagelar definido a priori.
#
# Lee solo el archivo de anotacion del NCBI (datos/referencia/genomic.gff).
# No consulta ningun dato de expresion: el conjunto flagelar se fija antes
# de mirar los resultados.
set -euo pipefail
source scripts/rutas.sh

GFF=$REF/genomic.gff
A=$DIR_RES/anotacion
[ -f "$GFF" ] || { echo "ERROR: falta $GFF"; exit 1; }

# Una fila por CDS: locus, inicio, fin, hebra, producto, proteina
awk -F'\t' '$3=="CDS"{
  lt="";pr="";pid="";
  if (match($9,/locus_tag=[^;]+/))  lt=substr($9,RSTART+10,RLENGTH-10);
  if (match($9,/product=[^;]+/))    pr=substr($9,RSTART+8,RLENGTH-8);
  if (match($9,/protein_id=[^;]+/)) pid=substr($9,RSTART+11,RLENGTH-11);
  gsub(/%2C/,",",pr); gsub(/%3B/,";",pr);
  print lt"\t"$4"\t"$5"\t"$7"\t"pr"\t"pid
}' "$GFF" | sort -u > "$DIR_INT/cds_completo.tsv"

cut -f1-5 "$DIR_INT/cds_completo.tsv"                                   > "$A/anotacion_cds.tsv"
awk -F'\t' '$6!=""{print $1"\t"$6}' "$DIR_INT/cds_completo.tsv" | sort -u > "$A/locus_a_proteina.tsv"
grep -Ei 'flagell|chemotax|hook|motility|basal.?body' "$A/anotacion_cds.tsv" | sort -u > "$A/conjunto_flagelar.tsv"
cut -f2 "$A/locus_a_proteina.tsv" | sort | uniq -d                       > "$A/proteinas_duplicadas.txt"

echo "CDS totales            : $(wc -l < "$A/anotacion_cds.tsv")"
echo "Con protein_id         : $(wc -l < "$A/locus_a_proteina.tsv")"
echo "Accesiones unicas      : $(cut -f2 "$A/locus_a_proteina.tsv" | sort -u | wc -l)"
echo "Accesiones duplicadas  : $(wc -l < "$A/proteinas_duplicadas.txt")"
echo "Genes flagelares       : $(wc -l < "$A/conjunto_flagelar.tsv")"
echo
echo "--- conjunto flagelar ---"
cut -f1,5 "$A/conjunto_flagelar.tsv"
