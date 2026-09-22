#!/bin/bash
set -euo pipefail
cd ~/bbrna
GFF=ref/ncbi_dataset/data/GCF_000015445.1/genomic.gff
[ -f "$GFF" ] || { echo "FALTA $GFF"; exit 1; }

awk -F'\t' '$3=="CDS"{
  lt="";pr="";pid="";
  if (match($9,/locus_tag=[^;]+/))  lt=substr($9,RSTART+10,RLENGTH-10);
  if (match($9,/product=[^;]+/))    pr=substr($9,RSTART+8,RLENGTH-8);
  if (match($9,/protein_id=[^;]+/)) pid=substr($9,RSTART+11,RLENGTH-11);
  gsub(/%2C/,",",pr); gsub(/%3B/,";",pr);
  print lt"\t"$4"\t"$5"\t"$7"\t"pr"\t"pid
}' "$GFF" | sort -u > v2/cds_full.tsv

cut -f1-5 v2/cds_full.tsv > v2/cds_todos.tsv
awk -F'\t' '$6!=""{print $1"\t"$6}' v2/cds_full.tsv | sort -u > v2/lt2prot.tsv
grep -Ei 'flagell|chemotax|hook|motility|basal.?body' v2/cds_todos.tsv | sort -u > v2/flag_gff.tsv
cut -f2 v2/lt2prot.tsv | sort | uniq -d > v2/prot_dup.txt

echo "CDS totales            : $(wc -l < v2/cds_todos.tsv)"
echo "Con protein_id         : $(wc -l < v2/lt2prot.tsv)"
echo "Accesiones unicas      : $(cut -f2 v2/lt2prot.tsv | sort -u | wc -l)"
echo "Accesiones duplicadas  : $(wc -l < v2/prot_dup.txt)"
echo "Genes flagelares       : $(wc -l < v2/flag_gff.tsv)"
echo
echo "--- conjunto flagelar ---"
cut -f1,5 v2/flag_gff.tsv
echo
echo "--- comparacion con lo anterior ---"
for f in cds_todos.tsv lt2prot.tsv flag_gff.tsv; do
  if [ -f "$f" ]; then
    diff -q "$f" "v2/$f" >/dev/null \
      && echo "IGUAL    $f" \
      || echo "DISTINTO $f   viejo=$(wc -l < $f)  nuevo=$(wc -l < v2/$f)"
  else echo "AUSENTE  $f"; fi
done
