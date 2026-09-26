#!/bin/bash
# 00. Descarga de las 23 corridas del SRA y cuantificacion con salmon.
#
# Solo hace falta para la replica completa (nivel 2). El analisis normal parte
# de las cuantificaciones ya incluidas en cuantificacion/.
#
# Escribe en cuantificacion_desde_sra/ (nunca toca cuantificacion/).
# Necesita internet y unos 20 GB libres: cada corrida se descarga, se
# cuantifica y se borra antes de pasar a la siguiente.
# Tiempo de referencia: unos 50 minutos para las 23 corridas.
set -uo pipefail
source scripts/rutas.sh
SALIDA=${DIR_CUANT_SRA:-cuantificacion_desde_sra}
HILOS=${HILOS:-6}
TMP=$DIR_INT/sra_temporal
mkdir -p "$SALIDA" "$TMP"

# ---- 1. Indice decoy-aware (gentrome = CDS + genoma; senuelo = cromosoma) ----
IDX=$DIR_INT/indice_salmon
if [ ! -f "$IDX/info.json" ]; then
  cat "$REF/cds_from_genomic.fna" "$REF/GCF_000015445.1_ASM1544v1_genomic.fna" > "$DIR_INT/gentrome.fna"
  grep '^>' "$REF/GCF_000015445.1_ASM1544v1_genomic.fna" | cut -d' ' -f1 | tr -d '>' > "$DIR_INT/senuelo.txt"
  echo "gentrome md5: $(md5sum < "$DIR_INT/gentrome.fna" | cut -d' ' -f1)  (esperado bc9bf0e33d39363355d2e013eaf31e44)"
  salmon index -t "$DIR_INT/gentrome.fna" -d "$DIR_INT/senuelo.txt" -i "$IDX" -k 31 -p "$HILOS" > "$DIR_INT/indice_salmon.log" 2>&1
fi
# El indice debe coincidir con el original en sus seis huellas
declare -A ESPERADO=(
  [seq_hash]=12753c8cd9dedf605a0f039c3c962431626522909d3a96f22d5dab18e723464f
  [name_hash]=41f23fc4821c3f6399f8fe37baf0d789fb0c4111766c6de7fb6d294621af9875
  [decoy_seq_hash]=d55bc36f256de6ffcf09122e72f0c0899e016c99a834e1c2104357b906310e5f
  [decoy_name_hash]=007b37bd5bd10390ec03620647723773955d1896e6bc004c2edc8fe6cfe0d9fb
  [num_refs]=1190 [first_decoy_index]=1189 )
for k in "${!ESPERADO[@]}"; do
  v=$(grep "\"$k\"" "$IDX/info.json" | sed 's/.*: *//; s/[",]//g')
  [ "$v" = "${ESPERADO[$k]}" ] || { echo "ERROR: el indice no coincide con el original ($k)"; exit 1; }
done
echo "Indice de salmon: las seis huellas coinciden con el original"

# ---- 2. Descarga y cuantificacion, una corrida a la vez ----
for R in $(tail -n +2 datos/muestras.tsv | cut -f1); do
  if [ -f "$SALIDA/$R/quant.sf" ]; then echo "[ya estaba] $R"; continue; fi
  echo "[$(date +%H:%M)] $R descargando"
  prefetch $R -O "$TMP" --max-size 50G > "$TMP/$R.log" 2>&1 || { echo "  FALLO prefetch $R"; continue; }
  fasterq-dump "$TMP/$R/$R.sra" -O "$TMP" -t "$TMP" -e "$HILOS" --split-3 >> "$TMP/$R.log" 2>&1 \
    || { echo "  FALLO fasterq-dump $R"; rm -rf "${TMP:?}/$R"; continue; }
  if [ -f "$TMP/${R}_1.fastq" ] && [ -f "$TMP/${R}_2.fastq" ]; then
    echo "[$(date +%H:%M)] $R cuantificando"
    salmon quant -i "$IDX" -l A -1 "$TMP/${R}_1.fastq" -2 "$TMP/${R}_2.fastq" -p "$HILOS" \
      -o "$SALIDA/$R" >> "$TMP/$R.log" 2>&1 || echo "  FALLO salmon $R"
  else
    echo "  SIN PAREJA $R"
  fi
  rm -rf "${TMP:?}/$R" "$TMP"/${R}*.fastq
done
echo "=== TERMINADO: $(ls "$SALIDA"/*/quant.sf 2>/dev/null | wc -l) de 23 en $SALIDA ==="
