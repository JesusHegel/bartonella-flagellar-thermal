#!/bin/bash
# 00. Descarga de las 23 corridas del SRA, control de calidad y cuantificacion.
#
# Solo hace falta para la replica completa (nivel 2). El analisis normal parte
# de las cuantificaciones ya incluidas en cuantificacion/.
#
# Escribe en cuantificacion_desde_sra/ (nunca toca cuantificacion/).
# Necesita internet y unos 20 GB libres: cada corrida se descarga, se
# cuantifica y se borra antes de pasar a la siguiente.
# Tiempo de referencia: unos 50 minutos para las 23 corridas.
#
# El informe de calidad (fastqc + multiqc) se escribe en control_de_calidad/
# y es opcional: ver entorno/control_de_calidad.yml. Se hace sobre las
# lecturas sin recortar, que es como entran a salmon.
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

# ---- 2. Control de calidad de las lecturas (opcional) ----
# fastqc y multiqc van en un entorno aparte (entorno/control_de_calidad.yml)
# para no tocar el entorno del analisis, que ya esta verificado. Si no estan
# instalados el script avisa y sigue: ningun resultado depende de ellos.
# La salida va fuera de resultados/, asi que verificar.sh no la ve.
DIR_QC=${DIR_QC:-control_de_calidad}
FASTQC=$(command -v fastqc || true)
if [ -n "$FASTQC" ]; then
  mkdir -p "$DIR_QC"
  echo "Control de calidad: activado ($(fastqc --version 2>/dev/null | grep -o 'FastQC.*'))"
else
  echo "AVISO: fastqc no esta instalado; no se generara el informe de calidad."
  echo "       Para generarlo: conda env create -f entorno/control_de_calidad.yml"
fi

# ---- 3. Descarga y cuantificacion, una corrida a la vez ----
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
    # Sobre las lecturas sin recortar, antes de que se borren mas abajo.
    # -t 2 porque son dos archivos por corrida (fastqc pide ~250 MB por hilo).
    if [ -n "$FASTQC" ]; then
      echo "[$(date +%H:%M)] $R control de calidad"
      fastqc -t 2 -o "$DIR_QC" --quiet "$TMP/${R}_1.fastq" "$TMP/${R}_2.fastq" \
        >> "$TMP/$R.log" 2>&1 || echo "  FALLO fastqc $R"
    fi
  else
    echo "  SIN PAREJA $R"
  fi
  rm -rf "${TMP:?}/$R" "$TMP"/${R}*.fastq
done
N_QUANT=$(ls "$SALIDA"/*/quant.sf 2>/dev/null | wc -l)
echo "=== TERMINADO: $N_QUANT de 23 en $SALIDA ==="

# ---- 4. Informe combinado de calidad ----
if [ -n "$FASTQC" ]; then
  N_QC=$(ls "$DIR_QC"/*_fastqc.html 2>/dev/null | wc -l)
  echo "=== Calidad: $N_QC informes en $DIR_QC (2 por corrida, se esperan $((N_QUANT * 2))) ==="
  if [ "$N_QC" -ne "$((N_QUANT * 2))" ]; then
    echo "AVISO: faltan informes de calidad. Las corridas que ya estaban cuantificadas"
    echo "       se saltan enteras. Para generarlos todos: borra $SALIDA y vuelve a correr."
  fi
  if command -v multiqc > /dev/null 2>&1; then
    if multiqc "$DIR_QC" -o "$DIR_QC" -n informe_multiqc --force --quiet > "$DIR_QC/multiqc.log" 2>&1; then
      echo "Informe combinado: $DIR_QC/informe_multiqc.html"
      # Tabla por muestra, que es la que va al material suplementario.
      # Se copia a un nombre fijo porque el resto de informe_multiqc_data/
      # no se sube (ver .gitignore).
      if [ -f "$DIR_QC/informe_multiqc_data/multiqc_fastqc.txt" ]; then
        cp "$DIR_QC/informe_multiqc_data/multiqc_fastqc.txt" "$DIR_QC/resumen_calidad.tsv"
        echo "Resumen por muestra : $DIR_QC/resumen_calidad.tsv"
      fi
    else
      echo "  FALLO multiqc (quedan los informes individuales de fastqc)"
    fi
  else
    echo "AVISO: multiqc no esta instalado; quedan los $N_QC informes individuales."
  fi
fi
