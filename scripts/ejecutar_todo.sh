#!/bin/bash
# Ejecuta el analisis completo, del script 01 al 12, y al final lo verifica.
#
# Uso, desde la carpeta raiz del repositorio y con el entorno activado:
#   conda activate bartonella
#   bash scripts/ejecutar_todo.sh
#
# Antes de empezar borra resultados/, figuras/ e intermedios/, para que
# ningun archivo viejo pueda hacerse pasar por uno recien generado.
# (Si hace falta recuperarlos: git checkout -- resultados figuras registros)
set -uo pipefail
source scripts/rutas.sh

echo "== Comprobando programas =="
falta=0
for p in Rscript python3 blastp makeblastdb datasets RNAfold unzip md5sum curl; do
  command -v $p >/dev/null 2>&1 && echo "  ok      $p" || { echo "  FALTA   $p"; falta=1; }
done
[ $falta -eq 0 ] || { echo "Faltan programas. ¿Activaste el entorno? (conda activate bartonella)"; exit 1; }

echo "== Borrando salidas anteriores =="
rm -rf "${DIR_RES:?}" "${DIR_FIG:?}" "${DIR_INT:?}" "${DIR_REG:?}"
source scripts/rutas.sh          # vuelve a crear las carpetas vacias

correr() {                        # correr <script>  -> registros/<script>.log
  local s=$1 log=$DIR_REG/${1%.*}.log
  echo "== $s  ($(date +%H:%M:%S))"
  case $s in
    *.sh) bash scripts/$s    > "$log" 2>&1 ;;
    *.R)  Rscript scripts/$s > "$log" 2>&1 ;;
  esac
  local e=$?
  [ $e -eq 0 ] || { echo "ERROR en $s (codigo $e). Revisa $log"; exit 1; }
}
correr 01_anotacion.sh
correr 02_expresion_diferencial.R
correr 03_enriquecimiento.R
correr 04_robustez.R
correr 05_pca_y_tasa_mapeo.R
correr 06_panel_genero_y_ortologia.sh
correr 07_quimiotaxis_y_utr.sh
correr 08_tabla_de_genes.R
correr 09_tabla_de_muestras.R
correr 10_descargar_string.sh
correr 11_analisis_string.R
correr 12_figuras.R

Rscript -e 'suppressMessages({library(tximport); library(DESeq2); library(ashr); library(fgsea); library(ggplot2)}); sessionInfo()' > "$DIR_REG/versiones_R.txt" 2>&1
{ salmon --version; prefetch --version; blastp -version | head -1; RNAfold --version; datasets --version; python3 --version; } \
  > "$DIR_REG/versiones_programas.txt" 2>&1
echo "== Terminado ($(date +%H:%M:%S)) =="
echo
bash scripts/verificar.sh
