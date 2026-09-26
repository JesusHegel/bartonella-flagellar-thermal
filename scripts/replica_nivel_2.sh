#!/bin/bash
# Replica completa (nivel 2): desde las lecturas crudas del SRA hasta la
# comparacion con los resultados publicados.
#
# Uso, desde la carpeta raiz del repositorio y con el entorno activado:
#   conda activate bartonella
#   bash scripts/replica_nivel_2.sh
#
# 1. Descarga y cuantifica las 23 corridas (script 00) en cuantificacion_desde_sra/.
#    Si se interrumpe, al volver a correrlo retoma donde quedo.
# 2. Corre el analisis completo (01 a 12) sobre esas cuantificaciones, en
#    carpetas aparte que terminan en _desde_sra (nunca toca las publicadas).
# 3. Compara las cifras principales con las publicadas (comparar_nivel_2.R).
set -uo pipefail
source scripts/rutas.sh

bash scripts/00_descargar_y_cuantificar.sh || exit 1
n=$(ls cuantificacion_desde_sra/*/quant.sf 2>/dev/null | wc -l)
[ "$n" -eq 23 ] || { echo "ERROR: hay $n de 23 corridas cuantificadas. Vuelve a correr este script: retoma donde quedo."; exit 1; }

echo; echo "== Analisis completo sobre las cuantificaciones desde el SRA =="
DIR_CUANT=cuantificacion_desde_sra DIR_RES=resultados_desde_sra DIR_INT=intermedios_desde_sra \
DIR_FIG=figuras_desde_sra DIR_REG=registros_desde_sra VERIFICAR=0 \
  bash scripts/ejecutar_todo.sh || exit 1

echo; echo "== Comparacion con los resultados publicados =="
Rscript scripts/comparar_nivel_2.R 2>&1 | tee registros_desde_sra/comparar_nivel_2.log
