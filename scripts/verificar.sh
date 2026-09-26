#!/bin/bash
# Compara cada tabla regenerada con la version publicada en el repositorio (git).
#
#   bash scripts/verificar.sh
#
# - La mayoria de tablas deben ser IDENTICAS byte a byte.
# - Los contrastes de DESeq2, fgsea y la PCA pueden variar en el ultimo decimal
#   segun el procesador; si no son identicas se comparan con tolerancia
#   (scripts/comparar_con_tolerancia.R) y se aceptan como EQUIVALENTES solo si
#   no cambia ninguna decision de significancia.
# - Las figuras no se comparan: sus archivos guardan la fecha de creacion.
#   Se revisan a la vista.
set -uo pipefail
source scripts/rutas.sh
git rev-parse --git-dir >/dev/null 2>&1 || { echo "ERROR: esta carpeta no es un repositorio git"; exit 1; }

con_tolerancia() {
  case $1 in
    resultados/expresion_diferencial/*|resultados/enriquecimiento/enriquecimiento_fgsea.csv|resultados/muestras/pca_y_tasa_mapeo.csv|\
    resultados/genes/*|resultados/string/red_cambio_2x_nodos.csv) return 0 ;;
    *) return 1 ;;
  esac
}

echo "===== VERIFICACION ($DIR_RES frente a la version publicada) ====="
git diff --quiet HEAD -- datos cuantificacion \
  && echo "Datos de entrada: sin cambios" \
  || echo "AVISO: los datos de entrada (datos/ o cuantificacion/) no coinciden con los publicados"

TMP=$(mktemp)
iguales=0; equivalentes=(); distintos=(); faltan=()
for ruta in $(git ls-files resultados); do
  m=""
  real="$DIR_RES/${ruta#resultados/}"
  if [ ! -f "$real" ]; then faltan+=("$real"); continue; fi
  git show "HEAD:$ruta" > "$TMP"
  if cmp -s "$real" "$TMP"; then
    iguales=$((iguales+1))
  elif con_tolerancia "$ruta" && m=$(Rscript scripts/comparar_con_tolerancia.R "$real" "$TMP" 2>&1); then
    equivalentes+=("$real  ->  $m")
  else
    distintos+=("$real  $m")
  fi
done
rm -f "$TMP"
nuevos=$(comm -13 <(git ls-files resultados | sed "s#^resultados/#$DIR_RES/#" | sort) <(find "$DIR_RES" -type f | sort))

echo "Identicos    : $iguales"
echo "Equivalentes : ${#equivalentes[@]}"; for f in "${equivalentes[@]}"; do echo "   $f"; done
echo "Distintos    : ${#distintos[@]}";    for f in "${distintos[@]}";    do echo "   DISTINTO  $f"; done
echo "Faltan       : ${#faltan[@]}";       for f in "${faltan[@]}";       do echo "   FALTA     $f"; done
[ -n "$nuevos" ] && { echo "Archivos nuevos, sin version publicada:"; echo "$nuevos" | sed 's/^/   NUEVO     /'; }

if [ ${#distintos[@]} -gt 0 ] || [ ${#faltan[@]} -gt 0 ]; then
  echo "RESULTADO: hay diferencias. No sigas hasta entender la causa."; exit 1
elif [ ${#equivalentes[@]} -gt 0 ]; then
  echo "RESULTADO: correcto. $iguales tablas identicas y ${#equivalentes[@]} equivalentes dentro de la tolerancia."
else
  echo "RESULTADO: todo identico."
fi
