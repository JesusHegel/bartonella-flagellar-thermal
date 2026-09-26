#!/bin/bash
# Compara cada archivo con la huella md5 registrada en verificacion/sumas_md5.txt.
#
#   bash scripts/verificar.sh             verifica y resume
#   bash scripts/verificar.sh --generar   rehace la lista (solo el autor, tras una
#                                          corrida validada)
#
# La lista cubre los datos de entrada (datos/, cuantificacion/) y todas las
# tablas de resultados/. Las figuras no se comparan por md5 porque los PNG y
# PDF guardan la fecha de creacion dentro del archivo: se revisan a la vista.
set -uo pipefail
source scripts/rutas.sh
LISTA=verificacion/sumas_md5.txt

if [ "${1:-}" = "--generar" ]; then
  mkdir -p verificacion
  find datos cuantificacion resultados -type f | sort | xargs md5sum > "$LISTA"
  echo "Lista rehecha: $(wc -l < "$LISTA") archivos en $LISTA"
  exit 0
fi

[ -f "$LISTA" ] || { echo "ERROR: falta $LISTA"; exit 1; }
iguales=0; distintos=(); faltan=()
while read -r suma ruta; do
  case $ruta in
    resultados/*)     real="$DIR_RES/${ruta#resultados/}" ;;
    cuantificacion/*) [ "$DIR_CUANT" = cuantificacion ] || continue; real=$ruta ;;
    *)                real=$ruta ;;
  esac
  if [ ! -f "$real" ]; then faltan+=("$real")
  elif [ "$(md5sum < "$real" | cut -d' ' -f1)" = "$suma" ]; then iguales=$((iguales+1))
  else distintos+=("$real"); fi
done < "$LISTA"
nuevos=$(comm -13 <(grep ' resultados/' "$LISTA" | sed "s#  resultados/#  $DIR_RES/#" | cut -d' ' -f3- | sort) \
                  <(find "$DIR_RES" -type f | sort))

echo "===== VERIFICACION ($DIR_RES) ====="
echo "Identicos : $iguales"
echo "Distintos : ${#distintos[@]}";  for f in "${distintos[@]}"; do echo "   DISTINTO  $f"; done
echo "Faltan    : ${#faltan[@]}";     for f in "${faltan[@]}";    do echo "   FALTA     $f"; done
[ -n "$nuevos" ] && { echo "Archivos nuevos, sin huella registrada:"; echo "$nuevos" | sed 's/^/   NUEVO     /'; }
if [ ${#distintos[@]} -eq 0 ] && [ ${#faltan[@]} -eq 0 ]; then
  echo "RESULTADO: todo identico."
else
  echo "RESULTADO: hay diferencias. No sigas hasta entender la causa."; exit 1
fi
