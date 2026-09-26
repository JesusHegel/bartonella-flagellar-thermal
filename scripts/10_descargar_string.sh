#!/bin/bash
# 10. Descarga de la base STRING para KC583 (especie 360095).
#
# Solo descarga si datos/string/ no existe todavia. La copia descargada se
# guarda en el repositorio con su fecha y su version, de modo que cualquier
# replica usa exactamente los mismos datos, sin internet y aunque STRING cambie.
# Para volver a descargar a proposito:  FORZAR_DESCARGA_STRING=1 bash scripts/10_descargar_string.sh
#
# Se usan los locus antiguos (BARBAKC583_xxxx), que son los identificadores de
# STRING para esta cepa. Se descarga, para todas las proteinas del genoma:
#   - la correspondencia locus -> identificador STRING
#   - la red de asociaciones con puntaje combinado >= 0,4 (confianza media)
#   - la anotacion funcional de cada proteina (KEGG, GO, palabras clave, etc.)
# El analisis se hace despues, en local (script 11), con nuestro universo.
set -uo pipefail
source scripts/rutas.sh
D=datos/string
ESPECIE=360095
QUIEN=bartonella_unmsm_reanalisis

if [ -f "$D/version.tsv" ] && [ "${FORZAR_DESCARGA_STRING:-0}" != 1 ]; then
  echo "STRING ya descargado el $(cat "$D/fecha_descarga.txt"); se usa la copia guardada:"
  cat "$D/version.tsv"
  exit 0
fi
command -v curl >/dev/null || { echo "ERROR: falta curl"; exit 1; }
mkdir -p "$D"

# Version y direccion estable (la que no cambia aunque salga una version nueva)
curl -sS --fail "https://string-db.org/api/tsv/version" > "$D/version.tsv" || { echo "ERROR: no responde STRING"; exit 1; }
BASE=$(awk -F'\t' 'NR==2{print $2}' "$D/version.tsv")
[ -n "$BASE" ] || { echo "ERROR: no se pudo leer la direccion estable"; cat "$D/version.tsv"; exit 1; }
echo "STRING: $(awk -F'\t' 'NR==2{print $1}' "$D/version.tsv")  ->  $BASE"

# Locus antiguos de todas las CDS, separados por retorno de carro (formato de STRING)
Rscript -e 'g <- read.csv("resultados/genes/tabla_de_genes.csv", stringsAsFactors = FALSE)
            x <- sort(unique(na.omit(g$locus_antiguo))); cat(paste(x, collapse = "\r"))' > "$DIR_INT/string_consulta.txt"
echo "Locus consultados: $(tr '\r' '\n' < "$DIR_INT/string_consulta.txt" | wc -l)"

consulta() {   # consulta <metodo> <archivo de salida> [parametros extra...]
  local metodo=$1 salida=$2; shift 2
  curl -sS --fail -X POST "$BASE/api/tsv/$metodo" \
       --data-urlencode "identifiers@$IDS" -d "species=$ESPECIE" -d "caller_identity=$QUIEN" "$@" > "$salida" \
    || { echo "ERROR en $metodo"; rm -f "$salida"; exit 1; }
  echo "  $metodo: $(($(wc -l < "$salida") - 1)) filas"
}

IDS=$DIR_INT/string_consulta.txt
consulta get_string_ids "$D/identificadores.tsv" -d "limit=1" -d "echo_query=1"
sleep 2

# Para la red y la anotacion se usan ya los identificadores de STRING
Rscript -e 'm <- read.delim("datos/string/identificadores.tsv", quote = "", stringsAsFactors = FALSE)
            cat(paste(sort(unique(m$stringId)), collapse = "\r"))' > "$DIR_INT/string_ids.txt"
IDS=$DIR_INT/string_ids.txt
consulta network "$D/red.tsv" -d "required_score=400" -d "network_type=functional"
sleep 2
consulta functional_annotation "$D/anotacion_funcional.tsv" -d "allow_pubmed=0"

date -u +%F > "$D/fecha_descarga.txt"
( cd "$D" && md5sum identificadores.tsv red.tsv anotacion_funcional.tsv version.tsv > md5.txt )
echo "Descarga guardada en $D/ ($(du -sh "$D" | cut -f1))"
