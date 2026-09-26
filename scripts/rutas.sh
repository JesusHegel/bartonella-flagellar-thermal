# Carpetas del proyecto. Todos los scripts de bash cargan este archivo.
# Por defecto se trabaja con las cuantificaciones incluidas en el repositorio.
# Para la replica desde el SRA se cambian con variables de entorno (ver README).

[ -f CIFRAS_CONFIRMADAS.md ] || { echo "ERROR: ejecuta los scripts desde la carpeta raiz del repositorio."; exit 1; }

export LC_ALL=C                              # mismo orden de 'sort' en cualquier computadora
DIR_CUANT=${DIR_CUANT:-cuantificacion}       # cuantificaciones de salmon (entrada)
DIR_RES=${DIR_RES:-resultados}               # tablas de resultados
DIR_INT=${DIR_INT:-intermedios}              # objetos intermedios (no se suben a git)
DIR_FIG=${DIR_FIG:-figuras}                  # figuras
DIR_DESC=${DIR_DESC:-descargas}              # proteomas descargados del NCBI (no se suben a git)
DIR_REG=${DIR_REG:-registros}                # salida de pantalla de cada script
REF=datos/referencia                         # genoma de referencia KC583 (GCF_000015445.1)

mkdir -p "$DIR_RES"/{anotacion,expresion_diferencial,enriquecimiento,robustez,muestras,genomica_comparada} \
         "$DIR_INT" "$DIR_FIG" "$DIR_DESC" "$DIR_REG"
