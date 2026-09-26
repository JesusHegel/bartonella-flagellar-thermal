# Carpetas del proyecto. Todos los scripts de bash cargan este archivo.
# Por defecto se trabaja con las cuantificaciones incluidas en el repositorio.
# Para la replica desde el SRA se cambian con variables de entorno (ver README).

[ -f CIFRAS_CONFIRMADAS.md ] || { echo "ERROR: ejecuta los scripts desde la carpeta raiz del repositorio."; exit 1; }

# Idioma: texto en UTF-8 (acentos correctos en las figuras) y orden alfabetico
# estricto "C" (mayusculas antes que minusculas), igual en cualquier computadora.
# El orden importa: R ordena con el las condiciones de DESeq2 y las listas de genes.
unset LC_ALL
export LANG=C.UTF-8 LC_CTYPE=C.UTF-8 LC_COLLATE=C
DIR_CUANT=${DIR_CUANT:-cuantificacion}       # cuantificaciones de salmon (entrada)
DIR_RES=${DIR_RES:-resultados}               # tablas de resultados
DIR_INT=${DIR_INT:-intermedios}              # objetos intermedios (no se suben a git)
DIR_FIG=${DIR_FIG:-figuras}                  # figuras
DIR_DESC=${DIR_DESC:-descargas}              # proteomas descargados del NCBI (no se suben a git)
DIR_REG=${DIR_REG:-registros}                # salida de pantalla de cada script
REF=datos/referencia                         # genoma de referencia KC583 (GCF_000015445.1)

mkdir -p "$DIR_RES"/{anotacion,expresion_diferencial,enriquecimiento,robustez,muestras,genomica_comparada,genes,string} \
         "$DIR_INT" "$DIR_FIG" "$DIR_DESC" "$DIR_REG"
