# Carpetas del proyecto. Todos los scripts de R cargan este archivo.
# Por defecto se trabaja con las cuantificaciones incluidas en el repositorio.
# Para la replica desde el SRA se cambian con variables de entorno (ver README).

if (!file.exists("CIFRAS_CONFIRMADAS.md"))
  stop("Ejecuta los scripts desde la carpeta raiz del repositorio.")

leer_variable <- function(nombre, por_defecto) {
  valor <- Sys.getenv(nombre)
  if (nzchar(valor)) valor else por_defecto
}
DIR_CUANT <- leer_variable("DIR_CUANT", "cuantificacion")
DIR_RES   <- leer_variable("DIR_RES",   "resultados")
DIR_INT   <- leer_variable("DIR_INT",   "intermedios")
DIR_FIG   <- leer_variable("DIR_FIG",   "figuras")

for (d in c(file.path(DIR_RES, c("anotacion", "expresion_diferencial", "enriquecimiento",
                                 "robustez", "muestras", "genomica_comparada")),
            DIR_INT, DIR_FIG))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Atajos para las rutas mas usadas
res <- function(...) file.path(DIR_RES, ...)
int <- function(...) file.path(DIR_INT, ...)
fig <- function(...) file.path(DIR_FIG, ...)
