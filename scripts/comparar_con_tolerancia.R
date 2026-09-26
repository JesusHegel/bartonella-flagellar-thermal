# Compara una tabla regenerada con la publicada cuando no son identicas byte a byte.
#
# Uso: Rscript scripts/comparar_con_tolerancia.R <tabla_nueva> <tabla_publicada>
#
# Algunas tablas (contrastes de DESeq2, fgsea, PCA) dependen de calculos
# numericos cuyo ultimo decimal varia segun el procesador y la libreria de
# algebra lineal. Se aceptan como equivalentes si se cumplen TODAS estas
# condiciones:
#   1. mismas filas, mismas columnas y mismos textos
#   2. valores faltantes (NA) en las mismas celdas
#   3. cada valor numerico difiere en menos de 1e-4 (absoluto) o menos de
#      1 % (relativo)
#   4. en los contrastes: ningun gen cambia de significativo (padj < 0,05) a
#      no significativo ni al reves, y ningun log2FoldChange mayor que 1e-4
#      en valor absoluto cambia de signo
# Referencia: entre dos instalaciones distintas en la misma computadora se
# observaron diferencias maximas de 1,4e-5 en log2FoldChange y 3,9e-5 en
# pvalue, sin ningun cambio de significancia (CIFRAS_CONFIRMADAS.md, 12c).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Uso: Rscript scripts/comparar_con_tolerancia.R <nueva> <publicada>")
x <- read.csv(args[1], check.names = FALSE); y <- read.csv(args[2], check.names = FALSE)

falla <- function(motivo) { cat("NO EQUIVALENTE:", motivo, "\n"); quit(status = 1) }
if (!identical(dim(x), dim(y)) || !identical(names(x), names(y))) falla("filas o columnas distintas")
num <- vapply(y, is.numeric, logical(1))
if (!identical(vapply(x, is.numeric, logical(1)), num)) falla("tipos de columna distintos")
if (!identical(x[!num], y[!num])) falla("textos distintos")
X <- as.matrix(x[num]); Y <- as.matrix(y[num])
if (!identical(is.na(X), is.na(Y))) falla("valores faltantes en celdas distintas")

d  <- abs(X - Y); r <- d / pmax(abs(Y), 1e-300)
ok <- is.na(d) | d < 1e-4 | r < 0.01
if (!all(ok)) falla(sprintf("%d valores fuera de tolerancia", sum(!ok)))

detalle <- ""
if (all(c("log2FoldChange", "padj") %in% names(y))) {
  sx <- !is.na(x$padj) & x$padj < 0.05; sy <- !is.na(y$padj) & y$padj < 0.05
  if (any(sx != sy)) falla(sprintf("%d genes cambian de significancia", sum(sx != sy)))
  # El signo solo cuenta cuando el cambio no es practicamente cero (ashr deja
  # muchos genes en valores del orden de 1e-7, cuyo signo es ruido).
  claro <- !is.na(y$log2FoldChange) & abs(y$log2FoldChange) >= 1e-4
  if (any(sign(x$log2FoldChange[claro]) != sign(y$log2FoldChange[claro]), na.rm = TRUE))
    falla("algun log2FoldChange cambia de signo")
  detalle <- sprintf(", %d significativos iguales", sum(sy))
}
cat(sprintf("equivalente (dif. absoluta max %.1e, relativa max %.1e%s)\n",
            max(d, na.rm = TRUE), max(r, na.rm = TRUE), detalle))
