# Categorias funcionales usadas en el enriquecimiento (script 03) y en la tabla
# de genes (script 08). Se definen por palabras clave sobre el producto anotado,
# salvo la flagelar, que es el conjunto fijado a priori por el script 01.
# Estar en un solo archivo garantiza que ambos scripts usen exactamente la misma
# definicion.
PATRONES_CATEGORIAS <- c(
  Ribosoma   = "ribosomal protein",
  Chaperona  = "chaperon|heat shock|GroE|DnaJ|DnaK|ClpB|HtpX",
  TranspABC  = "ABC transporter",
  HierroHemo = "hemin|heme|iron|TonB|ferr",
  Traduccion = "tRNA|elongation factor|translation",
  Hipotetica = "hypothetical")

# Devuelve una lista: para cada categoria, los locus que pertenecen a ella.
definir_categorias <- function(locus, producto, flagelares) {
  c(list(Flagelar = locus[locus %in% flagelares]),
    lapply(PATRONES_CATEGORIAS, function(p) locus[grepl(p, producto, ignore.case = TRUE)]))
}
