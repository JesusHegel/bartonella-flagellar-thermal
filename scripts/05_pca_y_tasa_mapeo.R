# 05. Relacion entre la PCA y la tasa de asignacion de lecturas.
#
# La tasa de cada muestra se toma de los informes de salmon
# (aux_info/meta_info.json, campo percent_mapped). Como el indice es un
# gentrome, esa tasa es la fraccion de lecturas asignadas a CDS bacterianas,
# no el contenido bacteriano total.
source("scripts/rutas.R")

# Se conserva el texto original del JSON para no perder decimales.
js <- Sys.glob(file.path(DIR_CUANT, "*", "aux_info", "meta_info.json"))
if (!length(js)) stop("No encuentro ", DIR_CUANT, "/*/aux_info/meta_info.json")
mp <- data.frame(
  run   = basename(dirname(dirname(js))),
  mapeo = vapply(js, function(f) {
    l <- grep('"percent_mapped"', readLines(f, warn = FALSE), value = TRUE)
    sub('.*:[[:space:]]*([0-9.eE+-]+).*', '\\1', l[1])
  }, character(1)), stringsAsFactors = FALSE)
write.table(mp, res("muestras", "tasa_mapeo.tsv"), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("tasa_mapeo.tsv generado desde", length(js), "informes de salmon\n")

p <- read.csv(int("pca.csv"), row.names = 1); p$run <- rownames(p)
m <- read.delim(res("muestras", "tasa_mapeo.tsv"), header = FALSE); names(m) <- c("run", "mapeo")
d <- merge(p, m, by = "run"); d <- d[order(d$PC1), ]
d$matriz <- ifelse(d$cond %in% c("HB37", "HBBG"), "sangre humana",
            ifelse(d$cond == "HUVE", "celulas endoteliales",
            ifelse(d$cond == "PlBG", "placa + atmosfera sanguinea", "placa / caldo")))
print(d[, c("run", "cond", "matriz", "mapeo", "PC1", "PC2")])
cat("\n*** Correlacion PC1 ~ tasa de mapeo:", round(cor(d$PC1, d$mapeo), 3), "***\n")
cat("*** Correlacion PC2 ~ tasa de mapeo:", round(cor(d$PC2, d$mapeo), 3), "***\n\n")
print(aggregate(cbind(PC1, mapeo) ~ matriz, data = d, FUN = function(x) round(mean(x), 1)))
write.csv(d, res("muestras", "pca_y_tasa_mapeo.csv"), row.names = FALSE)
