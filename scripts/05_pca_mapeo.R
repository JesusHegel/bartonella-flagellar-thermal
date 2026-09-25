setwd("~/bbrna")

# La tasa de mapeo de cada muestra se extrae de los informes de salmon.
# Se conserva el texto original del JSON para no perder decimales.
if (!file.exists("salida/mapeo.tsv")) {
  js <- Sys.glob("quants/*/aux_info/meta_info.json")
  if (!length(js)) stop("No encuentro quants/*/aux_info/meta_info.json")
  mp <- data.frame(
    run   = basename(dirname(dirname(js))),
    mapeo = vapply(js, function(f) {
      l <- grep('"percent_mapped"', readLines(f, warn=FALSE), value=TRUE)
      sub('.*:[[:space:]]*([0-9.eE+-]+).*', '\\1', l[1])
    }, character(1)), stringsAsFactors=FALSE)
  write.table(mp, "salida/mapeo.tsv", sep="\t",
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  cat("mapeo.tsv generado desde", length(js), "informes de salmon\n")
}

p <- read.csv("salida/pca.csv", row.names=1); p$run <- rownames(p)
m <- read.delim("salida/mapeo.tsv", header=FALSE); names(m) <- c("run","mapeo")
d <- merge(p, m, by="run"); d <- d[order(d$PC1),]
d$matriz <- ifelse(d$cond %in% c("HB37","HBBG"),"sangre humana",
            ifelse(d$cond=="HUVE","celulas endoteliales",
            ifelse(d$cond=="PlBG","placa + atmosfera sanguinea","placa / caldo")))
print(d[,c("run","cond","matriz","mapeo","PC1","PC2")])
cat("\n*** Correlacion PC1 ~ tasa de mapeo:", round(cor(d$PC1,d$mapeo),3), "***\n")
cat("*** Correlacion PC2 ~ tasa de mapeo:", round(cor(d$PC2,d$mapeo),3), "***\n\n")
print(aggregate(cbind(PC1,mapeo) ~ matriz, data=d, FUN=function(x) round(mean(x),1)))
write.csv(d, "salida/chk_pca_mapeo.csv", row.names=FALSE)
