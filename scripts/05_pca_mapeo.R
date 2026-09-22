setwd("~/bbrna")
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
