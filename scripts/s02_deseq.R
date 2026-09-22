suppressMessages({library(tximport); library(DESeq2); library(ashr)})
setwd("~/bbrna")

dirs <- list.dirs("quants", recursive=FALSE, full.names=FALSE)
cat("Carpetas en quants/:", length(dirs), "\n")

sm <- read.delim("samples.tsv", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
cat("samples.tsv:", nrow(sm), "filas |", paste(names(sm), collapse=" | "), "\n")

ci <- which(sapply(sm, function(x) mean(dirs %in% as.character(x)) > 0.8))[1]
cc <- which(sapply(sm, function(x) any(grepl("^Pl(25|30|37)$", as.character(x)))))[1]
if (is.na(ci) || is.na(cc)) { print(head(sm)); stop("No detecto las columnas; pegame esta salida.") }
cat("Columna muestra:", names(sm)[ci], " | Columna condicion:", names(sm)[cc], "\n")

sm    <- sm[match(dirs, as.character(sm[[ci]])), ]
files <- file.path("quants", dirs, "quant.sf"); names(files) <- dirs
stopifnot(all(file.exists(files)))

txi <- tximport(files, type="salmon", txOut=TRUE, dropInfReps=TRUE)
cd  <- data.frame(row.names=dirs, cond=factor(as.character(sm[[cc]])))
cat("\nMuestras por condicion:\n"); print(table(cd$cond))

dds <- DESeqDataSetFromTximport(txi, colData=cd, design=~cond)
cat("\nTranscritos totales    :", nrow(dds), "\n")
dds <- dds[rowSums(counts(dds)) > 0, ]
cat("Transcritos analizables:", nrow(dds), "\n")
dds <- DESeq(dds, quiet=TRUE); saveRDS(dds, "v2/dds_v2.rds")

vsd <- vst(dds, blind=TRUE)
p   <- plotPCA(vsd, intgroup="cond", ntop=500, returnData=TRUE)
pv  <- round(100*attr(p,"percentVar"))
cat("\n*** PC1 =", pv[1], "%   PC2 =", pv[2], "% ***\n")
write.csv(p, "v2/pca_v2.csv")
print(p[order(p$PC1), c("cond","PC1","PC2")])

yale <- cd$cond %in% c("Pl25","Pl30","Pl37","pH06","pH07","pH08")
cat("\nMuestras en el lote restringido:", sum(yale), "\n")
ddsY <- dds[, yale]; colData(ddsY)$cond <- droplevels(colData(ddsY)$cond)
ddsY <- DESeq(ddsY, quiet=TRUE); saveRDS(ddsY, "v2/ddsY_v2.rds")

wp <- sub("_[0-9]+$","", sub(".*_cds_","", rownames(dds)))
ct <- function(o,a,b){ d <- as.data.frame(lfcShrink(o, contrast=c("cond",a,b), type="ashr", quiet=TRUE)); d$wp <- wp; d }

write.csv(ct(dds ,"Pl37","Pl25"), "v2/res_rG.csv")
write.csv(ct(ddsY,"Pl37","Pl25"), "v2/res_rY.csv")
write.csv(ct(ddsY,"Pl30","Pl25"), "v2/res_r30.csv")
write.csv(ct(ddsY,"Pl37","Pl30"), "v2/res_r3730.csv")
write.csv(ct(ddsY,"pH06","pH08"), "v2/res_rpH.csv")
cat("\nCinco contrastes escritos en v2/\n")
