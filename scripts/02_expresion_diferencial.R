# Expresion diferencial: modelo global y modelo restringido a un lote tecnico
suppressMessages({library(DESeq2); library(ashr)})
setwd("~/bbrna")

dds <- readRDS("dds_global.rds")
map <- read.delim("lt2prot.tsv", header=FALSE, stringsAsFactors=FALSE)
names(map) <- c("lt","prot"); map$lt <- sub("BARBAKC583_","",map$lt)
wp <- sub("_[0-9]+$","", sub(".*_cds_","", rownames(dds)))
fl <- sub("BARBAKC583_","", read.delim("flag_gff.tsv", header=FALSE)[,1])
g  <- match(map$prot[match(fl, map$lt)], wp)

# Modelo global (23 muestras)
rG <- lfcShrink(dds, contrast=c("cond","Pl37","Pl25"), type="ashr")

# Modelo restringido al lote tecnico homogeneo (12 muestras)
yale <- colData(dds)$cond %in% c("Pl25","Pl30","Pl37","pH06","pH07","pH08")
ddsY <- dds[, yale]; colData(ddsY)$cond <- droplevels(colData(ddsY)$cond)
ddsY <- DESeq(ddsY, quiet=TRUE)
rY <- lfcShrink(ddsY, contrast=c("cond","Pl37","Pl25"), type="ashr")

write.csv(data.frame(locus=fl,
  lfc_global=round(rG$log2FoldChange[g],2), padj_global=signif(rG$padj[g],2),
  lfc_yale  =round(rY$log2FoldChange[g],2), padj_yale  =signif(rY$padj[g],2)),
  "chk_yale.csv", row.names=FALSE)
