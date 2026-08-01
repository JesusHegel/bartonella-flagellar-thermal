# Dos verificaciones de robustez del resultado principal:
#  (a) concordancia entre particiones independientes por replica
#  (b) especificidad frente al contraste de pH
suppressMessages({library(DESeq2); library(ashr)})
setwd("~/bbrna")

dds <- readRDS("dds_global.rds")
map <- read.delim("lt2prot.tsv", header=FALSE, stringsAsFactors=FALSE)
names(map) <- c("lt","prot"); map$lt <- sub("BARBAKC583_","",map$lt)
wp <- sub("_[0-9]+$","", sub(".*_cds_","", rownames(dds)))
fl <- sub("BARBAKC583_","", read.delim("flag_gff.tsv", header=FALSE)[,1])
g  <- match(map$prot[match(fl, map$lt)], wp)

# (a) Particiones independientes: una sola replica por condicion en cada una
nc <- counts(dds, normalized=TRUE)
cn <- as.character(colData(dds)$cond)
i25 <- which(cn=="Pl25"); i37 <- which(cn=="Pl37")
lr <- function(a,b) log2((nc[,b]+1)/(nc[,a]+1))
A <- lr(i25[1], i37[1])   # particion A
B <- lr(i25[2], i37[2])   # particion B
cat("Genes flagelares con misma direccion en ambas particiones:",
    sum(sign(A[g])==sign(B[g])), "de", length(g), "\n")
cat("Correlacion entre particiones (genoma completo):",
    round(cor(A, B, use="complete.obs"), 3), "\n")

# (b) Especificidad: respuesta a temperatura frente a respuesta a pH
te <- lfcShrink(dds, contrast=c("cond","Pl37","Pl25"), type="ashr")
ph <- lfcShrink(dds, contrast=c("cond","pH06","pH08"), type="ashr")
cat("\nFlagelares significativos a temperatura:", sum(te$padj[g]<0.05, na.rm=TRUE), "de 31\n")
cat("Flagelares significativos a pH         :", sum(ph$padj[g]<0.05, na.rm=TRUE), "de 31\n")
cat("Fondo del genoma  temperatura:", sum(te$padj<0.05, na.rm=TRUE),
    " pH:", sum(ph$padj<0.05, na.rm=TRUE), "\n")
