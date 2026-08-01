# Enriquecimiento por dos metodos: hipergeometrico y sobre el ranking completo
suppressMessages({library(DESeq2); library(ashr); library(fgsea)})
setwd("~/bbrna")

dds <- readRDS("dds_global.rds")
map <- read.delim("lt2prot.tsv", header=FALSE, stringsAsFactors=FALSE)
names(map) <- c("lt","prot"); map$lt <- sub("BARBAKC583_","",map$lt)
wp <- sub("_[0-9]+$","", sub(".*_cds_","", rownames(dds)))
yale <- colData(dds)$cond %in% c("Pl25","Pl30","Pl37","pH06","pH07","pH08")
ddsY <- dds[, yale]; colData(ddsY)$cond <- droplevels(colData(ddsY)$cond)
ddsY <- DESeq(ddsY, quiet=TRUE)
b <- lfcShrink(ddsY, contrast=c("cond","Pl37","Pl25"), type="ashr")

cds <- read.delim("cds_todos.tsv", header=FALSE, stringsAsFactors=FALSE)
names(cds) <- c("lt","ini","fin","hebra","prod")
cds$lt <- sub("BARBAKC583_","",cds$lt)
cds$wp <- map$prot[match(cds$lt, map$lt)]

# Universo con identificadores de proteina unicos (evita el doble recuento
# de la region duplicada de ~28 kb presente en la referencia)
u <- cds[!duplicated(cds$wp) & !is.na(cds$wp) & cds$wp!="", ]
u$i <- match(u$wp, wp)
u <- u[!is.na(u$i) & !is.na(b$padj[u$i]), ]

fl <- sub("BARBAKC583_","", read.delim("flag_gff.tsv", header=FALSE)[,1])
rk <- b$log2FoldChange[u$i]; names(rk) <- u$lt; rk <- sort(rk, decreasing=TRUE)

vias <- list(
  Flagelar   = intersect(fl, names(rk)),
  Ribosoma   = u$lt[grepl("ribosomal protein", u$prod, ignore.case=TRUE)],
  Chaperona  = u$lt[grepl("chaperon|heat shock|GroE|DnaJ|DnaK|ClpB|HtpX", u$prod, ignore.case=TRUE)],
  TranspABC  = u$lt[grepl("ABC transporter", u$prod, ignore.case=TRUE)],
  HierroHemo = u$lt[grepl("hemin|heme|iron|TonB|ferr", u$prod, ignore.case=TRUE)],
  Traduccion = u$lt[grepl("tRNA|elongation factor|translation", u$prod, ignore.case=TRUE)])

set.seed(42)
res <- fgsea(pathways=vias, stats=rk, minSize=3, nPermSimple=100000)
print(res[order(res$pval), c("pathway","size","NES","pval","padj")])
