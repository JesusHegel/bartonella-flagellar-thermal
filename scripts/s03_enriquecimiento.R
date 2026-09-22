suppressMessages({library(fgsea)})
setwd("~/bbrna")

b   <- read.csv("v2/res_rY.csv", row.names=1); wp <- b$wp
cds <- read.delim("v2/cds_todos.tsv", header=FALSE, stringsAsFactors=FALSE)
names(cds) <- c("lt","ini","fin","hebra","prod")
map <- read.delim("v2/lt2prot.tsv", header=FALSE, stringsAsFactors=FALSE)
names(map) <- c("lt","prot")
cds$wp <- map$prot[match(cds$lt, map$lt)]

u <- cds[!duplicated(cds$wp) & !is.na(cds$wp) & cds$wp!="", ]
u$i <- match(u$wp, wp)
u <- u[!is.na(u$i) & !is.na(b$padj[u$i]), ]
N <- nrow(u)

fl  <- read.delim("v2/flag_gff.tsv", header=FALSE)[,1]
lfc <- b$log2FoldChange[u$i]; pad <- b$padj[u$i]
sig <- !is.na(pad) & pad < 0.05
ind <- sig & lfc > 0; rep <- sig & lfc < 0

cat("*** Universo deduplicado:", N, " inducidos:", sum(ind), " reprimidos:", sum(rep), "***\n\n")

vias <- list(
  Flagelar   = u$lt[u$lt %in% fl],
  Ribosoma   = u$lt[grepl("ribosomal protein", u$prod, ignore.case=TRUE)],
  Chaperona  = u$lt[grepl("chaperon|heat shock|GroE|DnaJ|DnaK|ClpB|HtpX", u$prod, ignore.case=TRUE)],
  TranspABC  = u$lt[grepl("ABC transporter", u$prod, ignore.case=TRUE)],
  HierroHemo = u$lt[grepl("hemin|heme|iron|TonB|ferr", u$prod, ignore.case=TRUE)],
  Traduccion = u$lt[grepl("tRNA|elongation factor|translation", u$prod, ignore.case=TRUE)],
  Hipotetica = u$lt[grepl("hypothetical", u$prod, ignore.case=TRUE)])

hy <- do.call(rbind, lapply(names(vias), function(v){
  g <- u$lt %in% vias[[v]]; K <- sum(g); ki <- sum(g&ind); kr <- sum(g&rep)
  data.frame(categoria=v, genes=K, inducidos=ki, reprimidos=kr,
             esp_red=round(K*sum(rep)/N,1),
             p_ind=phyper(ki-1,K,N-K,sum(ind),lower.tail=FALSE),
             p_rep=phyper(kr-1,K,N-K,sum(rep),lower.tail=FALSE))}))
hy$q_ind <- p.adjust(hy$p_ind,"BH"); hy$q_rep <- p.adjust(hy$p_rep,"BH")
cat("--- HIPERGEOMETRICO ---\n"); print(hy)
write.csv(hy, "v2/chk_categorias.csv", row.names=FALSE)

rk <- lfc; names(rk) <- u$lt; rk <- sort(rk[!is.na(rk)], decreasing=TRUE)
set.seed(42)
fg <- fgsea(pathways=vias, stats=rk, minSize=3, nPermSimple=100000)
fg <- as.data.frame(fg[order(fg$pval), c("pathway","size","NES","pval","padj")])
cat("\n--- FGSEA (semilla 42, 100000 permutaciones) ---\n"); print(fg)
write.csv(fg, "v2/chk_fgsea.csv", row.names=FALSE)
write.csv(data.frame(universo=N, inducidos=sum(ind), reprimidos=sum(rep)),
          "v2/chk_universo.csv", row.names=FALSE)
