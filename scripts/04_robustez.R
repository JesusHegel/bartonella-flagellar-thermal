suppressMessages({library(DESeq2)})
setwd("~/bbrna")
dds <- readRDS("salida/dds.rds"); ddsY <- readRDS("salida/dds_lote.rds")
rG<-read.csv("salida/res_rG.csv",row.names=1); rY<-read.csv("salida/res_rY.csv",row.names=1)
r30<-read.csv("salida/res_r30.csv",row.names=1); r33<-read.csv("salida/res_r3730.csv",row.names=1)
rp<-read.csv("salida/res_rpH.csv",row.names=1); wp<-rG$wp

map<-read.delim("salida/lt2prot.tsv",header=FALSE); names(map)<-c("lt","prot")
fl <-read.delim("salida/flag_gff.tsv",header=FALSE); names(fl)<-c("lt","ini","fin","hebra","prod")
g  <- match(map$prot[match(fl$lt,map$lt)], wp)
cat("Flagelares localizados en la matriz:", sum(!is.na(g)), "de", nrow(fl), "\n\n")

cnt <- function(r,lab) data.frame(definicion=lab,
  genoma_padj05      = sum(r$padj<0.05, na.rm=TRUE),
  genoma_padj05_lfc1 = sum(r$padj<0.05 & abs(r$log2FoldChange)>1, na.rm=TRUE),
  flag_padj05        = sum(r$padj[g]<0.05, na.rm=TRUE),
  flag_padj05_lfc1   = sum(r$padj[g]<0.05 & abs(r$log2FoldChange[g])>1, na.rm=TRUE))
tb <- rbind(cnt(rG,"global 37vs25"), cnt(rY,"restringido 37vs25"),
            cnt(r30,"restringido 30vs25"), cnt(rp,"restringido pH06vspH08"))
cat("--- CONTEOS BAJO CADA DEFINICION ---\n"); print(tb)
write.csv(tb,"salida/chk_conteos.csv",row.names=FALSE)

fr <- as.data.frame(results(ddsY, contrast=c("cond","Pl37","Pl25"), lfcThreshold=1))
est <- data.frame(criterio=c("post hoc padj<0.05 & |LFC|>1","lfcThreshold=1 formal"),
  n_significativos=c(sum(rY$padj<0.05 & abs(rY$log2FoldChange)>1,na.rm=TRUE), sum(fr$padj<0.05,na.rm=TRUE)),
  n_flagelares=c(sum(rY$padj[g]<0.05 & abs(rY$log2FoldChange[g])>1,na.rm=TRUE), sum(fr$padj[g]<0.05,na.rm=TRUE)))
cat("\n--- SENSIBILIDAD AL UMBRAL ---\n"); print(est)
write.csv(est,"salida/chk_estricto.csv",row.names=FALSE)

yl <- data.frame(locus=fl$lt, producto=fl$prod,
  lfc_global=round(rG$log2FoldChange[g],2), padj_global=signif(rG$padj[g],2),
  lfc_restr =round(rY$log2FoldChange[g],2), padj_restr =signif(rY$padj[g],2))
yl$dif <- round(yl$lfc_global-yl$lfc_restr,2); yl <- yl[order(yl$lfc_restr),]
cat("\n--- GLOBAL vs RESTRINGIDO ---\n"); print(yl[,c("locus","lfc_global","lfc_restr","dif")])
cat("\n*** Diferencia maxima entre modelos:", max(abs(yl$dif),na.rm=TRUE), "***\n")
write.csv(yl,"salida/chk_yale.csv",row.names=FALSE)

c31 <- data.frame(locus=fl$lt, producto=fl$prod,
  L30v25=round(r30$log2FoldChange[g],2), p30v25=signif(r30$padj[g],2),
  L37v25=round(rY$log2FoldChange[g],2),  p37v25=signif(rY$padj[g],2),
  L37v30=round(r33$log2FoldChange[g],2), p37v30=signif(r33$padj[g],2))
c31 <- c31[order(c31$L37v25),]
cat("\n--- CURVA 25 / 30 / 37 ---\n"); print(c31[,c("locus","L30v25","L37v25","L37v30","p37v30")])
write.csv(c31,"salida/chk_curva31.csv",row.names=FALSE)

cat("\n--- ESPECIFICIDAD: temperatura vs pH ---\n")
ef <- data.frame(contraste=c("temperatura 37vs25","pH 06vs08"),
  flagelares_sig=c(sum(rY$padj[g]<0.05,na.rm=TRUE), sum(rp$padj[g]<0.05,na.rm=TRUE)),
  genoma_sig    =c(sum(rY$padj<0.05,na.rm=TRUE),    sum(rp$padj<0.05,na.rm=TRUE)))
print(ef); write.csv(ef,"salida/chk_especificidad.csv",row.names=FALSE)

nc <- counts(dds, normalized=TRUE); cn <- as.character(colData(dds)$cond)
i25<-which(cn=="Pl25"); i37<-which(cn=="Pl37")
cat("\nReplicas Pl25:",length(i25)," Pl30:",sum(cn=="Pl30")," Pl37:",length(i37),"\n")
if (length(i25)>=2 && length(i37)>=2){
  lr<-function(a,b) log2((nc[,b]+1)/(nc[,a]+1))
  A<-lr(i25[1],i37[1]); B<-lr(i25[2],i37[2])
  pt<-data.frame(locus=fl$lt, particionA=round(A[g],2), particionB=round(B[g],2))
  pt$misma_direccion <- sign(pt$particionA)==sign(pt$particionB)
  cat("\n--- PARTICIONES INDEPENDIENTES ---\n")
  cat("Misma direccion:", sum(pt$misma_direccion,na.rm=TRUE), "de", nrow(pt), "\n")
  cat("Correlacion genoma completo:", round(cor(A,B,use="complete.obs"),3), "\n")
  write.csv(pt,"salida/chk_particiones.csv",row.names=FALSE)
}

mn <- sapply(c("Pl25","Pl30","Pl37"), function(k) round(rowMeans(nc[,cn==k,drop=FALSE])))
lk <- data.frame(wp=wp, lt=map$lt[match(wp,map$prot)])
lk$prod <- read.delim("salida/cds_todos.tsv",header=FALSE)[,5][match(lk$lt, read.delim("salida/cds_todos.tsv",header=FALSE)[,1])]
ab <- cbind(lk, mn)
cat("\n--- 15 MAS ABUNDANTES A 37 C ---\n")
print(head(ab[order(-ab$Pl37), c("lt","prod","Pl25","Pl30","Pl37")], 15))
write.csv(head(ab[order(-ab$Pl37),],15), "salida/chk_ranking.csv", row.names=FALSE)
write.csv(head(ab[order(-rowMeans(mn)),],15), "salida/chk_top_global.csv", row.names=FALSE)
cat("\n--- FLAGELINA RS05045 a 25/30/37 ---\n"); print(ab[ab$lt=="BARBAKC583_RS05045", c("Pl25","Pl30","Pl37")])
cat("--- GroEL / chaperonas ---\n"); print(head(ab[grepl("GroEL|chaperonin", ab$prod, ignore.case=TRUE), c("lt","prod","Pl25","Pl30","Pl37")],5))
