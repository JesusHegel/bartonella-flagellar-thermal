# 08. Figuras del proyecto (version de trabajo; se rehacen en la Fase 3).
suppressMessages({library(DESeq2); library(ggplot2)})
source("scripts/rutas.R")

ROJO <- "#B2182B"; AZUL <- "#2166AC"; GRIS <- "#BDBDBD"; TINTA <- "#1A1A1A"
w  <- function(x, n=74) paste(strwrap(x, width=n), collapse="\n")
co <- function(x) sub("\\.", ",", x)

tema <- theme_bw(base_size=13) + theme(
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(colour="grey93", linewidth=.35),
  panel.border     = element_rect(colour="grey30", linewidth=.6),
  axis.title       = element_text(face="bold", colour=TINTA),
  axis.text        = element_text(colour="grey20"),
  plot.title       = element_text(face="bold", size=14, colour=TINTA, lineheight=1.1),
  plot.subtitle    = element_text(size=10.5, colour="grey35", lineheight=1.18),
  plot.caption     = element_text(size=8.6, colour="grey45", hjust=0),
  plot.title.position = "plot", plot.caption.position = "plot",
  legend.position  = "top", legend.title = element_blank(),
  legend.key.height= grid::unit(4,"mm"),
  plot.margin      = margin(10,16,8,10))

dds <- readRDS(int("dds_23_muestras.rds"))
rY  <- read.csv(res("expresion_diferencial", "37C_vs_25C.csv"), row.names=1)
r30 <- read.csv(res("expresion_diferencial", "30C_vs_25C.csv"), row.names=1)
r33 <- read.csv(res("expresion_diferencial", "37C_vs_30C.csv"), row.names=1)
wp  <- rY$wp
map <- read.delim(res("anotacion", "locus_a_proteina.tsv"), header=FALSE); names(map) <- c("lt","prot")
cds <- read.delim(res("anotacion", "anotacion_cds.tsv"), header=FALSE)
names(cds) <- c("lt","ini","fin","hebra","prod")
fl  <- read.delim(res("anotacion", "conjunto_flagelar.tsv"), header=FALSE)
names(fl) <- c("lt","ini","fin","hebra","prod")
g   <- match(map$prot[match(fl$lt, map$lt)], wp)

nc <- counts(dds, normalized=TRUE); cn <- as.character(colData(dds)$cond)
mn <- sapply(c("Pl25","Pl30","Pl37"), function(k) rowMeans(nc[, cn==k, drop=FALSE]))

# ============== FIGURA 5: DISOCIACION TERMICA ==============
ch_lt <- cds$lt[grepl("chaperon|heat shock|GroE|DnaJ|DnaK|ClpB|HtpX", cds$prod, ignore.case=TRUE)]
gch   <- match(map$prot[match(ch_lt, map$lt)], wp); gch <- gch[!is.na(gch)]
idx <- function(i, grupo){
  m <- mn[i, , drop=FALSE]; m <- m[m[,"Pl25"] > 10 | m[,"Pl37"] > 10, , drop=FALSE]
  do.call(rbind, lapply(c(25,30,37), function(tt){
    col <- c("25"="Pl25","30"="Pl30","37"="Pl37")[as.character(tt)]
    data.frame(gen=rownames(m), temp=tt,
               lfc=log2((m[,col]+1)/(m[,"Pl25"]+1)), grupo=grupo, row.names=NULL)}))}
nf <- sum(!is.na(g)); nch <- length(gch)
GF <- sprintf("Aparato flagelar (n=%d)", nf); GC <- sprintf("Chaperonas (n=%d)", nch)
d5  <- rbind(idx(g[!is.na(g)], GF), idx(gch, GC))
med <- aggregate(lfc ~ grupo + temp, d5, median)
lab <- med[med$temp==37,]; lab$corto <- ifelse(lab$grupo==GF, "Flagelar", "Chaperonas")
yr  <- range(d5$lfc)

p5 <- ggplot(d5, aes(temp, lfc, group=gen, colour=grupo)) +
  geom_hline(yintercept=0, colour="grey55", linewidth=.5) +
  geom_line(linewidth=.32, alpha=.28) +
  geom_line(data=med, aes(temp, lfc, group=grupo, colour=grupo), linewidth=2.2, inherit.aes=FALSE) +
  geom_point(data=med, aes(temp, lfc, colour=grupo), size=3.8, inherit.aes=FALSE) +
  geom_text(data=lab, aes(x=37.22, y=lfc, label=corto, colour=grupo),
            hjust=0, fontface="bold", size=4.3, inherit.aes=FALSE, show.legend=FALSE) +
  annotate("text", x=27.5, y=yr[1]*0.96, hjust=.5, size=3.9, fontface="italic",
           colour=ROJO, label="El flagelo se apaga aquí") +
  annotate("text", x=33.5, y=yr[2]*0.90, hjust=.5, size=3.9, fontface="italic",
           colour=AZUL, label="Las chaperonas siguen subiendo") +
  scale_colour_manual(values=setNames(c(ROJO,AZUL), c(GF,GC))) +
  scale_x_continuous(breaks=c(25,30,37), labels=c("25 °C","30 °C","37 °C"),
                     expand=expansion(mult=c(.05,.20))) +
  coord_cartesian(clip="off") +
  labs(title="Dos programas térmicos con perfiles distintos",
       subtitle=w("Cada línea fina es un gen; la línea gruesa es la mediana del grupo. La represión flagelar se completa a 30 °C y no avanza más allá; la respuesta de chaperonas sigue aumentando hasta 37 °C."),
       x="Temperatura de cultivo",
       y=expression(bold(log[2]~"(conteo / conteo a 25 °C)")),
       caption="Conteos normalizados por DESeq2, indexados a 25 °C. PRJNA647605, modelo restringido (12 muestras).") +
  tema + theme(plot.margin=margin(10,34,8,10))
ggsave(fig("curva_termica.png"), p5, width=235, height=145, units="mm", dpi=300, bg="white")
ggsave(fig("curva_termica.pdf"), p5, width=235, height=145, units="mm", bg="white")
cat("\n--- Fig5: medianas por grupo ---\n"); print(med)

# ============== FIGURA 1: PCA ==============
vsd <- vst(dds, blind=TRUE)
pv  <- plotPCA(vsd, intgroup="cond", ntop=500, returnData=TRUE)
pct <- round(100*attr(pv,"percentVar"))
pm  <- read.csv(res("muestras", "pca_y_tasa_mapeo.csv"))
pv$run <- rownames(pv); pv <- merge(pv, pm[,c("run","mapeo","matriz")], by="run")
rho <- co(sprintf("%.2f", cor(pv$PC1, pv$mapeo)))
niv <- c("placa / caldo","placa + atmosfera sanguinea","sangre humana","celulas endoteliales")
eti <- c("Placa / caldo  (n=12)","Placa + atmósfera sanguínea  (n=3)",
         "Sangre humana  (n=6)","Células endoteliales  (n=2)")
pv$mat <- factor(eti[match(pv$matriz, niv)], levels=eti)
p1 <- ggplot(pv, aes(PC1, PC2, colour=mat, shape=mat)) +
  geom_point(size=4.3, alpha=.92, stroke=1.1) +
  scale_colour_manual(values=setNames(c(AZUL,"#6BAED6",ROJO,"#7B3294"), eti)) +
  scale_shape_manual(values=setNames(c(16,17,15,18), eti)) +
  labs(title="La mayor fuente de variación es la matriz biológica",
       subtitle=w(sprintf("PC1 correlaciona %s con la fracción de lecturas asignadas a genes bacterianos. Por eso los contrastes se restringen a las 12 muestras de matriz homogénea.", rho)),
       x=sprintf("PC1 (%d %% de la varianza)", pct[1]),
       y=sprintf("PC2 (%d %% de la varianza)", pct[2]),
       caption="23 corridas, un solo centro de secuenciación (University of Montana, HiSeq 2500).") +
  tema + theme(legend.position="right", plot.margin=margin(10,10,8,10))
ggsave(fig("pca_y_tasa_mapeo.png"), p1, width=245, height=135, units="mm", dpi=300, bg="white")
ggsave(fig("pca_y_tasa_mapeo.pdf"), p1, width=245, height=135, units="mm", bg="white")

# ============== FIGURA 2: VOLCANO ==============
vv <- data.frame(lfc=rY$log2FoldChange, padj=rY$padj, flag=FALSE)
vv$flag[g[!is.na(g)]] <- TRUE
vv <- vv[!is.na(vv$padj) & vv$padj > 0, ]
vv$y <- -log10(vv$padj); ymax <- max(vv$y)
p2 <- ggplot(vv[!vv$flag,], aes(lfc, y)) +
  geom_point(colour=GRIS, size=1.3, alpha=.55) +
  geom_point(data=vv[vv$flag,], aes(lfc, y), colour=ROJO, size=2.7, alpha=.95) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", colour="grey60", linewidth=.45) +
  geom_hline(yintercept=-log10(.05), linetype="dashed", colour="grey60", linewidth=.45) +
  annotate("text", x=min(vv$lfc)*.95, y=ymax*.97, hjust=0, size=4.2,
           fontface="bold", colour=ROJO, label=sprintf("Genes flagelares (n=%d)", nf)) +
  annotate("text", x=min(vv$lfc)*.95, y=ymax*.88, hjust=0, size=3.6,
           fontface="italic", colour="grey40", label="Líneas discontinuas: q = 0,05 y |cambio| = 1") +
  labs(title="37 °C frente a 25 °C: el aparato flagelar ocupa el extremo reprimido",
       subtitle=w("Modelo restringido a las 12 muestras de matriz homogénea, con encogimiento ashr."),
       x=expression(bold(log[2]~"del cambio  (37 °C / 25 °C)")),
       y=expression(bold(-log[10]*"(q)"))) + tema
ggsave(fig("volcan_37C_vs_25C.png"), p2, width=225, height=145, units="mm", dpi=300, bg="white")
ggsave(fig("volcan_37C_vs_25C.pdf"), p2, width=225, height=145, units="mm", bg="white")

# ============== FIGURA 3: MAPA DE CALOR ==============
r  <- regexpr("(Fli[A-Z]|Flg[A-Z]|Flh[A-Z]|Flb[A-Z]|Fla[A-Z]|Mot[A-Z]|Che[A-Z])", fl$prod)
sy <- rep(NA_character_, nrow(fl)); sy[r>0] <- regmatches(fl$prod, r)
sy[is.na(sy) & grepl("flagellin", fl$prod, ignore.case=TRUE)] <- "Flagelina"
sy[is.na(sy)] <- "Gancho-asoc."
nom <- paste0(sy, "  ·  ", sub("BARBAKC583_","", fl$lt))
h <- rbind(data.frame(gen=nom, ctr="30 vs 25", lfc=r30$log2FoldChange[g]),
           data.frame(gen=nom, ctr="37 vs 25", lfc=rY$log2FoldChange[g]),
           data.frame(gen=nom, ctr="37 vs 30", lfc=r33$log2FoldChange[g]))
h$gen <- factor(h$gen, levels=nom[order(rY$log2FoldChange[g])])
h$ctr <- factor(h$ctr, levels=c("30 vs 25","37 vs 25","37 vs 30"))
p3 <- ggplot(h, aes(ctr, gen, fill=lfc)) +
  geom_tile(colour="white", linewidth=1) +
  geom_text(aes(label=sprintf("%.1f", lfc)), size=2.9,
            colour=ifelse(abs(h$lfc) > 1.5, "white", "grey25")) +
  scale_fill_gradient2(low=ROJO, mid="grey96", high=AZUL, midpoint=0,
                       name=expression(log[2]~"cambio")) +
  labs(title="Los 31 genes flagelares, gen por gen",
       subtitle=w("La tercera columna está casi vacía: entre 30 y 37 °C ya no ocurre nada.", 60),
       x=NULL, y=NULL) +
  tema + theme(legend.position="right",
               legend.title=element_text(size=9.5, face="bold"),
               axis.text.y=element_text(size=8.6, family="mono"),
               axis.text.x=element_text(size=11, face="bold"),
               panel.grid.major=element_blank(),
               plot.margin=margin(10,12,8,10))
ggsave(fig("mapa_calor_flagelares.png"), p3, width=185, height=225, units="mm", dpi=300, bg="white")
ggsave(fig("mapa_calor_flagelares.pdf"), p3, width=185, height=225, units="mm", bg="white")

# ============== FIGURA 4: PANEL DE GENERO ==============
pn <- read.delim(res("genomica_comparada", "panel_genero_blastp.tsv"), header=FALSE); names(pn) <- c("gen","especie","ident")
cnt <- as.data.frame(table(unique(pn[,c("gen","especie")])$especie)); names(cnt) <- c("especie","n")
cnt$especie <- factor(cnt$especie, levels=cnt$especie[order(cnt$n)])
p4 <- ggplot(cnt, aes(n, especie)) +
  geom_col(fill=AZUL, width=.66) +
  geom_text(aes(label=paste0(n, " / 31")), hjust=-0.2, size=4.1, fontface="bold", colour=TINTA) +
  scale_x_continuous(limits=c(0,37), expand=expansion(mult=c(0,.02))) +
  labs(title=w("El regulón flagelar está presente o ausente, sin estados intermedios", 52),
       subtitle=w("Homólogos de los 31 genes por BLASTp (e < 1e-5) en siete especies del género."),
       x="Genes flagelares detectados", y=NULL) +
  tema + theme(axis.text.y=element_text(face="italic", size=11.5),
               panel.grid.major.y=element_blank())
ggsave(fig("panel_genero.png"), p4, width=215, height=125, units="mm", dpi=300, bg="white")
ggsave(fig("panel_genero.pdf"), p4, width=215, height=125, units="mm", bg="white")

cat("\nFiguras escritas:\n"); print(list.files(DIR_FIG))
