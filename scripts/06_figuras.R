# ---------------------------------------------------------------
# Figuras del manuscrito
#
# Cuatro figuras principales (Fig1-Fig4) y tres versiones
# simplificadas para presentacion oral. Las versiones simplificadas
# representan los mismos datos con menor densidad de informacion:
# no aplican ningun filtro ni transformacion adicional.
# ---------------------------------------------------------------
suppressMessages({library(DESeq2); library(ashr); library(ggplot2)
                  library(ggrepel); library(pheatmap); library(RColorBrewer)})
setwd("~/bbrna"); dir.create("figuras", showWarnings=FALSE)

dds <- readRDS("dds_global.rds")
map <- read.delim("lt2prot.tsv", header=FALSE, stringsAsFactors=FALSE)
names(map) <- c("lt","prot"); map$lt <- sub("BARBAKC583_","",map$lt)
wp  <- sub("_[0-9]+$","", sub(".*_cds_","", rownames(dds)))
cds <- read.delim("cds_todos.tsv", header=FALSE, stringsAsFactors=FALSE)
names(cds) <- c("lt","ini","fin","hebra","prod")
cds$lt <- sub("BARBAKC583_","",cds$lt)
fl <- sub("BARBAKC583_","", read.delim("flag_gff.tsv", header=FALSE)[,1])

# ============ FIGURA 1: PCA ============
vsd <- vst(dds, blind=TRUE); saveRDS(vsd, "vsd.rds")
pv  <- plotPCA(vsd, intgroup="cond", returnData=TRUE)
pct <- round(100*attr(pv,"percentVar"))
yale <- c("Pl25","Pl30","Pl37","pH06","pH07","pH08")
pv$lote <- ifelse(pv$cond %in% yale, "Group A", "Group B")
pv$cond <- factor(as.character(pv$cond),
  levels=c("Pl25","Pl30","Pl37","pH06","pH07","pH08","PlBG","HB37","HBBG","HUVE"))
brecha <- (min(pv$PC1[pv$lote=="Group A"]) + max(pv$PC1[pv$lote=="Group B"]))/2

p1 <- ggplot(pv, aes(PC1, PC2, colour=cond, shape=lote)) +
  geom_vline(xintercept=brecha, linetype="dashed", colour="grey45", linewidth=.7) +
  geom_point(size=4.4, alpha=.95, stroke=1.1) +
  scale_shape_manual(values=c("Group A"=16,"Group B"=15),
    labels=c("Group A"="Group A (single centre)",
             "Group B"="Group B (blood / cell samples)"), name=NULL) +
  scale_colour_manual(values=c(
    "Pl25"="#08306b","Pl30"="#4292c6","Pl37"="#9ecae1",
    "pH06"="#00441b","pH07"="#41ab5d","pH08"="#a1d99b",
    "PlBG"="#a63603","HB37"="#e6550d","HBBG"="#fdae6b","HUVE"="#54278f"),
    name="Condition") +
  annotate("text", x=brecha, y=max(pv$PC2)*1.13, label="No overlap between groups",
           size=4.3, fontface="bold", colour="grey25", hjust=0.5) +
  annotate("text", x=min(pv$PC1)*0.97, y=min(pv$PC2)*1.15,
           label="PC2 orders plate cultures by temperature",
           size=4, fontface="italic", colour="grey35", hjust=0) +
  labs(x=paste0("PC1 (", pct[1], "% of variance)  \u2014  technical axis"),
       y=paste0("PC2 (", pct[2], "% of variance)")) +
  theme_bw(base_size=14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour="grey93"),
        panel.border=element_rect(colour="grey25", linewidth=.7),
        axis.title=element_text(size=13, face="bold", colour="grey10"),
        axis.text=element_text(size=11.5, colour="grey15"),
        legend.text=element_text(size=10.5),
        legend.title=element_text(size=11, face="bold"),
        legend.key.height=grid::unit(4.6,"mm"),
        plot.margin=margin(10,6,8,8))
ggsave("figuras/Fig1_PCA.png", p1, width=235, height=125, units="mm", dpi=300, bg="white")

# Version simplificada para presentacion: mismos 23 puntos, dos categorias
pv$grupo <- ifelse(pv$lote=="Group A","Group A (single centre)",
                                      "Group B (blood / cell samples)")
p1b <- ggplot(pv, aes(PC1, PC2, colour=grupo, shape=grupo)) +
  geom_vline(xintercept=brecha, linetype="dashed", colour="grey40", linewidth=.9) +
  geom_point(size=6, alpha=.9, stroke=1.3) +
  scale_colour_manual(values=c("Group A (single centre)"="#08519c",
                               "Group B (blood / cell samples)"="#d94801"), name=NULL) +
  scale_shape_manual(values=c("Group A (single centre)"=16,
                              "Group B (blood / cell samples)"=15), name=NULL) +
  annotate("text", x=brecha, y=max(pv$PC2)*1.15, label="No overlap",
           size=6, fontface="bold", colour="grey20", hjust=0.5) +
  labs(x=paste0("PC1 (", pct[1], "% of variance)  \u2014  technical axis"),
       y=paste0("PC2 (", pct[2], "%)")) +
  theme_bw(base_size=18) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour="grey93"),
        panel.border=element_rect(colour="grey25", linewidth=.8),
        axis.title=element_text(size=16, face="bold", colour="grey10"),
        axis.text=element_text(size=14, colour="grey15"),
        legend.position="top", legend.text=element_text(size=15),
        plot.margin=margin(12,14,10,10))
ggsave("figuras/Fig1_PCA_simple.png", p1b, width=215, height=130, units="mm", dpi=300, bg="white")

# ============ FIGURA 2: volcano ============
r <- lfcShrink(dds, contrast=c("cond","Pl37","Pl25"), type="ashr", quiet=TRUE)
d <- data.frame(lt=map$lt[match(wp, map$prot)], lfc=r$log2FoldChange, padj=r$padj)
d <- d[!is.na(d$padj) & !is.na(d$lt), ]
d$prod <- cds$prod[match(d$lt, cds$lt)]
d$grp <- "Other"
d$grp[grepl("chaperon|heat shock|GroE|DnaJ|DnaK|ClpB|HtpX", d$prod, ignore.case=TRUE)] <- "Chaperone"
d$grp[d$lt %in% fl] <- "Flagellar"
# Se etiquetan los transcritos mencionados en el texto
lab <- c("RS05045","RS05430","RS04185","RS05735","RS06465","RS05920")
nom <- c("flagellin","FlbT","FliO","GroEL","DnaK","Pap31")
d$eti <- nom[match(d$lt, lab)]
pap <- d[d$lt == "RS05920", ]   # Pap31 no es flagelar: se marca su borde

p2 <- ggplot(d[order(match(d$grp, c("Other","Chaperone","Flagellar"))), ],
             aes(lfc, -log10(padj))) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey50", linewidth=.55) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", colour="grey50", linewidth=.55) +
  geom_point(aes(colour=grp, size=grp, alpha=grp)) +
  geom_point(data=pap, colour="grey30", fill="grey85", shape=21,
             size=3.4, stroke=1.0) +
  geom_text_repel(aes(label=ifelse(is.na(eti),"",eti), hjust=ifelse(lfc<0,1,0)),
    size=4.3, colour="grey10", fontface="bold", box.padding=1.2, point.padding=.5,
    nudge_x=ifelse(is.na(d$eti),0, ifelse(d$eti=="flagellin",-0.55,0)),
    force=3.5, force_pull=.25, min.segment.length=0, segment.colour="grey30",
    max.overlaps=Inf, max.iter=20000, seed=11, na.rm=TRUE) +
  annotate("text", x=-Inf, y=-Inf, label="\u2190 Repressed at 37 \u00b0C",
           hjust=-0.04, vjust=-0.75, size=3.9, fontface="bold", colour="grey20") +
  annotate("text", x=Inf, y=-Inf, label="Induced at 37 \u00b0C \u2192",
           hjust=1.04, vjust=-0.75, size=3.9, fontface="bold", colour="grey20") +
  annotate("text", x=Inf, y=-log10(0.05), label="adjusted p = 0.05",
           hjust=1.05, vjust=-0.65, size=3.3, colour="grey40") +
  scale_colour_manual(name="Transcript class", breaks=c("Flagellar","Chaperone","Other"),
    values=c("Flagellar"="#b2182b","Chaperone"="#2166ac","Other"="grey78"),
    labels=c("Flagellar"=sprintf("Flagellar genes (n = %d)", sum(d$grp=="Flagellar")),
             "Chaperone"=sprintf("Chaperones (n = %d)", sum(d$grp=="Chaperone")),
             "Other"=sprintf("Other transcripts (n = %s)",
                             format(sum(d$grp=="Other"), big.mark=",")))) +
  scale_size_manual(values=c("Other"=1.35,"Chaperone"=3.15,"Flagellar"=3.15), guide="none") +
  scale_alpha_manual(values=c("Other"=.55,"Chaperone"=.95,"Flagellar"=.95), guide="none") +
  guides(colour=guide_legend(override.aes=list(size=3.4, alpha=1))) +
  labs(x=expression(log[2]~"fold change (37 \u00b0C vs 25 \u00b0C)"),
       y=expression(-log[10]~"(adjusted "*italic(p)*")")) +
  theme_bw(base_size=13) +
  theme(panel.grid.major=element_line(colour="grey92", linewidth=.35),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(colour="grey25", linewidth=.7),
        axis.title=element_text(size=13, face="bold", colour="grey10"),
        axis.text=element_text(size=11.5, colour="grey15"),
        legend.position="inside", legend.position.inside=c(.03,.97),
        legend.justification=c(0,1),
        legend.title=element_text(size=10.5, face="bold"),
        legend.text=element_text(size=9.8),
        legend.background=element_rect(fill="white", colour="grey70", linewidth=.45),
        plot.margin=margin(8,12,8,8))
ggsave("figuras/Fig2_volcano.png", p2, width=200, height=130, units="mm", dpi=300, bg="white")

# ============ FIGURA 3: heatmap del regulon flagelar ============
nc <- counts(dds, normalized=TRUE); rownames(nc) <- wp
cn <- as.character(colData(dds)$cond)
keep <- which(cn %in% c("Pl25","Pl30","Pl37")); keep <- keep[order(cn[keep])]
m <- log2(nc[, keep] + 1)
colnames(m) <- paste0(sub("Pl","",cn[keep]), " \u00b0C\nrep ", rep(1:2,3))
mm <- m[match(map$prot[match(fl, map$lt)], rownames(m)), ]
pr <- cds$prod[match(fl, cds$lt)]
nm <- sub("^flagellar ","",pr); nm <- sub(" domain-containing protein","",nm)
nm <- sub("N-terminal helical ","",nm)
rownames(mm) <- paste0(fl,"  ",substr(nm,1,30))

# El orden de asignacion importa: motor se aplica despues de gancho
mod <- rep("Other", length(pr))
mod[grepl("FlgD|FlgE|FlgK|FlgL|hook", pr, ignore.case=TRUE)] <- "Hook"
mod[grepl("FlgB|FlgC|FlgF|FlgG|FliE|basal.?body|L-ring|P-ring|MS.?ring|FliF", pr, ignore.case=TRUE)] <- "Basal body"
mod[grepl("FliI|FliP|FliQ|FliR|FlhA|FlhB|export|secretion", pr, ignore.case=TRUE)] <- "T3SS export"
mod[grepl("MotA|MotB|MotC|FliG|FliM|FliN|motor|switch|stator", pr, ignore.case=TRUE)] <- "Motor"
mod[grepl("flagellin", pr, ignore.case=TRUE)] <- "Filament"
mod[grepl("FlbT|FlaF|regulator|repressor", pr, ignore.case=TRUE)] <- "Regulator"
ann <- data.frame(Module=mod); rownames(ann) <- rownames(mm)
o <- order(factor(mod, levels=c("Regulator","T3SS export","Basal body","Hook","Motor","Filament","Other")))
mm2 <- mm[o,]; ann2 <- ann[o,,drop=FALSE]
cols <- list(Module=c("Regulator"="#762a83","T3SS export"="#1b7837","Basal body"="#2166ac",
  "Hook"="#f4a582","Motor"="#b2182b","Filament"="#000000","Other"="grey70"))

png("figuras/Fig3_heatmap.png", width=195, height=225, units="mm", res=300)
pheatmap(mm2, scale="row", cluster_rows=FALSE, cluster_cols=FALSE,
  annotation_row=ann2, annotation_colors=cols, gaps_col=c(2,4),
  annotation_names_row=FALSE,
  color=colorRampPalette(rev(brewer.pal(9,"RdBu")))(60),
  fontsize_row=8.2, fontsize_col=10, fontsize=10,
  border_color="white", angle_col=0,
  main="Flagellar regulon across temperature")
dev.off()

# Version simplificada: media por temperatura, referida al valor a 25 C
prom <- t(apply(mm2, 1, function(x) c("25"=mean(x[1:2]),"30"=mean(x[3:4]),"37"=mean(x[5:6]))))
prom <- prom - prom[,1]
lin <- data.frame(gen=rep(rownames(prom),3), temp=rep(c(25,30,37), each=nrow(prom)),
                  val=as.vector(prom))
lin$tipo <- ifelse(grepl("RS04185", lin$gen), "FliO", "Other flagellar genes (n = 30)")
lin$tipo <- factor(lin$tipo, levels=c("FliO","Other flagellar genes (n = 30)"))
pl <- ggplot(lin, aes(temp, val, group=gen, colour=tipo)) +
  geom_hline(yintercept=0, colour="grey80", linewidth=.4) +
  geom_line(linewidth=1.0, alpha=.55) + geom_point(size=2.4, alpha=.75) +
  scale_colour_manual(values=c("FliO"="#1b7837",
    "Other flagellar genes (n = 30)"="#b2182b"), name=NULL) +
  scale_x_continuous(breaks=c(25,30,37), labels=c("25 \u00b0C","30 \u00b0C","37 \u00b0C"),
                     expand=expansion(mult=.09)) +
  labs(x=NULL, y=expression("Change in expression, log"[2]*" (vs 25 \u00b0C)")) +
  annotate("text", x=24.75, y=-2.4,
    label="Most of the change\noccurs between\n25 and 30 \u00b0C",
    size=4.4, fontface="bold", colour="grey25", hjust=0, lineheight=1.05) +
  theme_minimal(base_size=17) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        axis.text.x=element_text(size=19, face="bold", colour="grey10"),
        axis.text.y=element_text(size=14),
        axis.title.y=element_text(size=14, face="bold"),
        legend.position="top", legend.text=element_text(size=15),
        plot.margin=margin(14,18,10,10))
ggsave("figuras/Fig3_lineas.png", pl, width=195, height=135, units="mm", dpi=300, bg="white")

# ============ FIGURA 4: panel comparativo del genero ============
# Requiere haber ejecutado 04_panel_genero.sh, que genera panel.tsv
setwd("~/bbrna/genero")
pn <- read.delim("panel.tsv", stringsAsFactors=FALSE)
pn$lt <- map$lt[match(pn$gen, map$prot)]
sp <- c("bacilliformis","clarridgeiae","schoenbuchensis","ancashensis",
        "tribocorum","quintana","henselae")
g <- expand.grid(lt=unique(pn$lt), especie=sp, stringsAsFactors=FALSE)
g$ident <- pn$ident[match(paste(g$lt,g$especie), paste(pn$lt,pn$especie))]
g$prod <- cds$prod[match(g$lt, cds$lt)]
g$nom <- sub("^flagellar ","",g$prod)
g$nom <- sub(" domain-containing protein","",g$nom); g$nom <- substr(g$nom,1,34)
g$eti <- factor(paste0(g$lt,"  ",g$nom),
                levels=rev(unique(paste0(g$lt,"  ",g$nom)[order(g$lt)])))
# B. bacilliformis es la referencia (100 % por definicion): se omite la columna
g2 <- g[g$especie != "bacilliformis", ]
g2$especie <- factor(as.character(g2$especie), levels=sp[-1])
levels(g2$especie) <- paste0("B. ", levels(g2$especie))
nf <- nlevels(g2$eti)
exc <- g2[g2$lt %in% c("RS04185","RS05600") &
          g2$especie %in% c("B. tribocorum","B. quintana","B. henselae"), ]

p4 <- ggplot(g2, aes(especie, eti, fill=ident)) +
  geom_tile(colour="white", linewidth=.45) +
  geom_tile(data=exc, fill=NA, colour="grey10", linewidth=.75, show.legend=FALSE) +
  annotate("segment", x=3.5, xend=3.5, y=.5, yend=nf+1.9, colour="grey15", linewidth=.7) +
  annotate("segment", x=.5, xend=3.45, y=nf+1.15, yend=nf+1.15, colour="grey15", linewidth=.7) +
  annotate("segment", x=3.55, xend=6.5, y=nf+1.15, yend=nf+1.15, colour="grey15", linewidth=.7) +
  annotate("text", x=2, y=nf+1.75, label="FLAGELLUM PRESENT",
           fontface="bold", colour="grey15", size=4.1) +
  annotate("text", x=5, y=nf+1.75, label="FLAGELLUM ABSENT",
           fontface="bold", colour="grey15", size=4.1) +
  scale_fill_gradientn(colours=c("#F3B56A","#D46A55","#A53447","#68001F"),
    na.value="#D9D9D9",
    name=expression("BLASTp identity to "*italic("B. bacilliformis")*" (%)"),
    guide=guide_colourbar(direction="horizontal",
      theme=theme(legend.key.width=grid::unit(70,"mm"),
                  legend.key.height=grid::unit(4.5,"mm"),
                  legend.title.position="top", legend.text.position="bottom"))) +
  scale_y_discrete(expand=expansion(add=c(.15,2.4))) +
  labs(x=NULL, y=NULL,
    caption="Gray cells: no homolog recovered by BLASTp (E < 1e-5). Boxed cells: FliI (ATPase-domain hit) and FliO (genuinely retained).") +
  coord_cartesian(clip="off") +
  theme_minimal(base_size=11) +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(angle=30, hjust=1, vjust=1, face="italic",
                                 size=10, colour="grey15", margin=margin(t=6)),
        axis.text.y=element_text(size=8.6, colour="grey15", margin=margin(r=7)),
        axis.ticks=element_blank(),
        legend.position="top", legend.justification="left",
        legend.title=element_text(size=10.5, face="bold", colour="grey15"),
        legend.text=element_text(size=9.5, colour="grey20"),
        plot.caption=element_text(size=9, colour="grey30", hjust=0, margin=margin(t=8)),
        plot.caption.position="plot", plot.margin=margin(8,14,10,8))
ggsave("../figuras/Fig4_genero.png", p4, width=185, height=235, units="mm", dpi=300, bg="white")

# Version simplificada: recuento de genes recuperados por especie
res <- as.data.frame(table(pn$especie)); names(res) <- c("especie","n")
res <- res[res$especie != "bacilliformis", ]
res$flagelo <- ifelse(res$especie %in% c("clarridgeiae","schoenbuchensis","ancashensis"),
                      "Flagellated","Non-flagellated")
res$especie <- factor(res$especie,
  levels=c("henselae","quintana","tribocorum","ancashensis","schoenbuchensis","clarridgeiae"))
pb <- ggplot(res, aes(especie, n, fill=flagelo)) +
  geom_col(width=.68) +
  geom_text(aes(label=n), hjust=-0.35, size=6.5, fontface="bold", colour="grey15") +
  coord_flip(ylim=c(0,36)) +
  scale_fill_manual(values=c("Flagellated"="#b2182b","Non-flagellated"="grey65"), name=NULL) +
  labs(x=NULL, y="Flagellar genes recovered (of 31)") +
  theme_minimal(base_size=17) +
  theme(axis.text.y=element_text(face="italic", size=17, colour="grey10"),
        axis.text.x=element_text(size=14),
        axis.title.x=element_text(size=15, face="bold", margin=margin(t=10)),
        panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        panel.grid.major.x=element_line(colour="grey90"),
        legend.position="top", legend.text=element_text(size=15),
        plot.margin=margin(14,20,10,10))
ggsave("../figuras/Fig4_barras.png", pb, width=220, height=125, units="mm", dpi=300, bg="white")

cat("Figuras generadas en ~/bbrna/figuras/\n")
