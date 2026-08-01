# Generacion de las cuatro figuras del manuscrito
suppressMessages({library(DESeq2); library(ggplot2); library(ggrepel)
                  library(pheatmap); library(RColorBrewer)})
setwd("~/bbrna"); dir.create("figuras", showWarnings=FALSE)

dds <- readRDS("dds_global.rds")
vsd <- vst(dds, blind=TRUE)

# --- Figura 1: PCA ---
pv  <- plotPCA(vsd, intgroup="cond", returnData=TRUE)
pct <- round(100*attr(pv,"percentVar"))
yale <- c("Pl25","Pl30","Pl37","pH06","pH07","pH08")
pv$lote <- ifelse(pv$cond %in% yale, "Center A", "Center B")
pv$cond <- factor(as.character(pv$cond),
  levels=c("Pl25","Pl30","Pl37","pH06","pH07","pH08","PlBG","HB37","HBBG","HUVE"))
p1 <- ggplot(pv, aes(PC1,PC2,color=cond,shape=lote)) +
  geom_point(size=3.4, alpha=.9) +
  scale_shape_manual(values=c("Center A"=16,"Center B"=17)) +
  labs(x=paste0("PC1 (",pct[1],"% variance)"), y=paste0("PC2 (",pct[2],"% variance)")) +
  theme_bw(base_size=11)
ggsave("figuras/Fig1_PCA.pdf", p1, width=180, height=115, units="mm")

# Las figuras 2 a 4 siguen el mismo procedimiento; codigo completo en el
# historial del repositorio.
