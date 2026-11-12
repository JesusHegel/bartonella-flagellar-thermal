# -*- coding: utf-8 -*-
# Contraste placa frente a liquido a temperatura y pH iguales.
#
# Pl30 y pH07 comparten temperatura (30 C), pH (7,4), centro de secuenciacion
# y tiempo de traslado. Difieren solo en la matriz: agar frente a caldo. El
# script 14 mostro que el conjunto flagelar queda en el percentil ~19 en Pl30
# y ~52 en pH07, lo que sugiere que parte de la represion atribuida a la
# temperatura puede deberse al soporte. Aqui se mide con DESeq2.
#
#   Rscript scripts/16_solido_vs_liquido.R
#
# Salida: resultados/condiciones/solido_vs_liquido.csv
suppressPackageStartupMessages({ library(DESeq2); library(ashr) })

source("scripts/rutas.R")
r <- function(...) {
  p <- file.path(...)
  for (nd in list(c("resultados", DIR_RES), c("intermedios", DIR_INT),
                  c("figuras", DIR_FIG))) {
    if (grepl(paste0("^", nd[1], "(/|$)"), p))
      return(file.path(nd[2], sub(paste0("^", nd[1], "/?"), "", p)))
  }
  p
}

ARCH <- r("intermedios", "dds_12_muestras.rds")
if (!file.exists(ARCH)) {
  cat("No se encuentra", ARCH, "\nContenido de intermedios/:\n")
  print(list.files(r("intermedios"))); quit(status = 1)
}
dds <- readRDS(ARCH)
cat("Muestras:", ncol(dds), " transcritos:", nrow(dds), "\n")
cat("Niveles de cond:", paste(levels(dds$cond), collapse = ", "), "\n")
if (!all(c("Pl30", "pH07") %in% levels(dds$cond))) {
  cat("Faltan niveles Pl30 o pH07\n"); quit(status = 1)
}

# ---- conjunto flagelar: locus -> proteina -> transcrito ------------------
l2p <- read.delim(r("resultados/anotacion/locus_a_proteina.tsv"),
                  header = FALSE, stringsAsFactors = FALSE)
flag_loci <- read.delim(r("resultados/anotacion/conjunto_flagelar.tsv"),
                        header = FALSE, stringsAsFactors = FALSE)[[1]]
flag_prot <- l2p[[2]][l2p[[1]] %in% flag_loci]

prot_de <- function(x) {
  vapply(strsplit(x, "_cds_", fixed = TRUE), function(q) {
    if (length(q) != 2) return(NA_character_)
    paste(strsplit(q[2], "_", fixed = TRUE)[[1]][1:2], collapse = "_")
  }, character(1))
}
prot_fila <- prot_de(rownames(dds))
es_flag <- !is.na(prot_fila) & prot_fila %in% flag_prot
cat("Genes flagelares por anotacion:", length(flag_loci),
    " transcritos emparejados:", sum(es_flag), "\n\n")
if (sum(es_flag) < 25) { cat("Emparejamiento insuficiente\n"); quit(status = 1) }

# ---- contraste ----------------------------------------------------------
# Positivo = mas alto en liquido (pH07) que en placa (Pl30).
res <- lfcShrink(dds, contrast = c("cond", "pH07", "Pl30"), type = "ashr",
                 quiet = TRUE)
res <- as.data.frame(res)
res$flagelar <- es_flag
res$locus <- l2p[[1]][match(prot_fila, l2p[[2]])]
con_p <- res[!is.na(res$padj), ]
cat("Transcritos con valor q:", nrow(con_p), "\n\n")

PADJ <- 0.05
sube <- con_p$padj < PADJ & con_p$log2FoldChange > 0
baja <- con_p$padj < PADJ & con_p$log2FoldChange < 0
cat("Mas alto en liquido :", sum(sube), "transcritos\n")
cat("Mas alto en placa   :", sum(baja), "transcritos\n\n")

# ---- enriquecimiento hipergeometrico del conjunto flagelar --------------
hiper <- function(sel, etiq) {
  N <- nrow(con_p); K <- sum(con_p$flagelar)
  n <- sum(sel);    k <- sum(sel & con_p$flagelar)
  p <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  cat(sprintf("%-22s %2d de %2d flagelares en %4d de %4d  p = %.4g\n",
              etiq, k, K, n, N, p))
  invisible(c(k = k, K = K, n = n, N = N, p = p))
}
cat("=== ENRIQUECIMIENTO DEL CONJUNTO FLAGELAR ===\n")
h_sube <- hiper(sube, "mas alto en liquido")
h_baja <- hiper(baja, "mas alto en placa")

# ---- distribucion completa, sin umbral ----------------------------------
lf <- con_p$log2FoldChange[con_p$flagelar]
lr <- con_p$log2FoldChange[!con_p$flagelar]
w <- wilcox.test(lf, lr, alternative = "greater")
cat(sprintf("\nlog2FC mediano  flagelares: %+.3f   resto: %+.3f\n",
            median(lf), median(lr)))
cat(sprintf("Wilcoxon una cola (flagelares mas altos): p = %.4g\n", w$p.value))

# ---- los 31, uno por uno ------------------------------------------------
fl <- con_p[con_p$flagelar, c("locus", "log2FoldChange", "padj")]
fl <- fl[order(-fl$log2FoldChange), ]
cat("\n=== CONJUNTO FLAGELAR, LIQUIDO FRENTE A PLACA ===\n")
cat(sprintf("%-22s %8s %11s\n", "locus", "log2FC", "padj"))
for (i in seq_len(nrow(fl)))
  cat(sprintf("%-22s %+8.3f %11.3g\n", fl$locus[i],
              fl$log2FoldChange[i], fl$padj[i]))

dir.create(r("resultados/condiciones"), showWarnings = FALSE, recursive = TRUE)
sal <- r("resultados/condiciones/solido_vs_liquido.csv")
write.csv(con_p[, c("locus", "log2FoldChange", "lfcSE", "pvalue", "padj", "flagelar")],
          sal, row.names = FALSE)
cat("\nTabla:", sal, "\n")
