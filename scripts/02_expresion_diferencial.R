# 02. Expresion diferencial con DESeq2.
#
# Importa las 23 cuantificaciones de salmon, ajusta dos modelos y escribe
# los cinco contrastes:
#   - modelo con las 23 muestras (solo para el control de calidad y la
#     comparacion de robustez)
#   - modelo restringido a las 12 muestras de placa/caldo de las series de
#     temperatura (Pl25, Pl30, Pl37) y de pH (pH06, pH07, pH08)
suppressMessages({library(tximport); library(DESeq2); library(ashr)})
source("scripts/rutas.R")

muestras_dir <- list.dirs(DIR_CUANT, recursive = FALSE, full.names = FALSE)
cat("Carpetas en", DIR_CUANT, ":", length(muestras_dir), "\n")

tabla <- read.delim("datos/muestras.tsv", stringsAsFactors = FALSE)
stopifnot(all(c("run", "cond") %in% names(tabla)),
          setequal(tabla$run, muestras_dir))
tabla <- tabla[match(muestras_dir, tabla$run), ]

archivos <- file.path(DIR_CUANT, muestras_dir, "quant.sf"); names(archivos) <- muestras_dir
stopifnot(all(file.exists(archivos)))

txi <- tximport(archivos, type = "salmon", txOut = TRUE, dropInfReps = TRUE)
diseno <- data.frame(row.names = muestras_dir, cond = factor(tabla$cond))
cat("\nMuestras por condicion:\n"); print(table(diseno$cond))

# ---- Modelo con las 23 muestras ----
dds <- DESeqDataSetFromTximport(txi, colData = diseno, design = ~cond)
cat("\nTranscritos totales    :", nrow(dds), "\n")
dds <- dds[rowSums(counts(dds)) > 0, ]
cat("Transcritos analizables:", nrow(dds), "\n")
dds <- DESeq(dds, quiet = TRUE)
saveRDS(dds, int("dds_23_muestras.rds"))

vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "cond", ntop = 500, returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
cat("\n*** PC1 =", pv[1], "%   PC2 =", pv[2], "% ***\n")
write.csv(pca, int("pca.csv"))
print(pca[order(pca$PC1), c("cond", "PC1", "PC2")])

# ---- Modelo restringido a 12 muestras de placa/caldo ----
en_modelo_restringido <- diseno$cond %in% c("Pl25", "Pl30", "Pl37", "pH06", "pH07", "pH08")
cat("\nMuestras en el modelo restringido:", sum(en_modelo_restringido), "\n")
dds12 <- dds[, en_modelo_restringido]
colData(dds12)$cond <- droplevels(colData(dds12)$cond)
dds12 <- DESeq(dds12, quiet = TRUE)
saveRDS(dds12, int("dds_12_muestras.rds"))

# ---- Contrastes, con encogimiento ashr ----
proteina <- sub("_[0-9]+$", "", sub(".*_cds_", "", rownames(dds)))
contraste <- function(objeto, a, b) {
  d <- as.data.frame(lfcShrink(objeto, contrast = c("cond", a, b), type = "ashr", quiet = TRUE))
  d$wp <- proteina
  d
}
ed <- function(nombre) res("expresion_diferencial", nombre)
write.csv(contraste(dds,   "Pl37", "Pl25"), ed("37C_vs_25C_modelo_23_muestras.csv"))
write.csv(contraste(dds12, "Pl37", "Pl25"), ed("37C_vs_25C.csv"))
write.csv(contraste(dds12, "Pl30", "Pl25"), ed("30C_vs_25C.csv"))
write.csv(contraste(dds12, "Pl37", "Pl30"), ed("37C_vs_30C.csv"))
write.csv(contraste(dds12, "pH06", "pH08"), ed("pH6_vs_pH8.csv"))
cat("\nCinco contrastes escritos en", res("expresion_diferencial"), "\n")
