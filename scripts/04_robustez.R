# 04. Robustez y perfil termico de los 31 genes flagelares.
#
#   - conteo de genes significativos bajo cada definicion de umbral
#   - sensibilidad al umbral de magnitud (post hoc frente a lfcThreshold formal)
#   - modelo de 23 muestras frente a modelo restringido de 12, gen por gen
#   - curva 25 / 30 / 37 C de los 31 genes flagelares
#   - especificidad: respuesta a temperatura frente a respuesta a pH
#   - concordancia entre dos particiones independientes de las replicas
#   - medias de conteos normalizados a 25, 30 y 37 C de todos los genes
#     (sirven para comparar un mismo gen entre temperaturas; NO para comparar
#     genes distintos entre si, porque no estan corregidas por longitud)
suppressMessages({library(DESeq2)})
source("scripts/rutas.R")

dds   <- readRDS(int("dds_23_muestras.rds"))
dds12 <- readRDS(int("dds_12_muestras.rds"))
ed <- function(nombre) read.csv(res("expresion_diferencial", nombre), row.names = 1)
r_global <- ed("37C_vs_25C_modelo_23_muestras.csv")
r37_25   <- ed("37C_vs_25C.csv")
r30_25   <- ed("30C_vs_25C.csv")
r37_30   <- ed("37C_vs_30C.csv")
r_ph     <- ed("pH6_vs_pH8.csv")
wp <- r_global$wp

map <- read.delim(res("anotacion", "locus_a_proteina.tsv"), header = FALSE); names(map) <- c("lt", "prot")
fl  <- read.delim(res("anotacion", "conjunto_flagelar.tsv"), header = FALSE)
names(fl) <- c("lt", "ini", "fin", "hebra", "prod")
g <- match(map$prot[match(fl$lt, map$lt)], wp)       # fila de cada gen flagelar
cat("Flagelares localizados en la matriz:", sum(!is.na(g)), "de", nrow(fl), "\n\n")

# ---- Conteos bajo cada definicion ----
contar <- function(r, etiqueta) data.frame(definicion = etiqueta,
  genoma_padj05      = sum(r$padj < 0.05, na.rm = TRUE),
  genoma_padj05_lfc1 = sum(r$padj < 0.05 & abs(r$log2FoldChange) > 1, na.rm = TRUE),
  flag_padj05        = sum(r$padj[g] < 0.05, na.rm = TRUE),
  flag_padj05_lfc1   = sum(r$padj[g] < 0.05 & abs(r$log2FoldChange[g]) > 1, na.rm = TRUE))
tb <- rbind(contar(r_global, "global 37vs25"), contar(r37_25, "restringido 37vs25"),
            contar(r30_25, "restringido 30vs25"), contar(r_ph, "restringido pH06vspH08"))
cat("--- CONTEOS BAJO CADA DEFINICION ---\n"); print(tb)
write.csv(tb, res("robustez", "genes_significativos_por_umbral.csv"), row.names = FALSE)

# ---- Sensibilidad al umbral de magnitud ----
fr <- as.data.frame(results(dds12, contrast = c("cond", "Pl37", "Pl25"), lfcThreshold = 1))
est <- data.frame(criterio = c("post hoc padj<0.05 & |LFC|>1", "lfcThreshold=1 formal"),
  n_significativos = c(sum(r37_25$padj < 0.05 & abs(r37_25$log2FoldChange) > 1, na.rm = TRUE),
                       sum(fr$padj < 0.05, na.rm = TRUE)),
  n_flagelares = c(sum(r37_25$padj[g] < 0.05 & abs(r37_25$log2FoldChange[g]) > 1, na.rm = TRUE),
                   sum(fr$padj[g] < 0.05, na.rm = TRUE)))
cat("\n--- SENSIBILIDAD AL UMBRAL ---\n"); print(est)
write.csv(est, res("robustez", "sensibilidad_umbral_magnitud.csv"), row.names = FALSE)

# ---- Modelo de 23 muestras frente a modelo restringido ----
mg <- data.frame(locus = fl$lt, producto = fl$prod,
  lfc_global = round(r_global$log2FoldChange[g], 2), padj_global = signif(r_global$padj[g], 2),
  lfc_restr  = round(r37_25$log2FoldChange[g], 2),   padj_restr  = signif(r37_25$padj[g], 2))
mg$dif <- round(mg$lfc_global - mg$lfc_restr, 2); mg <- mg[order(mg$lfc_restr), ]
cat("\n--- GLOBAL vs RESTRINGIDO ---\n"); print(mg[, c("locus", "lfc_global", "lfc_restr", "dif")])
cat("\n*** Diferencia maxima entre modelos:", max(abs(mg$dif), na.rm = TRUE), "***\n")
write.csv(mg, res("robustez", "modelo_global_vs_restringido.csv"), row.names = FALSE)

# ---- Curva 25 / 30 / 37 C ----
c31 <- data.frame(locus = fl$lt, producto = fl$prod,
  L30v25 = round(r30_25$log2FoldChange[g], 2), p30v25 = signif(r30_25$padj[g], 2),
  L37v25 = round(r37_25$log2FoldChange[g], 2), p37v25 = signif(r37_25$padj[g], 2),
  L37v30 = round(r37_30$log2FoldChange[g], 2), p37v30 = signif(r37_30$padj[g], 2))
c31 <- c31[order(c31$L37v25), ]
cat("\n--- CURVA 25 / 30 / 37 ---\n"); print(c31[, c("locus", "L30v25", "L37v25", "L37v30", "p37v30")])
write.csv(c31, res("robustez", "flagelares_curva_termica.csv"), row.names = FALSE)

# ---- Especificidad: temperatura frente a pH ----
cat("\n--- ESPECIFICIDAD: temperatura vs pH ---\n")
ef <- data.frame(contraste = c("temperatura 37vs25", "pH 06vs08"),
  flagelares_sig = c(sum(r37_25$padj[g] < 0.05, na.rm = TRUE), sum(r_ph$padj[g] < 0.05, na.rm = TRUE)),
  genoma_sig     = c(sum(r37_25$padj < 0.05, na.rm = TRUE),    sum(r_ph$padj < 0.05, na.rm = TRUE)))
print(ef); write.csv(ef, res("robustez", "especificidad_temperatura_vs_ph.csv"), row.names = FALSE)

# ---- Particiones independientes de las replicas ----
nc <- counts(dds, normalized = TRUE); cn <- as.character(colData(dds)$cond)
i25 <- which(cn == "Pl25"); i37 <- which(cn == "Pl37")
cat("\nReplicas Pl25:", length(i25), " Pl30:", sum(cn == "Pl30"), " Pl37:", length(i37), "\n")
if (length(i25) >= 2 && length(i37) >= 2) {
  lr <- function(a, b) log2((nc[, b] + 1) / (nc[, a] + 1))
  A <- lr(i25[1], i37[1]); B <- lr(i25[2], i37[2])
  pt <- data.frame(locus = fl$lt, particionA = round(A[g], 2), particionB = round(B[g], 2))
  pt$misma_direccion <- sign(pt$particionA) == sign(pt$particionB)
  cat("\n--- PARTICIONES INDEPENDIENTES ---\n")
  cat("Misma direccion:", sum(pt$misma_direccion, na.rm = TRUE), "de", nrow(pt), "\n")
  cat("Correlacion genoma completo:", round(cor(A, B, use = "complete.obs"), 3), "\n")
  write.csv(pt, res("robustez", "concordancia_particiones.csv"), row.names = FALSE)
}

# ---- Medias de conteos normalizados por temperatura (todos los genes) ----
cds <- read.delim(res("anotacion", "anotacion_cds.tsv"), header = FALSE)
names(cds) <- c("lt", "ini", "fin", "hebra", "prod")
medias <- sapply(c("Pl25", "Pl30", "Pl37"), function(k) round(rowMeans(nc[, cn == k, drop = FALSE])))
tm <- data.frame(wp = wp, lt = map$lt[match(wp, map$prot)], row.names = NULL)
tm$prod <- cds$prod[match(tm$lt, cds$lt)]
tm <- cbind(tm, medias, row.names = NULL)
write.csv(tm, res("robustez", "medias_normalizadas_25_30_37C.csv"), row.names = FALSE)

cat("\n--- Genes citados en CIFRAS, seccion 8 (medias de conteos normalizados) ---\n")
citados <- c("BARBAKC583_RS05045", "BARBAKC583_RS05735", "BARBAKC583_RS01390", "BARBAKC583_RS06465",
             "BARBAKC583_RS02960", "BARBAKC583_RS05740", "BARBAKC583_RS01060", "BARBAKC583_RS02640",
             "BARBAKC583_RS05920")
print(tm[match(citados, tm$lt), c("lt", "prod", "Pl25", "Pl30", "Pl37")], row.names = FALSE)
