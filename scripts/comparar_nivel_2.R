# Compara la replica desde el SRA (nivel 2) con los resultados publicados.
#
# Uso (lo llama scripts/replica_nivel_2.sh al final):
#   Rscript scripts/comparar_nivel_2.R
#
# salmon con varios hilos no es determinista, asi que desde el SRA no se espera
# una copia byte a byte, sino las mismas conclusiones. Se compara en dos niveles:
#   A. Cuantificacion: lecturas y asignaciones de cada corrida.
#   B. Cifras principales de CIFRAS_CONFIRMADAS.md, recalculadas en las dos
#      carpetas de resultados, y genes que cambian de significancia.
# Escribe la tabla B en resultados_desde_sra/comparacion_con_publicado.csv.
source("scripts/rutas.R")
options(width = 200)
PUB_C <- "cuantificacion"; SRA_C <- leer_variable("DIR_CUANT_SRA", "cuantificacion_desde_sra")
PUB_R <- "resultados";     SRA_R <- leer_variable("DIR_RES_SRA",   "resultados_desde_sra")
for (d in c(SRA_C, SRA_R)) if (!dir.exists(d)) stop("No existe ", d, ". Corre antes scripts/replica_nivel_2.sh")

# ------------------------------------------------ A. cuantificacion ---------
meta <- function(dir, run) {
  j <- readLines(file.path(dir, run, "aux_info", "meta_info.json"), warn = FALSE)
  num <- function(k) as.numeric(sub(".*: *([0-9.]+).*", "\\1", grep(paste0('"', k, '"'), j, value = TRUE)[1]))
  c(procesadas = num("num_processed"), asignadas = num("num_mapped"))
}
lecturas <- function(dir, run) {
  q <- read.delim(file.path(dir, run, "quant.sf"), stringsAsFactors = FALSE)
  setNames(q$NumReads, q$Name)
}
runs <- read.delim("datos/muestras.tsv", stringsAsFactors = FALSE)[, 1]
A <- do.call(rbind, lapply(runs, function(r) {
  mp <- meta(PUB_C, r); ms <- meta(SRA_C, r)
  lp <- lecturas(PUB_C, r); ls <- lecturas(SRA_C, r)[names(lp)]
  data.frame(run = r, procesadas_pub = mp[["procesadas"]], procesadas_sra = ms[["procesadas"]],
             asignadas_pub = mp[["asignadas"]], asignadas_sra = ms[["asignadas"]],
             dif_asignadas_pct = round(100 * (ms[["asignadas"]] - mp[["asignadas"]]) / mp[["asignadas"]], 3),
             transcritos_distintos = sum(abs(ls - lp) > 0.5),
             correlacion_log = round(cor(log2(lp + 1), log2(ls + 1)), 6))
}))
cat("=== A. CUANTIFICACION (", nrow(A), "corridas ) ===\n")
print(A, row.names = FALSE)
cat("\nLecturas procesadas identicas en", sum(A$procesadas_pub == A$procesadas_sra), "de", nrow(A), "corridas\n")
cat("Diferencia en lecturas asignadas: maxima", max(abs(A$dif_asignadas_pct)), "%\n")
cat("Correlacion minima de NumReads (log2):", min(A$correlacion_log), "\n\n")

# ------------------------------------------------ B. cifras principales -----
flag <- read.delim(file.path(PUB_R, "anotacion", "conjunto_flagelar.tsv"), header = FALSE)[, 1]
genes <- read.csv(file.path(PUB_R, "genes", "tabla_de_genes.csv"), stringsAsFactors = FALSE)
de <- function(dir, f) {        # contraste de DESeq2 con el locus de cada transcrito
  d <- read.csv(file.path(dir, "expresion_diferencial", f), row.names = 1)
  d$locus <- genes$locus[match(rownames(d), genes$transcrito)]
  d$sig <- !is.na(d$padj) & d$padj < 0.05
  d
}
leer <- function(dir, sub, f) read.csv(file.path(dir, sub, f), stringsAsFactors = FALSE)

cifras <- function(dir) {
  t  <- de(dir, "37C_vs_25C.csv"); t30 <- de(dir, "30C_vs_25C.csv"); t3730 <- de(dir, "37C_vs_30C.csv")
  ph <- de(dir, "pH6_vs_pH8.csv"); gl <- de(dir, "37C_vs_25C_modelo_23_muestras.csv")
  fl <- t$locus %in% flag
  hy <- leer(dir, "enriquecimiento", "enriquecimiento_hipergeometrico.csv")
  fg <- leer(dir, "enriquecimiento", "enriquecimiento_fgsea.csv")
  un <- leer(dir, "enriquecimiento", "universo_genes.csv")
  su <- leer(dir, "robustez", "sensibilidad_umbral_magnitud.csv")
  st <- leer(dir, "string", "enriquecimiento_string_reprimidos.csv")
  kg <- st[st$termino == "bbk02040", ]
  c("37 vs 25: significativos (padj<0,05)"          = sum(t$sig),
    "37 vs 25: significativos y |LFC|>1"            = sum(t$sig & abs(t$log2FoldChange) > 1),
    "37 vs 25: flagelares significativos"           = sum(t$sig & fl),
    "37 vs 25: flagelares reprimidos"               = sum(t$sig & fl & t$log2FoldChange < 0),
    "37 vs 25: flagelares inducidos"                = sum(t$sig & fl & t$log2FoldChange > 0),
    "Modelo de 23 muestras: significativos"         = sum(gl$sig),
    "Modelo de 23 muestras: flagelares significativos" = sum(gl$sig & gl$locus %in% flag),
    "30 vs 25: significativos"                      = sum(t30$sig),
    "30 vs 25: flagelares significativos"           = sum(t30$sig & t30$locus %in% flag),
    "37 vs 30: flagelares no significativos"        = sum(!t3730$sig & t3730$locus %in% flag),
    "pH 6 vs 8: significativos"                     = sum(ph$sig),
    "pH 6 vs 8: flagelares significativos"          = sum(ph$sig & ph$locus %in% flag),
    "Universo: inducidos"                           = un$inducidos,
    "Universo: reprimidos"                          = un$reprimidos,
    "Hipergeometrico Flagelar: reprimidos"          = hy$reprimidos[hy$categoria == "Flagelar"],
    "Hipergeometrico Flagelar: q"                   = signif(hy$q_rep[hy$categoria == "Flagelar"], 3),
    "Hipergeometrico Chaperona: q de inducidos"     = signif(hy$q_ind[hy$categoria == "Chaperona"], 3),
    "fgsea Flagelar: NES"                           = round(fg$NES[fg$pathway == "Flagelar"], 3),
    "fgsea Flagelar: q"                             = signif(fg$padj[fg$pathway == "Flagelar"], 3),
    "Umbral formal lfcThreshold=1: significativos"  = su$n_significativos[grepl("formal", su$criterio)],
    "Umbral formal lfcThreshold=1: flagelares"      = su$n_flagelares[grepl("formal", su$criterio)],
    "STRING KEGG bbk02040: reprimidos"              = kg$genes_conjunto,
    "STRING KEGG bbk02040: q"                       = signif(kg$q, 3))
}
cp <- cifras(PUB_R); cs <- cifras(SRA_R)[names(cp)]
# Conteos: se comparan exactos. Valores q y NES: diferencia relativa en %.
es_conteo <- !grepl(": (q|NES)$|: q de ", names(cp))
B <- data.frame(cifra = names(cp), publicado = unname(cp), desde_sra = unname(cs),
                diferencia = ifelse(es_conteo, unname(cs - cp), round(100 * unname((cs - cp) / cp), 1)),
                unidad = ifelse(es_conteo, "genes", "% relativo"))
escribir <- function(x) vapply(x, function(v) if (v != 0 && abs(v) < 1e-3) formatC(v, format = "e", digits = 2) else format(v), "")
cat("=== B. CIFRAS PRINCIPALES ===\n")
print(data.frame(cifra = B$cifra, publicado = escribir(B$publicado), desde_sra = escribir(B$desde_sra),
                 diferencia = paste(B$diferencia, B$unidad)), row.names = FALSE, right = FALSE)
write.csv(B, file.path(SRA_R, "comparacion_con_publicado.csv"), row.names = FALSE)

cat("\n=== GENES QUE CAMBIAN DE SIGNIFICANCIA (padj<0,05) ===\n")
for (f in c("37C_vs_25C.csv", "30C_vs_25C.csv", "37C_vs_30C.csv", "pH6_vs_pH8.csv", "37C_vs_25C_modelo_23_muestras.csv")) {
  p <- de(PUB_R, f); s <- de(SRA_R, f)[rownames(p), ]
  cambian <- p$locus[p$sig != s$sig]
  ok <- !is.na(p$log2FoldChange) & !is.na(s$log2FoldChange)
  cat(sprintf("%-36s %3d genes cambian (%d flagelares); correlacion de log2FC %.5f; diferencia maxima %.3f\n",
              f, length(cambian), sum(cambian %in% flag), cor(p$log2FoldChange[ok], s$log2FoldChange[ok]),
              max(abs(p$log2FoldChange[ok] - s$log2FoldChange[ok]))))
  if (length(cambian) && length(cambian) <= 15)
    cat("     ", paste0(sub("BARBAKC583_", "", cambian), " (padj ", signif(p$padj[p$sig != s$sig], 2), " -> ",
                        signif(s$padj[p$sig != s$sig], 2), ")", collapse = "; "), "\n")
}

cat("\nRESULTADO: conteos iguales en", sum(B$diferencia[es_conteo] == 0), "de", sum(es_conteo),
    "; valores q y NES con diferencia relativa maxima de", max(abs(B$diferencia[!es_conteo])), "%.\n")
cat("Desde el SRA no se espera identidad exacta (salmon con varios hilos no es determinista):\n",
    "lo que importa es que no cambie ninguna conclusion. Revisa las tablas A y B.\n")
