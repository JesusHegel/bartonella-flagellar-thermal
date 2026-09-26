# 08. Tabla de genes: que gen es cada transcrito y cuanto se expresa.
#
# Una fila por cada CDS del genoma (1 216 secuencias de cds_from_genomic.fna):
#   - identidad: locus actual y antiguo, nombre del gen, producto, proteina,
#     coordenadas, pseudogen, copias con secuencia identica
#   - si entro en la cuantificacion, si es analizable y si forma parte del
#     universo del enriquecimiento; categoria funcional
#   - TPM medio por condicion (salmon; corregido por longitud, permite comparar
#     genes distintos entre si)
#   - log2 del cambio y padj de los cuatro contrastes del modelo restringido
#
# Tambien escribe los 20 genes mas inducidos y los 20 mas reprimidos a 37 C
# frente a 25 C, con los conteos normalizados de cada replica.
suppressMessages({library(tximport); library(DESeq2)})
source("scripts/rutas.R")
source("scripts/categorias.R")

# ---- Identidad de cada CDS, desde las cabeceras del FASTA del NCBI ----
fa  <- readLines("datos/referencia/cds_from_genomic.fna")
cab <- grep("^>", fa)
seqs <- vapply(seq_along(cab), function(i) {
  fin <- if (i < length(cab)) cab[i + 1] - 1 else length(fa)
  paste(fa[(cab[i] + 1):fin], collapse = "")
}, character(1))
atributo <- function(clave) {
  m <- regmatches(fa[cab], regexpr(paste0("\\[", clave, "=[^]]*\\]"), fa[cab]))
  v <- rep(NA_character_, length(cab))
  v[grepl(paste0("\\[", clave, "="), fa[cab])] <- sub(paste0("^\\[", clave, "=(.*)\\]$"), "\\1", m)
  v
}
g <- data.frame(
  transcrito = sub("^>(\\S+).*", "\\1", fa[cab]),
  locus      = atributo("locus_tag"),
  gen        = atributo("gene"),
  producto   = atributo("protein"),
  proteina   = atributo("protein_id"),
  ubicacion  = atributo("location"),
  pseudogen  = grepl("\\[pseudo=true\\]", fa[cab]),
  longitud_nt = nchar(seqs),
  stringsAsFactors = FALSE)
g$hebra <- ifelse(grepl("^complement", g$ubicacion), "-", "+")

# Locus antiguo (usado por STRING y UniProt), desde los genes del GFF
lineas <- strsplit(grep("^#", readLines("datos/referencia/genomic.gff"), value = TRUE, invert = TRUE), "\t")
gff <- vapply(Filter(function(x) length(x) >= 9 && x[3] == "gene", lineas), `[`, character(1), 9)
extraer <- function(x, clave) {       # clave exacta: "locus_tag" no confunde con "old_locus_tag"
  patron <- paste0("(^|;)", clave, "=([^;]+)")
  ifelse(grepl(patron, x), sub(paste0(".*", patron, ".*"), "\\2", x), NA)
}
g$locus_antiguo <- extraer(gff, "old_locus_tag")[match(g$locus, extraer(gff, "locus_tag"))]

# Copias con secuencia identica (salmon conserva una sola y le asigna las lecturas)
grupos <- split(g$locus, seqs)
g$copias_identicas <- vapply(seq_len(nrow(g)), function(i) {
  otros <- setdiff(grupos[[seqs[i]]], g$locus[i])
  if (length(otros)) paste(otros, collapse = ";") else ""
}, character(1))

# ---- Cuantificacion y analisis ----
dds   <- readRDS(int("dds_23_muestras.rds"))
dds12 <- readRDS(int("dds_12_muestras.rds"))
muestras <- colnames(dds12)
archivos <- file.path(DIR_CUANT, muestras, "quant.sf"); names(archivos) <- muestras
txi <- tximport(archivos, type = "salmon", txOut = TRUE, dropInfReps = TRUE)
g$en_cuantificacion <- g$transcrito %in% rownames(txi$abundance)
g$analizable        <- g$transcrito %in% rownames(dds)

# Universo del enriquecimiento: misma regla que el script 03
ed  <- function(nombre) read.csv(res("expresion_diferencial", nombre), row.names = 1)
r37 <- ed("37C_vs_25C.csv")
cds <- read.delim(res("anotacion", "anotacion_cds.tsv"), header = FALSE, stringsAsFactors = FALSE)
names(cds) <- c("lt", "ini", "fin", "hebra", "prod")
map <- read.delim(res("anotacion", "locus_a_proteina.tsv"), header = FALSE, stringsAsFactors = FALSE)
names(map) <- c("lt", "prot")
cds$wp <- map$prot[match(cds$lt, map$lt)]
u <- cds[!duplicated(cds$wp) & !is.na(cds$wp) & cds$wp != "", ]
u$i <- match(u$wp, r37$wp)
u <- u[!is.na(u$i) & !is.na(r37$padj[u$i]), ]
g$en_universo <- g$locus %in% u$lt

fl <- read.delim(res("anotacion", "conjunto_flagelar.tsv"), header = FALSE)[, 1]
cats <- definir_categorias(g$locus, g$producto, fl)
g$categoria <- vapply(g$locus, function(l) paste(names(cats)[vapply(cats, function(v) l %in% v, logical(1))],
                                                   collapse = ";"), character(1), USE.NAMES = FALSE)
g$conjunto_flagelar <- g$locus %in% fl

# TPM medio por condicion (modelo restringido)
cond <- as.character(colData(dds12)$cond)
for (k in c("Pl25", "Pl30", "Pl37", "pH06", "pH07", "pH08")) {
  tpm <- rowMeans(txi$abundance[, cond == k, drop = FALSE])
  g[[paste0("TPM_", k)]] <- round(tpm[match(g$transcrito, names(tpm))], 1)
}

# Resultados de los cuatro contrastes del modelo restringido
for (ct in list(c("37C_vs_25C", "37vs25"), c("30C_vs_25C", "30vs25"),
                c("37C_vs_30C", "37vs30"), c("pH6_vs_pH8", "pH6vspH8"))) {
  r <- ed(paste0(ct[1], ".csv")); i <- match(g$transcrito, rownames(r))
  g[[paste0("log2FC_", ct[2])]] <- round(r$log2FoldChange[i], 3)
  g[[paste0("padj_", ct[2])]]   <- signif(r$padj[i], 3)
}
write.csv(g[, c("transcrito", "locus", "locus_antiguo", "gen", "producto", "proteina", "ubicacion",
                "hebra", "longitud_nt", "pseudogen", "copias_identicas", "en_cuantificacion",
                "analizable", "en_universo", "conjunto_flagelar", "categoria",
                grep("^TPM_|^log2FC_|^padj_", names(g), value = TRUE))],
          res("genes", "tabla_de_genes.csv"), row.names = FALSE)

cat("CDS en el genoma                 :", nrow(g), "\n")
cat("  pseudogenes                    :", sum(g$pseudogen), "\n")
cat("  con copia de secuencia identica:", sum(g$copias_identicas != ""), "\n")
cat("En la cuantificacion de salmon   :", sum(g$en_cuantificacion), "\n")
cat("Analizables en DESeq2            :", sum(g$analizable), "\n")
cat("Universo del enriquecimiento     :", sum(g$en_universo), "\n")
cat("Con nombre de gen                :", sum(!is.na(g$gen)), "\n")
cat("Con locus antiguo                :", sum(!is.na(g$locus_antiguo)), "\n")

# ---- Los 20 mas inducidos y los 20 mas reprimidos, 37 C frente a 25 C ----
nc <- counts(dds12, normalized = TRUE)
orden <- order(cond, colnames(dds12))
sel <- orden[cond[orden] %in% c("Pl25", "Pl30", "Pl37")]
sig <- g[g$analizable & !is.na(g$padj_37vs25) & g$padj_37vs25 < 0.05, ]
arriba <- head(sig[order(-sig$log2FC_37vs25), ], 20)
abajo  <- head(sig[order(sig$log2FC_37vs25), ], 20)
top <- rbind(cbind(direccion = "inducido", posicion = 1:nrow(arriba), arriba),
             cbind(direccion = "reprimido", posicion = 1:nrow(abajo), abajo))
top$baseMean <- round(r37$baseMean[match(top$transcrito, rownames(r37))], 1)
conteos <- round(nc[match(top$transcrito, rownames(nc)), sel, drop = FALSE])
colnames(conteos) <- paste0(cond[sel], "_", colnames(dds12)[sel])
top <- cbind(top[, c("direccion", "posicion", "locus", "gen", "producto", "categoria",
                     "baseMean", "log2FC_37vs25", "padj_37vs25")], conteos)
write.csv(top, res("genes", "top20_37C_vs_25C.csv"), row.names = FALSE)

cat("\n--- 20 mas inducidos (37 C frente a 25 C) ---\n")
print(top[top$direccion == "inducido", c("posicion", "locus", "gen", "producto", "log2FC_37vs25", "padj_37vs25")], row.names = FALSE)
cat("\n--- 20 mas reprimidos ---\n")
print(top[top$direccion == "reprimido", c("posicion", "locus", "gen", "producto", "log2FC_37vs25", "padj_37vs25")], row.names = FALSE)
cat("\nFlagelares entre los 20 mas reprimidos:", sum(top$direccion == "reprimido" & grepl("Flagelar", top$categoria)), "\n")
