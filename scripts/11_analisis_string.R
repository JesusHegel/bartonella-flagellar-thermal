# 11. Analisis con STRING: anotacion independiente y red de asociaciones.
#
# Usa la copia de STRING guardada en datos/string/ (script 10). Todo el calculo
# es local y con el mismo universo que el script 03 (proteinas deduplicadas con
# resultado en DESeq2), para que sea comparable con nuestro enriquecimiento.
#
#   1. Enriquecimiento hipergeometrico de los terminos de STRING (KEGG, GO,
#      palabras clave, dominios) en los genes reprimidos e inducidos a 37 C
#      frente a 25 C (padj < 0,05). Correccion BH dentro de cada categoria.
#   2. Concordancia entre nuestro conjunto flagelar (palabras clave) y la via
#      KEGG de ensamblaje flagelar (anotacion independiente).
#   3. Modulo flagelar en la red: aristas entre los genes flagelares frente a
#      10 000 conjuntos aleatorios del mismo tamano (semilla 42).
#   4. Red de los genes con cambio mayor a 2 veces, para la figura.
source("scripts/rutas.R")
D <- "datos/string"
if (!file.exists(file.path(D, "version.tsv"))) stop("Falta datos/string/. Ejecuta antes el script 10.")
leer <- function(f) read.delim(file.path(D, f), quote = "", stringsAsFactors = FALSE, check.names = FALSE)
version <- leer("version.tsv"); ids <- leer("identificadores.tsv")
red <- leer("red.tsv"); ann <- leer("anotacion_funcional.tsv")
cat("STRING", version[1, 1], "descargado el", readLines(file.path(D, "fecha_descarga.txt")), "\n")

g <- read.csv(res("genes", "tabla_de_genes.csv"), stringsAsFactors = FALSE)
g$string_id <- ids$stringId[match(g$locus_antiguo, ids$queryItem)]
U <- g[g$en_universo & !is.na(g$string_id), ]
N <- nrow(U)
rep <- U$string_id[U$padj_37vs25 < 0.05 & U$log2FC_37vs25 < 0]
ind <- U$string_id[U$padj_37vs25 < 0.05 & U$log2FC_37vs25 > 0]
cat("Universo:", sum(g$en_universo), "proteinas;", N, "con identificador STRING\n")
cat("Reprimidos:", length(rep), " Inducidos:", length(ind), "\n\n")

# ---- 1. Enriquecimiento local de los terminos de STRING ----
completar <- function(x) ifelse(grepl("^360095\\.", x), x, paste0("360095.", x))
ann$miembros <- lapply(strsplit(ann$inputGenes, ","), function(v) intersect(completar(trimws(v)), U$string_id))
ann$K <- lengths(ann$miembros)
if (!any(ann$K > 0)) stop("Ningun termino de STRING coincide con el universo: revisar el formato de inputGenes")
cats <- c("KEGG", "Process", "Function", "Component", "Keyword", "Pfam", "InterPro")
ann <- ann[ann$category %in% cats & ann$K >= 3, ]
nombre <- setNames(ifelse(is.na(U$gen), sub("BARBAKC583_", "", U$locus), U$gen), U$string_id)

enriquecer <- function(conjunto, etiqueta) {
  n <- length(conjunto)
  r <- do.call(rbind, lapply(seq_len(nrow(ann)), function(i) {
    m <- ann$miembros[[i]]; k <- sum(m %in% conjunto); K <- length(m)
    data.frame(categoria = ann$category[i], termino = ann$term[i], descripcion = ann$description[i],
               genes_universo = K, genes_conjunto = k, esperados = round(K * n / N, 1),
               p = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
               genes = paste(sort(nombre[m[m %in% conjunto]]), collapse = ";"))
  }))
  r$q <- ave(r$p, r$categoria, FUN = function(p) p.adjust(p, "BH"))
  r <- r[order(r$q, r$p), ]
  write.csv(r, res("string", paste0("enriquecimiento_string_", etiqueta, ".csv")), row.names = FALSE)
  cat("--- Terminos con q < 0,05 en genes", etiqueta, "---\n")
  print(head(r[r$q < 0.05, c("categoria", "termino", "descripcion", "genes_universo",
                              "genes_conjunto", "esperados", "q")], 25), row.names = FALSE)
  cat("\n")
  r
}
er <- enriquecer(rep, "reprimidos")
ei <- enriquecer(ind, "inducidos")

# ---- 2. Conjunto flagelar frente a la via KEGG de ensamblaje flagelar ----
flag <- U$string_id[U$conjunto_flagelar]
kf <- ann[ann$category == "KEGG" & grepl("flagellar assembly", ann$description, ignore.case = TRUE), ]
if (nrow(kf)) {
  mk <- kf$miembros[[1]]
  conc <- data.frame(termino_kegg = kf$term[1], genes_kegg = length(mk), genes_flagelares = length(flag),
                     en_ambos = length(intersect(mk, flag)),
                     solo_kegg = paste(sort(nombre[setdiff(mk, flag)]), collapse = ";"),
                     solo_flagelares = paste(sort(nombre[setdiff(flag, mk)]), collapse = ";"))
  write.csv(conc, res("string", "concordancia_flagelar_kegg.csv"), row.names = FALSE)
  cat("--- Conjunto flagelar frente a KEGG", kf$term[1], "---\n"); print(t(conc)); cat("\n")
} else cat("AVISO: no hay termino KEGG de ensamblaje flagelar en la anotacion descargada\n\n")

# ---- 3. Modulo flagelar en la red (puntaje combinado >= 0,4) ----
red$a <- pmin(red$stringId_A, red$stringId_B); red$b <- pmax(red$stringId_A, red$stringId_B)
red <- red[!duplicated(red[, c("a", "b")]) & red$a %in% U$string_id & red$b %in% U$string_id & red$score >= 0.4, ]
aristas_dentro <- function(v) sum(red$a %in% v & red$b %in% v)
obs <- aristas_dentro(flag)
set.seed(42)
azar <- replicate(10000, aristas_dentro(sample(U$string_id, length(flag))))
mod <- data.frame(genes_flagelares = length(flag), aristas_observadas = obs,
                  aristas_azar_media = round(mean(azar), 1), aristas_azar_max = max(azar),
                  p_empirico = (sum(azar >= obs) + 1) / (length(azar) + 1),
                  aristas_red_universo = nrow(red), proteinas_universo = N)
write.csv(mod, res("string", "modulo_flagelar_red.csv"), row.names = FALSE)
cat("--- Modulo flagelar en la red ---\n"); print(t(mod)); cat("\n")

# ---- 4. Red de los genes con cambio mayor a 2 veces (para la figura) ----
c2 <- U[U$padj_37vs25 < 0.05 & abs(U$log2FC_37vs25) > 1, ]
e2 <- red[red$a %in% c2$string_id & red$b %in% c2$string_id, c("a", "b", "score")]
e2$gen_a <- nombre[e2$a]; e2$gen_b <- nombre[e2$b]
write.csv(e2, res("string", "red_cambio_2x_aristas.csv"), row.names = FALSE)
nodos <- data.frame(string_id = c2$string_id, gen = nombre[c2$string_id], locus = c2$locus,
                    producto = c2$producto, log2FC_37vs25 = c2$log2FC_37vs25,
                    conjunto_flagelar = c2$conjunto_flagelar,
                    grado = vapply(c2$string_id, function(s) sum(e2$a == s | e2$b == s), numeric(1)))
write.csv(nodos, res("string", "red_cambio_2x_nodos.csv"), row.names = FALSE)
cat("Red de genes con cambio > 2 veces:", nrow(nodos), "genes,", nrow(e2), "aristas;",
    sum(nodos$grado == 0), "genes sin ninguna arista\n")
