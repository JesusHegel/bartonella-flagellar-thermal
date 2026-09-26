# 12. Figuras del manuscrito.
#
# Todas se generan desde las tablas de resultados/ y los objetos de intermedios/;
# ninguna se edita a mano. Si hace falta cambiar una figura, se cambia este script.
# Cada figura se guarda en PDF (vectorial, para la revista) y en PNG a 300 ppp.
#
#   figura_1  Flujo de trabajo
#   figura_2  Muestras: PCA y lecturas asignadas a genes bacterianos
#   figura_3  Volcanes de temperatura y de pH, con los 31 genes flagelares
#   figura_4  Genes flagelares y chaperonas: mapa de calor por replica y curva termica
#   figura_5  Red STRING de los genes con cambio mayor a 2 veces
#   figura_6  Ortologos reciprocos de los 31 genes flagelares en el genero
#
# Colores: flagelar naranja y chaperona aqua (par validado para daltonismo),
# el resto en gris; divergente azul-gris-rojo solo para puntuaciones Z.
suppressMessages({library(DESeq2); library(ggplot2); library(igraph); library(ggrepel); library(grid)})
source("scripts/rutas.R")
set.seed(42)

# ---------------------------------------------------------------- estilo ----
FLAG  <- "#eb6834"; CHAP <- "#1baf7a"; GRIS <- "#bdbcb8"; GRIS_O <- "#6f6e69"
AZUL  <- "#2a78d6"; ROJO <- "#e34948"; NEUTRO <- "#f0efec"
TINTA <- "#0b0b0b"; TINTA2 <- "#52514e"; EJE <- "#8a8984"
tema <- theme_classic(base_size = 8) + theme(
  axis.line   = element_line(linewidth = 0.3, colour = EJE),
  axis.ticks  = element_line(linewidth = 0.3, colour = EJE),
  axis.text   = element_text(colour = TINTA2, size = 7),
  axis.title  = element_text(colour = TINTA, size = 8),
  plot.tag    = element_text(face = "bold", size = 10, colour = TINTA),
  legend.title = element_text(size = 7, colour = TINTA),
  legend.text  = element_text(size = 7, colour = TINTA2),
  legend.key.size = unit(3, "mm"),
  strip.background = element_blank(),
  strip.text  = element_text(face = "bold", size = 7.5, colour = TINTA),
  plot.margin = margin(6, 8, 6, 6))
coma <- function(x) format(x, big.mark = " ", decimal.mark = ",", scientific = FALSE, trim = TRUE)
# El PDF guarda la fecha y hora de creacion; se fija una fecha constante (con el
# mismo numero de caracteres) para que cada corrida produzca el mismo archivo.
fijar_fecha_pdf <- function(f) {
  x <- readBin(f, "raw", file.info(f)$size)
  for (clave in c("/CreationDate (D:", "/ModDate (D:")) {
    i <- grepRaw(clave, x, fixed = TRUE)
    if (length(i) == 1) x[i + nchar(clave) + 0:13] <- charToRaw("20000101000000")
  }
  writeBin(x, f)
}
guardar <- function(g, nombre, ancho, alto) {
  for (ext in c("pdf", "png"))
    ggsave(fig(paste0(nombre, ".", ext)), g, width = ancho, height = alto, units = "mm",
           dpi = 300, bg = "white")
  fijar_fecha_pdf(fig(paste0(nombre, ".pdf")))
  cat("  ", nombre, "\n")
}
juntar <- function(..., anchos = NULL) {       # paneles lado a lado, sin paquetes extra
  gs <- lapply(list(...), ggplotGrob); n <- length(gs)
  if (is.null(anchos)) anchos <- rep(1, n)
  hijos <- lapply(seq_len(n), function(i) gTree(children = gList(gs[[i]]), vp = viewport(layout.pos.col = i)))
  gTree(children = do.call(gList, hijos), vp = viewport(layout = grid.layout(1, n, widths = unit(anchos, "null"))))
}

# ---------------------------------------------------------------- datos -----
g   <- read.csv(res("genes", "tabla_de_genes.csv"), stringsAsFactors = FALSE)
ids <- read.delim("datos/string/identificadores.tsv", quote = "", stringsAsFactors = FALSE)
g$pref <- ids$preferredName[match(g$locus_antiguo, ids$queryItem)]
g$flagelar  <- g$conjunto_flagelar
g$chaperona <- grepl("Chaperona", g$categoria) & !g$flagelar

# Nombre corto de cada gen, en estilo de proteina (FliG, GroEL...)
nombre_corto <- function(d) {
  patron <- "\\b(Fl[aighb][A-Z]|Mot[A-Z]|Gro(EL|ES)|Dna[JK]|Clp[A-Z]|Htp[A-Z]|Sec[A-Z]|Hem[A-Z]|Hsp[0-9]+|Hfq|ATP12)\\b"
  n <- rep(NA_character_, nrow(d)); hay <- grepl(patron, d$producto)
  n[hay] <- regmatches(d$producto[hay], regexpr(patron, d$producto[hay]))   # alineado fila a fila
  n[grepl("^flagellin", d$producto)] <- "FlaA"
  n[grepl("^flagellar hook-associated family", d$producto)] <- "HAP3"
  mayus <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))
  otro <- ifelse(!is.na(d$gen), mayus(d$gen),
          ifelse(!is.na(d$pref) & !grepl("^ABM|^BARBA", d$pref), mayus(d$pref), sub("BARBAKC583_", "", d$locus)))
  n <- ifelse(is.na(n), otro, n)
  dup <- n %in% n[duplicated(n)]
  ifelse(dup, paste0(n, " (", sub("BARBAKC583_", "", d$locus), ")"), n)
}
g$nombre <- NA
sel <- g$flagelar | g$chaperona
g$nombre[sel] <- nombre_corto(g[sel, ])
g$nombre[!sel] <- sub("BARBAKC583_", "", g$locus[!sel])

# ======================================================= Figura 1: flujo ====
cat("Figuras:\n")
caja <- function(x, y, texto, borde = EJE, relleno = "white", ancho = 52, alto = 11)
  list(annotate("rect", xmin = x - ancho / 2, xmax = x + ancho / 2, ymin = y - alto / 2, ymax = y + alto / 2,
                fill = relleno, colour = borde, linewidth = 0.35),
       annotate("text", x = x, y = y, label = texto, size = 2.35, colour = TINTA, lineheight = 0.95))
flecha <- function(x1, y1, x2, y2)
  annotate("segment", x = x1, y = y1, xend = x2, yend = y2, colour = EJE, linewidth = 0.35,
           arrow = arrow(length = unit(1.3, "mm"), type = "closed"))
col <- c(30, 90, 150)
f1 <- ggplot() + coord_cartesian(xlim = c(0, 180), ylim = c(0, 102), expand = FALSE) + theme_void() +
  annotate("text", x = col, y = 98, size = 2.6, fontface = "bold", colour = TINTA,
           label = c("Genoma de referencia", "Datos de expresión", "Genómica comparada")) +
  caja(col[1], 86, "NCBI GCF_000015445.1\nB. bacilliformis KC583") +
  caja(col[1], 66, "Anotación: 1 216 genes\ncodificantes") +
  caja(col[1], 46, "Conjunto flagelar a priori\n31 genes (sin datos de expresión)", borde = FLAG) +
  caja(col[2], 86, "SRA PRJNA647605\n23 corridas pareadas") +
  caja(col[2], 66, "salmon 2.4.1\níndice con señuelo (genoma)") +
  caja(col[2], 46, "tximport + DESeq2\n12 muestras de placa y caldo") +
  caja(col[2], 26, "5 contrastes\nencogimiento ashr") +
  caja(col[3], 86, "Proteomas de 7 especies\ny de la cepa USM-LMMB07") +
  caja(col[3], 66, "BLASTp: panel del género\ny ortólogos recíprocos") +
  caja(col[3], 46, "STRING 12.0\nKEGG, GO y red") +
  flecha(col[1], 80.5, col[1], 71.5) + flecha(col[1], 60.5, col[1], 51.5) +
  flecha(col[2], 80.5, col[2], 71.5) + flecha(col[2], 60.5, col[2], 51.5) + flecha(col[2], 40.5, col[2], 31.5) +
  flecha(col[3], 80.5, col[3], 71.5) +
  caja(60, 7, "Enriquecimiento: hipergeométrico + fgsea\n7 categorías funcionales", borde = TINTA2, ancho = 64) +
  caja(135, 7, "Anotación independiente (STRING)\ny módulo flagelar en la red", borde = TINTA2, ancho = 64) +
  flecha(col[1], 40.5, 45, 12.5) + flecha(col[2], 20.5, 70, 12.5) + flecha(col[2], 20.5, 120, 12.5) +
  flecha(col[3], 40.5, 145, 12.5)
guardar(f1, "figura_1_flujo_de_trabajo", 180, 105)

# ======================================================= Figura 2: muestras =
mu  <- read.csv(res("muestras", "tabla_de_muestras.csv"), stringsAsFactors = FALSE)
dds <- readRDS(int("dds_23_muestras.rds"))
pv  <- round(100 * attr(plotPCA(vst(dds, blind = TRUE), intgroup = "cond", ntop = 500, returnData = TRUE), "percentVar"))
pc  <- read.csv(res("muestras", "pca_y_tasa_mapeo.csv"), stringsAsFactors = FALSE)
niv <- c("placa / caldo", "placa + atmosfera sanguinea", "sangre humana", "celulas endoteliales")
eti <- c("Placa / caldo", "Placa + atmósfera sanguínea", "Sangre humana", "Células endoteliales")
col_m <- setNames(c(AZUL, "#86b6ef", ROJO, "#4a3aa7"), eti)
for_m <- setNames(c(21, 24, 22, 23), eti)
pc$matriz <- factor(eti[match(pc$matriz, niv)], levels = eti)
f2a <- ggplot(pc, aes(PC1, PC2, fill = matriz, shape = matriz)) +
  geom_hline(yintercept = 0, colour = "grey90", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey90", linewidth = 0.3) +
  geom_point(size = 2.4, colour = "white", stroke = 0.5) +
  scale_fill_manual(values = col_m, name = NULL) + scale_shape_manual(values = for_m, name = NULL) +
  labs(x = paste0("PC1 (", pv[1], " % de la varianza)"), y = paste0("PC2 (", pv[2], " % de la varianza)"), tag = "A") +
  tema + theme(legend.position = "bottom", legend.direction = "vertical")

mu$matriz <- factor(eti[match(mu$matriz, niv)], levels = eti)
mu$etiqueta <- paste0(mu$condicion, "  ", mu$run)
mu$etiqueta <- factor(mu$etiqueta, levels = rev(mu$etiqueta))
mu$grupo <- factor(ifelse(mu$en_modelo, paste0("En el modelo (", sum(mu$en_modelo), ")"), paste0("Excluidas (", sum(!mu$en_modelo), ")")),
                   levels = paste0(c("En el modelo (", "Excluidas ("), c(sum(mu$en_modelo), sum(!mu$en_modelo)), ")"))
f2b <- ggplot(mu, aes(asignadas_CDS / 1e6, etiqueta, fill = matriz)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 2, linetype = "22", colour = TINTA2, linewidth = 0.35) +
  geom_text(data = data.frame(grupo = factor(levels(mu$grupo)[1], levels = levels(mu$grupo)), x = 2.1, y = Inf),
            aes(x, y, label = "2 millones"), inherit.aes = FALSE, hjust = 0, vjust = 1.2, size = 2.2, colour = TINTA2) +
  facet_grid(grupo ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = col_m, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Lecturas asignadas a genes bacterianos (millones)", y = NULL, tag = "B") +
  tema + theme(axis.text.y = element_text(size = 5.8, family = "mono"), axis.ticks.y = element_blank(),
               strip.text.y = element_text(angle = -90))
guardar(juntar(f2a, f2b, anchos = c(1, 1.15)), "figura_2_muestras", 180, 95)

# ======================================================= Figura 3: volcanes =
volcan <- function(lfc, padj, etiqueta_x, tag) {
  d <- data.frame(lfc = g[[lfc]], q = g[[padj]], flag = g$flagelar)[g$analizable, ]
  d <- d[!is.na(d$q), ]; d$y <- -log10(pmax(d$q, 1e-300))
  n_ind <- sum(d$q < 0.05 & d$lfc > 0); n_rep <- sum(d$q < 0.05 & d$lfc < 0)
  n_fl  <- sum(d$flag & d$q < 0.05)
  list(d = d, g = ggplot(d, aes(lfc, y)) +
    geom_hline(yintercept = -log10(0.05), linetype = "22", colour = EJE, linewidth = 0.3) +
    geom_vline(xintercept = c(-1, 1), linetype = "22", colour = EJE, linewidth = 0.3) +
    geom_point(data = d[!d$flag, ], colour = GRIS, size = 0.8, alpha = 0.7) +
    geom_point(data = d[d$flag, ], fill = FLAG, colour = "white", shape = 21, size = 1.9, stroke = 0.4) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 2.3, colour = TINTA2,
             label = paste0("Reprimidos: ", coma(n_rep))) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3, size = 2.3, colour = TINTA2,
             label = paste0("Inducidos: ", coma(n_ind))) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 3.0, size = 2.3, colour = TINTA,
             label = paste0("Flagelares (naranja) significativos: ", n_fl, " de ", sum(g$flagelar))) +
    labs(x = etiqueta_x, y = expression(-log[10]~"(q)"), tag = tag) + tema)
}
vt <- volcan("log2FC_37vs25", "padj_37vs25", expression(log[2]~"del cambio, 37 °C / 25 °C"), "A")
vp <- volcan("log2FC_pH6vspH8", "padj_pH6vspH8", expression(log[2]~"del cambio, pH 6 / pH 8"), "B")
lx <- range(c(vt$d$lfc, vp$d$lfc)) * 1.05; ly <- c(0, max(c(vt$d$y, vp$d$y)) * 1.12)
guardar(juntar(vt$g + coord_cartesian(xlim = lx, ylim = ly), vp$g + coord_cartesian(xlim = lx, ylim = ly)),
        "figura_3_volcanes_temperatura_y_ph", 180, 80)

# ======================================================= Figura 4: flagelo vs chaperonas
dds12 <- readRDS(int("dds_12_muestras.rds"))
cond  <- as.character(colData(dds12)$cond)
temp  <- c(Pl25 = "25 °C", Pl30 = "30 °C", Pl37 = "37 °C")
cols  <- which(cond %in% names(temp)); cols <- cols[order(cond[cols], colnames(dds12)[cols])]
vs    <- assay(vst(dds12, blind = FALSE))[, cols]
h     <- g[(g$flagelar | g$chaperona) & g$analizable, ]
h$grupo <- ifelse(h$flagelar, paste0("Flagelar (n = ", sum(h$flagelar), ")"), paste0("Chaperona (n = ", sum(h$chaperona), ")"))
z <- t(scale(t(vs[h$transcrito, , drop = FALSE])))
hz <- data.frame(nombre = rep(h$nombre, ncol(z)), grupo = rep(h$grupo, ncol(z)),
                 lfc = rep(h$log2FC_37vs25, ncol(z)),
                 muestra = rep(colnames(z), each = nrow(z)), z = as.vector(z))
hz$temperatura <- factor(temp[cond[match(hz$muestra, colnames(dds12))]], levels = temp)
hz$replica <- ave(hz$muestra, hz$temperatura, FUN = function(m) paste0("R", match(m, sort(unique(m)))))
hz$nombre <- factor(hz$nombre, levels = unique(h$nombre[order(h$flagelar, -h$log2FC_37vs25)]))
hz$grupo  <- factor(hz$grupo, levels = unique(h$grupo[order(!h$flagelar)]))
f4a <- ggplot(hz, aes(replica, nombre, fill = pmax(pmin(z, 2), -2))) +
  geom_tile(colour = "white", linewidth = 0.4) +
  facet_grid(grupo ~ temperatura, scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = AZUL, mid = NEUTRO, high = ROJO, midpoint = 0, limits = c(-2, 2),
                       name = "Puntuación Z", breaks = c(-2, 0, 2),
                       guide = guide_colourbar(title.position = "top", barwidth = unit(25, "mm"), barheight = unit(2.5, "mm"))) +
  labs(x = NULL, y = NULL, tag = "A") + tema +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = 5.8),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 90),
        legend.position = "bottom", panel.spacing = unit(1, "mm"))

cv <- rbind(data.frame(nombre = h$nombre, grupo = h$grupo, t = 25, lfc = 0),
            data.frame(nombre = h$nombre, grupo = h$grupo, t = 30, lfc = h$log2FC_30vs25),
            data.frame(nombre = h$nombre, grupo = h$grupo, t = 37, lfc = h$log2FC_37vs25))
cv$grupo <- factor(cv$grupo, levels = levels(hz$grupo))
med <- aggregate(lfc ~ grupo + t, cv, median)
col_g <- setNames(c(FLAG, CHAP), levels(hz$grupo))
f4b <- ggplot(cv, aes(t, lfc, colour = grupo)) +
  geom_hline(yintercept = 0, colour = EJE, linewidth = 0.3) +
  geom_line(aes(group = nombre), linewidth = 0.3, alpha = 0.35) +
  geom_line(data = med, linewidth = 1.2) +
  geom_point(data = med, size = 2.2, shape = 21, fill = "white", stroke = 1) +
  scale_colour_manual(values = col_g, name = NULL) +
  scale_x_continuous(breaks = c(25, 30, 37), labels = temp) +
  labs(x = "Temperatura de cultivo", y = expression(log[2]~"del cambio respecto a 25 °C"), tag = "B") +
  tema + theme(legend.position = "bottom", legend.direction = "vertical")
guardar(juntar(f4a, f4b, anchos = c(1.25, 1)), "figura_4_flagelares_y_chaperonas", 180, 165)
cat("    medianas de la curva:\n"); print(med, row.names = FALSE)

# ======================================================= Figura 5: red STRING
nod <- read.csv(res("string", "red_cambio_2x_nodos.csv"), stringsAsFactors = FALSE)
ari <- read.csv(res("string", "red_cambio_2x_aristas.csv"), stringsAsFactors = FALSE)
nod <- nod[nod$grado > 0, ]
red <- graph_from_data_frame(ari[, c("a", "b", "score")], directed = FALSE, vertices = nod$string_id)
set.seed(42)
xy  <- layout_with_fr(red, weights = E(red)$score, niter = 2000)
nod <- cbind(nod, x = xy[, 1], y = xy[, 2])
m <- match(nod$locus, g$locus)
nod$nombre <- g$nombre[m]
nod$clase <- ifelse(g$flagelar[m], "Flagelar", ifelse(g$chaperona[m], "Chaperona",
             ifelse(nod$log2FC_37vs25 < 0, "Otro reprimido", "Otro inducido")))
nod$clase <- factor(nod$clase, levels = c("Flagelar", "Chaperona", "Otro reprimido", "Otro inducido"))
seg <- data.frame(x = nod$x[match(ari$a, nod$string_id)], y = nod$y[match(ari$a, nod$string_id)],
                  xend = nod$x[match(ari$b, nod$string_id)], yend = nod$y[match(ari$b, nod$string_id)],
                  score = ari$score)
f5 <- ggplot() +
  geom_segment(data = seg, aes(x, y, xend = xend, yend = yend, alpha = score), colour = "grey55", linewidth = 0.25) +
  geom_point(data = nod, aes(x, y, fill = clase, shape = clase), size = 2.1, colour = "white", stroke = 0.4) +
  geom_text_repel(data = nod[nod$clase %in% c("Flagelar", "Chaperona"), ], aes(x, y, label = nombre),
                  size = 2, colour = TINTA, segment.size = 0.2, segment.colour = EJE,
                  max.overlaps = Inf, box.padding = 0.15, min.segment.length = 0.2, seed = 42,
                  max.time = 60, max.iter = 10000) +   # termina por iteraciones, no por tiempo
  # Aun asi, ggrepel 0.9.8 puede colocar algunas etiquetas en otra posicion entre
  # corridas (comprobado: misma disposicion de la red, dos versiones del dibujo).
  scale_fill_manual(values = c(Flagelar = FLAG, Chaperona = CHAP, "Otro reprimido" = GRIS_O, "Otro inducido" = GRIS), name = NULL) +
  scale_shape_manual(values = c(Flagelar = 25, Chaperona = 24, "Otro reprimido" = 25, "Otro inducido" = 24), name = NULL) +
  scale_alpha_continuous(range = c(0.15, 0.7), name = "Puntaje STRING") +
  coord_equal() + theme_void(base_size = 8) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7, colour = TINTA2),
        legend.title = element_text(size = 7), plot.margin = margin(6, 6, 6, 6))
guardar(f5, "figura_5_red_string", 180, 160)
cat("    red: ", nrow(nod), "genes con al menos una arista;", nrow(ari), "aristas\n")

# ======================================================= Figura 6: genero ===
# Solo cuentan los ortologos reciprocos (script 06, busqueda de vuelta en KC583).
# Un acierto cuya mejor coincidencia de vuelta es otra proteina de KC583 es una
# paraloga: se marca con una cruz y no se cuenta.
pn  <- read.delim(res("genomica_comparada", "panel_genero_blastp.tsv"), header = FALSE, col.names = c("wp", "especie", "ident"))
vu  <- read.delim(res("genomica_comparada", "panel_genero_busqueda_de_vuelta.tsv"), stringsAsFactors = FALSE)
fl  <- g[g$flagelar, ]; fl <- fl[order(fl$locus), ]
esp <- c("bacilliformis", "ancashensis", "clarridgeiae", "schoenbuchensis", "tribocorum", "quintana", "henselae")
pan <- expand.grid(wp = fl$proteina, especie = esp, stringsAsFactors = FALSE)
k   <- match(paste(pan$wp, pan$especie), paste(vu$proteina_kc583, vu$especie))
pan$reciproco <- ifelse(pan$especie == "bacilliformis", !is.na(match(pan$wp, pn$wp[pn$especie == "bacilliformis"])),
                        !is.na(k) & vu$reciproco[k] %in% "si")
pan$paraloga  <- !is.na(k) & vu$reciproco[k] %in% "no"
pan$ident <- ifelse(pan$especie == "bacilliformis", pn$ident[match(paste(pan$wp, pan$especie), paste(pn$wp, pn$especie))],
                    vu$identidad[k])
pan$ident[!pan$reciproco] <- NA
pan$nombre <- factor(fl$nombre[match(pan$wp, fl$proteina)], levels = fl$nombre)
cuenta <- tapply(pan$reciproco, pan$especie, sum)
pan$especie <- factor(paste0("B. ", pan$especie, "  (", cuenta[pan$especie], "/", nrow(fl), ")"),
                      levels = rev(paste0("B. ", esp, "  (", cuenta[esp], "/", nrow(fl), ")")))
f6 <- ggplot(pan, aes(nombre, especie)) +
  geom_tile(aes(fill = ident), colour = "white", linewidth = 0.5) +
  geom_point(data = pan[pan$paraloga, ], aes(shape = "Acierto parálogo (no recíproco)"), size = 1.6, stroke = 0.5, colour = TINTA2) +
  scale_fill_gradient(low = "#cde2fb", high = "#0d366b", limits = c(20, 100), na.value = "#f4f3f1",
                      name = "Identidad (%)\ndel ortólogo recíproco") +
  scale_shape_manual(values = 4, name = NULL) +
  labs(x = NULL, y = NULL) + tema +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = element_text(face = "italic", size = 7), legend.position = "right",
        legend.title = element_text(size = 7), legend.text = element_text(size = 6.5))
guardar(f6, "figura_6_panel_genero", 180, 62)
