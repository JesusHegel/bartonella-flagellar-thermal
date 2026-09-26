# 09. Tabla de muestras: lecturas de cada corrida y por que entra o no al modelo.
#
# Datos de los informes de salmon (aux_info/meta_info.json):
#   num_processed         pares de lecturas procesados
#   num_mapped            asignados a CDS bacterianas (objetivos del indice)
#   num_decoy_fragments   asignados al senuelo (resto del cromosoma bacteriano)
# El indice es un gentrome, asi que lo que no se asigna a CDS ni al senuelo no
# se parece al genoma de KC583.
source("scripts/rutas.R")

leer <- function(f, campo) {
  l <- grep(paste0('"', campo, '"'), readLines(f, warn = FALSE), value = TRUE)[1]
  as.numeric(sub('.*:[[:space:]]*([0-9.eE+-]+).*', '\\1', l))
}
tabla <- read.delim("datos/muestras.tsv", stringsAsFactors = FALSE)
js <- file.path(DIR_CUANT, tabla$run, "aux_info", "meta_info.json")
stopifnot(all(file.exists(js)))
tabla$lecturas_procesadas <- vapply(js, leer, numeric(1), campo = "num_processed")
tabla$asignadas_CDS       <- vapply(js, leer, numeric(1), campo = "num_mapped")
tabla$asignadas_senuelo   <- vapply(js, leer, numeric(1), campo = "num_decoy_fragments")
tabla$pct_CDS     <- round(100 * tabla$asignadas_CDS / tabla$lecturas_procesadas, 2)
tabla$pct_senuelo <- round(100 * tabla$asignadas_senuelo / tabla$lecturas_procesadas, 2)
tabla$pct_sin_asignar <- round(100 - tabla$pct_CDS - tabla$pct_senuelo, 2)

tabla$matriz <- ifelse(tabla$cond %in% c("HB37", "HBBG"), "sangre humana",
                ifelse(tabla$cond == "HUVE", "celulas endoteliales",
                ifelse(tabla$cond == "PlBG", "placa + atmosfera sanguinea", "placa / caldo")))
serie <- c(Pl25 = "temperatura", Pl30 = "temperatura", Pl37 = "temperatura",
           pH06 = "pH", pH07 = "pH", pH08 = "pH")
tabla$en_modelo <- tabla$cond %in% names(serie)
tabla$motivo <- ifelse(tabla$en_modelo, paste0("serie de ", serie[tabla$cond]),
                ifelse(tabla$cond == "PlBG", "no pertenece a las series de temperatura ni de pH",
                       "matriz distinta (separada en PC1) y pocas lecturas asignadas a CDS"))

orden <- c("Pl25", "Pl30", "Pl37", "pH06", "pH07", "pH08", "PlBG", "HB37", "HBBG", "HUVE")
tabla <- tabla[order(match(tabla$cond, orden), tabla$run), ]
names(tabla)[names(tabla) == "cond"] <- "condicion"
write.csv(tabla, res("muestras", "tabla_de_muestras.csv"), row.names = FALSE)

print(tabla[, c("run", "condicion", "lecturas_procesadas", "asignadas_CDS", "pct_CDS", "pct_senuelo", "en_modelo")],
      row.names = FALSE)
cat("\nLecturas asignadas a CDS por grupo (millones):\n")
tabla$grupo <- ifelse(tabla$en_modelo, "modelo (12)", tabla$matriz)
print(aggregate(asignadas_CDS ~ grupo, tabla, function(x) round(range(x) / 1e6, 3)))
