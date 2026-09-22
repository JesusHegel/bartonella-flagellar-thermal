# Tablas de resultados

Regeneradas el 22 de septiembre de 2026 por los scripts de `../scripts/`.
Las cifras autorizadas están en [`../CIFRAS_CONFIRMADAS.md`](../CIFRAS_CONFIRMADAS.md).

| Archivo | Contenido | Script |
|---|---|---|
| `chk_universo.csv` | Universo deduplicado, inducidos y reprimidos | `s03` |
| `chk_conteos.csv` | Transcritos modificados bajo cada definición de umbral | `s04` |
| `chk_categorias.csv` | Enriquecimiento hipergeométrico por categoría funcional | `s03` |
| `chk_fgsea.csv` | Enriquecimiento sobre el ranking completo, sin umbrales | `s03` |
| `chk_estricto.csv` | Sensibilidad frente al umbral de magnitud | `s04` |
| `chk_yale.csv` | Modelo global frente a modelo restringido, gen por gen | `s04` |
| `chk_particiones.csv` | Concordancia entre dos particiones independientes | `s04` |
| `chk_curva31.csv` | Los 31 genes flagelares a 25, 30 y 37 °C | `s04` |
| `chk_especificidad.csv` | Respuesta a temperatura frente a respuesta a pH | `s04` |
| `chk_pca_mapeo.csv` | PC1, PC2, matriz biológica y tasa de mapeo por muestra | `s05` |
| `chk_ranking.csv` | Quince transcritos más abundantes a 37 °C | `s04` |
| `chk_top_global.csv` | Quince transcritos más abundantes por media global | `s04` |
| `chk_quimiotaxis.tsv` | Homólogos del sistema Che en KC583 | `s08` |
| `panel.tsv` | Identidad por BLASTp de los 31 flagelares en siete especies | `s07` |
| `rbh.tsv` | Ortólogos recíprocos KC583 ↔ USM-LMMB07 | `s07` |

## `sin_regenerar/`

Cuatro tablas de la corrida de julio de 2026 que **no se han vuelto a generar**
con el pipeline actual. Se conservan porque contienen resultados que no existen
en otra parte, pero sus cifras **no están validadas** y no deben citarse sin
recalcularlas.

| Archivo | Contenido | Nota |
|---|---|---|
| `chk_dichter.csv` | Mann-Whitney sobre reactividad serológica; resultado negativo (p = 1) | por regenerar |
| `chk_abundancia.csv` | Conteos de transcritos seleccionados por condición | por regenerar; además compara genes distintos entre sí, lo que los conteos de DESeq2 no permiten |
| `chk_replicas.csv` | Conteos por réplica individual | por regenerar |
| `chk_share.csv` | Porcentaje de librería cruda por transcrito | por regenerar |
