# Cifras confirmadas

Regeneración completa del análisis el 22 de septiembre de 2026, en carpeta
limpia, con comparación contra la versión anterior.
**Este archivo es la única fuente válida. Cualquier cifra que no esté aquí
no debe usarse.**

---

## 1. Origen de los datos

| Recurso | Identificador |
|---|---|
| RNA-seq | BioProject PRJNA647605 (Wachter et al. 2020, PLoS Negl Trop Dis 14(11):e0008671) |
| Genoma de referencia | GCF_000015445.1 — *B. bacilliformis* KC583, cromosoma NC_008783.1 |
| Cepa peruana (ortología) | GCF_001624625.1 — USM-LMMB07 |
| Panel de género | GCF_000253015.1 · GCF_039555305.1 · GCF_001281405.1 · GCF_000196435.1 · GCF_009936175.1 · GCF_019930925.1 |
| Quimiotaxis | GCF_037023865.1 · GCF_000022005.1 · GCF_000092025.1 · GCF_003324715.1 |

Las 23 corridas provienen de **un solo centro** (The University of Montana),
mismo instrumento (Illumina HiSeq 2500) y mismo diseño pareado, según los
metadatos de SRA. **No existen dos centros de secuenciación.**

---

## 2. Genoma y anotación

| Magnitud | Valor |
|---|---|
| CDS anotadas en el GFF (filas) | **1 217** |
| Genes codificantes (locus distintos) | **1 216** (prfB, RS03925, ocupa dos filas CDS por desplazamiento ribosómico programado) |
| Pseudogenes (sin `protein_id`) | 50 |
| CDS con `protein_id` | 1 166 |
| Accesiones de proteína únicas | **1 139** (coincide con `protein.faa` de NCBI) |
| Accesiones duplicadas (región de ~28 kb) | 27 |
| Genes con anotación flagelar | **31** |

El conjunto flagelar se define en `scripts/01_anotacion.sh`
por patrón sobre el campo `product` del GFF. **Ese script no consulta datos de
expresión.** Los 31 son la totalidad de los anotados, no una selección.

---

## 3. Muestras y calidad de cuantificación

23 corridas, 10 condiciones. Réplicas: 3 en HB37, HBBG y PlBG; 2 en las demás.

Tasa de mapeo por matriz biológica (media):

| Matriz | Condiciones | Mapeo | PC1 medio |
|---|---|---|---|
| Placa / caldo | Pl25, Pl30, Pl37, pH06, pH07, pH08 | 38,3 % | +14,7 |
| Placa + atmósfera sanguínea | PlBG | 45,5 % | −3,4 |
| Sangre humana | HB37, HBBG | 4,6 % | −23,2 |
| Células endoteliales | HUVE | **0,3 %** | −13,8 |

Entre el 0,2 y el 6 % de las lecturas de sangre y de células endoteliales se
asignaron a secuencias codificantes bacterianas, frente al 38 % en cultivo puro.
Esa diferencia de profundidad, junto con la PCA (sección 4), justifica excluir
esas muestras del modelo de temperatura y pH (lecturas por grupo más abajo).

El índice de salmon (`idx_kc583`) es un gentrome: las CDS como objetivo y el
cromosoma como señuelo. Por tanto `percent_mapped` mide **asignación a CDS**, no
contenido bacteriano total. El señuelo (cromosoma fuera de las CDS) recibe una
fracción grande **en todas las matrices bacterianas**: 39–70 % de los fragmentos
en placa y caldo, 56–66 % en sangre, y solo 0,8–1,0 % en células endoteliales
(`resultados/muestras/tabla_de_muestras.csv`). Por su tamaño es compatible con
ARN ribosomal no eliminado (el estudio original buscaba ARN pequeños), pero **no
se ha verificado**. No debe interpretarse como fracción humana.

| Transcritos | Valor |
|---|---|
| En el índice de salmon | 1 189 |
| Analizables (`rowSums > 0`) | **1 186** |

Lecturas asignadas a CDS por muestra (`num_mapped` de salmon):

| Grupo | Muestras | Lecturas asignadas a CDS |
|---|---|---|
| Placa / caldo (modelo restringido) | 12 | 1,04 – 7,89 millones |
| Placa + atmósfera sanguínea (PlBG) | 3 | 3,25 – 3,65 millones |
| Sangre humana (HB37, HBBG) | 6 | 0,21 – 0,88 millones |
| Células endoteliales (HUVE) | 2 | 0,045 – 0,059 millones |

Referencia: Haas et al. 2012 (BMC Genomics 13:734) indican que 2–3 millones de
fragmentos no ribosomales por muestra permiten detectar un número importante de
genes con cambios de 2 veces o más, y que 5–10 millones detectan casi todos los
genes salvo los de expresión más baja. Dos muestras incluidas quedan por debajo
de 2 millones (Pl25 SRR12653236: 1,63; Pl30 SRR12653239: 1,04). **La profundidad
sola no separa los grupos**: la muestra de sangre más profunda (0,88) está cerca
de la de cultivo más baja (1,04). La exclusión de sangre y células endoteliales
se apoya en la PCA (sección 4) y en el diseño; la de PlBG, solo en el diseño (no
pertenece a las series de temperatura ni de pH).

---

## 4. Estructura de los datos (PCA)

`vst(blind=TRUE)`, `plotPCA(ntop=500)`.

- **PC1 = 53 %** de la varianza · **PC2 = 16 %**
- Correlación PC1 ~ tasa de mapeo: **0,817**
- Correlación PC2 ~ tasa de mapeo: 0,06

PC1 separa por **matriz biológica**, que determina simultáneamente qué
fracción de las lecturas se asigna a genes bacterianos. Ambos factores son
inseparables en este diseño.

**Modelo restringido**: 12 muestras en matriz homogénea
(Pl25, Pl30, Pl37, pH06, pH07, pH08).

---

## 5. Expresión diferencial

`lfcShrink(type="ashr")`. Universo deduplicado: **1 138** identificadores de
proteína únicos.

| Contraste | padj<0,05 | padj<0,05 y \|LFC\|>1 | Flagelares sig. |
|---|---|---|---|
| Global 37 vs 25 | 538 | 214 | 30 / 31 |
| **Restringido 37 vs 25** | **488** | **207** | **31 / 31** |
| Restringido 30 vs 25 | 423 | 146 | 30 / 31 |
| Restringido pH06 vs pH08 | 375 | 74 | **6 / 31** |

Umbral formal (`results(lfcThreshold=1)`): 49 significativos, 7 flagelares.

De los 31 flagelares: **30 reprimidos y 1 inducido** (RS04185, FliO, +0,87).

---

## 6. Enriquecimiento

Modelo restringido. Universo 1 138 · 200 inducidos · 278 reprimidos.

### Hipergeométrico (corrección Benjamini-Hochberg)

| Categoría | Genes | Ind. | Repr. | Esperados repr. | q |
|---|---|---|---|---|---|
| **Flagelar** | 31 | 1 | **30** | 7,6 | **2,14 × 10⁻¹⁷** (repr.) |
| **Chaperona** | 15 | **8** | 2 | 3,7 | **1,17 × 10⁻²** (ind.) |
| **TransporteABC** | 53 | 3 | 22 | 12,9 | **1,31 × 10⁻²** (repr.) |
| Ribosoma | 56 | 7 | 3 | 13,7 | no significativo |
| HierroHemo | 24 | 4 | 8 | 5,9 | no significativo |
| Traducción | 59 | 8 | 12 | 14,4 | no significativo |
| Hipotética | 86 | 18 | 19 | 21,0 | no significativo |

### fgsea (semilla 42, 100 000 permutaciones, sin umbrales)

| Vía | Tamaño | NES | q |
|---|---|---|---|
| **Flagelar** | 31 | **−2,44** | **2,43 × 10⁻⁹** |
| TransporteABC | 53 | −1,87 | 5,00 × 10⁻⁴ |
| Chaperona | 15 | **+1,79** | 1,15 × 10⁻² |
| HierroHemo | 24 | −1,42 | 0,111 |
| Ribosoma | 56 | +1,19 | 0,187 |
| Hipotética | 86 | +1,19 | 0,187 |
| Traducción | 59 | −0,51 | 1,000 |

---

## 7. Robustez

| Prueba | Resultado |
|---|---|
| Diferencia máxima de LFC entre modelo global y restringido (31 genes) | **0,09** |
| Genes flagelares con la misma dirección en dos particiones independientes | **31 de 31** |
| Correlación entre particiones (genoma completo) | 0,581 |
| Especificidad: flagelares significativos a temperatura | **31 / 31** |
| Especificidad: flagelares significativos a pH | **6 / 31** |
| Fondo del genoma: temperatura / pH | 488 / 375 |

El resultado no depende del modelo, de la partición de réplicas ni del umbral,
y no aparece frente a un estímulo distinto de magnitud comparable.

---

## 8. Perfil térmico: dos programas con perfiles distintos

De los 31 genes flagelares, el contraste **37 vs 30 °C no es significativo en
30 de 31** (único con p<0,05: FliF RS05510, p = 0,012, que no sobrevive la
corrección por 31 pruebas). **23 de 31 alcanzan su mínimo a 30 °C.**

Conteos normalizados (media por condición; `resultados/robustez/medias_normalizadas_25_30_37C.csv`).
Sirven para comparar un mismo gen entre temperaturas, no genes distintos entre sí:

| Gen | Producto | 25 °C | 30 °C | 37 °C |
|---|---|---|---|---|
| RS05045 | **Flagelina** | **57 904** | **16 796** | **12 553** |
| RS05735 | **Chaperonina GroEL** | **41 002** | **43 266** | **199 585** |
| RS01390 | Factor sigma RpoH | 18 035 | 41 771 | 60 750 |
| RS06465 | Chaperona DnaK | 13 634 | 22 565 | 53 793 |
| RS02960 | Hsp20 | 6 681 | 20 740 | 43 968 |
| RS05740 | Co-chaperona GroES | 4 378 | 6 629 | 24 357 |
| RS01060 | Chaperona ClpB | 1 315 | 1 827 | 15 181 |
| RS02640 | Endopeptidasa La (Lon) | 3 992 | 6 089 | 12 609 |
| RS05920 | Proteína de membrana externa | 58 505 | 11 347 | 15 462 |

**La represión flagelar ocurre entre 25 y 30 °C** (flagelina −71 %) y no avanza
entre 30 y 37 °C. **Las chaperonas principales cambian sobre todo entre 30 y
37 °C**: GroEL +5 % de 25 a 30 °C y ×4,6 de 30 a 37 °C; ClpB ×1,4 y ×8,3.
No todo el sistema de choque térmico espera a 37 °C: RpoH (×2,3), Hsp20 (×3,1)
y DnaK (×1,7) ya suben entre 25 y 30 °C (razones entre las medias de la tabla).
Los dos programas tienen perfiles térmicos distintos; los datos no permiten
afirmar que uno preceda al otro.

Lectura correcta: la expresión flagelar es **máxima a 25 °C** y ya está
reprimida a 30 °C, sin descenso adicional a 37 °C.

---

## 9. Genómica comparada

### Panel de género (BLASTp, e<1e-5, de los 31 flagelares)

Cada acierto se buscó de vuelta en KC583 (`scripts/06_panel_genero_y_ortologia.sh`,
`resultados/genomica_comparada/panel_genero_busqueda_de_vuelta.tsv`): es ortólogo
recíproco si su mejor coincidencia de vuelta (mayor bitscore) es el mismo gen
flagelar.

| Especie | Aciertos BLASTp | **Ortólogos recíprocos** |
|---|---|---|
| *B. bacilliformis* | 31 / 31 | 31 / 31 |
| *B. clarridgeiae* | 31 / 31 | **31 / 31** |
| *B. schoenbuchensis* | 31 / 31 | **31 / 31** |
| *B. ancashensis* | 31 / 31 | **31 / 31** |
| *B. tribocorum* | 2 / 31 | **1 / 31** (FliO) |
| *B. quintana* | 2 / 31 | **1 / 31** (FliO) |
| *B. henselae* | 2 / 31 | **1 / 31** (FliO) |

- **FliO** (RS04185): 45,5 % (*B. tribocorum*), 46,1 % (*B. quintana*) y 49,4 %
  (*B. henselae*) de identidad, 100 % de cobertura, anotada como FliO en las tres.
- **El acierto de FliI (RS05600) es un parálogo:** la subunidad β de la ATP
  sintasa F0F1 (identidad 29,9–30,9 %, cobertura 54–57 %), cuya mejor
  coincidencia de vuelta en KC583 es la ATP sintasa β (WP_005765873.1), no FliI.

"2 de 31" queda retirado (sección 11): la cifra correcta es **1 de 31**.

### Ortología con la cepa peruana USM-LMMB07

Mejores coincidencias recíprocas (BLASTp, e<1e-10): **1 111 pares**
(1 137 proteínas en USM, 1 139 en KC583).
**Los 31 genes flagelares tienen ortólogo recíproco.**

### Quimiotaxis

94 proteínas Che de cuatro alfaproteobacterias con sistema completo, contra el
proteoma de KC583: **10 aciertos** (e<1e-3) sobre **7 proteínas**, con
identidades de 22–35 % y alineamientos parciales (el mejor cubre 194 de 1 160
aminoácidos). **KC583 no posee un sistema Che completo.**

---

## 10. Estructura del 5'UTR de la flagelina

Región intergénica de 141 nt aguas arriba de RS05045 (hebra negativa,
complemento reverso), plegada con RNAfold.

| | Valor |
|---|---|
| MFE del 5'UTR | **−39,5 kcal/mol** |
| 100 controles aleatorios del mismo genoma, 141 nt | media −26,6 · mediana −26,6 · percentil 5 −38,1 |
| Controles más estructurados que el UTR | **4 %** |

El 5'UTR está entre el 4 % más estructurado del genoma, lo que es compatible
con control post-transcripcional. Compatible, no demostrado.

---

## 11. Cifras retiradas

| Cifra anterior | Valor correcto | Motivo |
|---|---|---|
| Yale n=12 vs GENEWIZ n=11 | **un solo centro** | Metadatos de SRA: todas las corridas son de University of Montana |
| 1 216 CDS | **1 217 filas CDS = 1 216 genes** | prfB (RS03925) tiene dos segmentos; ambas cifras son correctas según se cuenten filas o genes |
| 1 176 analizables | **1 186** | Recuento anterior no reproducible |
| Universo 1 134 (o 1 175) | **1 138** | `lt2prot.tsv` anterior incluía 50 CDS sin `protein_id` |
| 29 de 31 flagelares | **31 de 31 significativos** (30 repr. + 1 ind.) | |
| 216 transcritos modificados | **207** (restringido, \|LFC\|>1) o **488** (padj<0,05) | Definición anterior no rastreable |
| q = 1,3 × 10⁻¹⁶ | **2,14 × 10⁻¹⁷** | |
| NES −2,46 · q 1,8 × 10⁻⁹ | **−2,44 · 2,43 × 10⁻⁹** | |
| 6 flagelares con umbral formal | **7** | |
| 2 de 31 flagelares en *B. tribocorum*, *B. quintana* y *B. henselae* | **1 de 31** (FliO) | El segundo acierto, de FliI, es la ATP sintasa β (parálogo); búsqueda de vuelta en KC583, 26-09-2026 |

---

## 12. Reproducibilidad

Todas las cifras de este documento se regeneran con los scripts de `scripts/`,
ejecutados en ese orden por `scripts/ejecutar_todo.sh`:

| Script | Produce | Tablas |
|---|---|---|
| `01_anotacion.sh` | sección 2 | `resultados/anotacion/` |
| `02_expresion_diferencial.R` | secciones 3, 4 y 5 | `resultados/expresion_diferencial/` |
| `03_enriquecimiento.R` | sección 6 | `resultados/enriquecimiento/` |
| `04_robustez.R` | secciones 5, 7 y 8 | `resultados/robustez/` |
| `05_pca_y_tasa_mapeo.R` | secciones 3 y 4 | `resultados/muestras/` |
| `06_panel_genero_y_ortologia.sh` | sección 9 (panel de género y ortología) | `resultados/genomica_comparada/` |
| `07_quimiotaxis_y_utr.sh` | secciones 9 (quimiotaxis) y 10 | `resultados/genomica_comparada/` |
| `08_tabla_de_genes.R` | tabla de genes y top 20 | `resultados/genes/` |
| `09_tabla_de_muestras.R` | sección 3 | `resultados/muestras/` |
| `10_descargar_string.sh` | copia de STRING | `datos/string/` |
| `11_analisis_string.R` | análisis con STRING | `resultados/string/` |
| `12_figuras.R` | sección 15 | `figuras/` |

`scripts/verificar.sh` compara cada tabla regenerada con la versión publicada en
git. La réplica en otra computadora se describe en
`GUIA_DE_REPLICA.md`.

Nombres anteriores al 25-09-2026 (para leer las secciones 12b y los registros
antiguos): `results/chk_*.csv` → `resultados/`; `chk_yale.csv` →
`modelo_global_vs_restringido.csv`; `chk_curva31.csv` →
`flagelares_curva_termica.csv`; `res_rY.csv` → `37C_vs_25C.csv`; `res_rG.csv` →
`37C_vs_25C_modelo_23_muestras.csv`; `panel.tsv` → `panel_genero_blastp.tsv`;
`rbh.tsv` → `ortologos_reciprocos_KC583_USM-LMMB07.tsv`; `quants/` →
`cuantificacion/`; los scripts 06, 07 y 08 anteriores son ahora 08, 06 y 07.

## 12b. Verificación de reproducibilidad (25 de septiembre de 2026)

La cadena publicada se ejecutó de principio a fin en una carpeta `salida/` vacía,
partiendo de `quants/` y de la referencia. Los scripts 01 a 05 regeneraron los
doce archivos de `results/` de forma **idéntica byte a byte** a los publicados.

Corrección aplicada en esa verificación: el script 05 requería `mapeo.tsv`, que
ningún script generaba. Ahora lo construye a partir de los informes
`quants/*/aux_info/meta_info.json` de salmon.

Secciones 9 y 10 reverificadas el mismo día en `salida/` limpia: los scripts 07
y 08 regeneraron `panel.tsv`, `rbh.tsv`, `chk_quimiotaxis.tsv` y
`utr_flagelina.fa` **idénticos byte a byte**, y los valores del UTR (−39,5
kcal/mol; controles media −26,6; 4 %) coinciden.

### Cadena desde el SRA

- Índice: `salmon index -t gentrome.fna -d decoys.txt -i idx_kc583 -k 31`.
  `gentrome.fna` = `cds_from_genomic.fna` + genoma completo concatenados (md5
  `bc9bf0e33d39363355d2e013eaf31e44`); `decoys.txt` = `NC_008783.1`. Reconstruido
  el 25-09-2026: las seis huellas de `info.json` (secuencias, nombres, señuelo)
  coinciden con el índice original.
- Cuantificación, igual en las 23 muestras: `salmon quant -i idx_kc583 -l A
  -1 R_1.fastq -2 R_2.fastq -p 6`, con lecturas obtenidas por `prefetch` y
  `fasterq-dump --split-3`.
- Prueba con SRR12653239 descargada de nuevo: mismas 2 786 511 lecturas;
  asignadas 1 040 550 frente a 1 040 336 (+0,02 %); 8 transcritos con NumReads
  distinto, ninguno flagelar; el mayor cambio en RS00010 (615 → 851). salmon con
  varios hilos no es determinista: desde el SRA los resultados se reproducen con
  diferencias mínimas, no byte a byte.

## 12c. Entorno de referencia (25 de septiembre de 2026)

Desde esta fecha, las tablas publicadas se generan con el entorno de conda
`entorno/bartonella.yml`, que reúne en uno solo los programas antes repartidos
en cuatro entornos. Comparado con la corrida anterior (entorno `de`), en la
misma computadora:

- 22 de 29 tablas **idénticas byte a byte**.
- Los 5 contrastes de DESeq2, fgsea y la PCA difieren solo en decimales lejanos:
  como máximo 1,4 × 10⁻⁵ en log2FoldChange, 3,9 × 10⁻⁵ en pvalue y 1 × 10⁻⁶ en
  el NES. **Ningún gen cambia de significancia ni de signo**; los valores p de
  fgsea son idénticos y ninguna cifra de este documento cambia.
- Dos corridas seguidas con el entorno `bartonella` dan tablas idénticas.

La causa es la librería de álgebra lineal de cada instalación. Por eso
`verificar.sh` compara esas 7 tablas con tolerancia (`comparar_con_tolerancia.R`:
diferencia menor a 1e-4 absoluta o 1 % relativa, sin cambios de significancia ni
de signo) y todas las demás byte a byte.

**Orden alfabético (26 de septiembre de 2026).** `rutas.sh` y `rutas.R` fijan el
orden alfabético "C" (`LC_COLLATE=C`, mayúsculas antes que minúsculas) y dejan
el texto en UTF-8 para los acentos de las figuras. Con el orden del idioma del
sistema, R ordena distinto las condiciones de DESeq2 (cambian los últimos
decimales) y las listas de genes de STRING (flgA antes que RS05045). Con el
orden "C", la corrida completa del 26-09-2026 dio **38 de 38 tablas idénticas
byte a byte** a las publicadas, y lo mismo desde un clon nuevo del repositorio.

## 13. Versiones de programas

```
salmon 2.4.1
R version 4.5.3 (2026-03-11) 
DESeq2 1.50.2 
tximport 1.38.2 
ashr 2.2.63 
fgsea 1.36.2 
blastp: 2.17.0+
RNAfold 2.7.0
datasets version: 18.33.1
sra-tools 3.4.1
Python 3.12 (en el entorno bartonella; 3.13 no es compatible con viennarna 2.7.0)
```

## 14. Tabla de genes y análisis con STRING (Fase 2, 25–26 de septiembre de 2026)

### Tabla de genes (`08_tabla_de_genes.R`, `resultados/genes/tabla_de_genes.csv`)

| Magnitud | Valor |
|---|---|
| CDS en `cds_from_genomic.fna` | 1 216 |
| Pseudogenes | 50 |
| Locus con copia de secuencia idéntica | 54 (27 pares; salmon conserva uno de cada par) |
| En la cuantificación de salmon | 1 189 |
| Analizables en DESeq2 | 1 186 |
| Universo del enriquecimiento | 1 138 |

**Identidad de la flagelina.** RS05045 = locus antiguo BARBAKC583_1040 = **flaA,
"Flagellin A"** en STRING. RS05445 = BARBAKC583_1120 es la proteína 3 asociada
al gancho (HAP3) en STRING; UniProt (A1UTS8) la rotula "Flagellin" por familia.
Las cifras de "flagelina" de este documento se refieren a RS05045.

**Top 20, 37 °C frente a 25 °C** (`resultados/genes/top20_37C_vs_25C.csv`,
padj < 0,05, ordenados por log2 del cambio):
- Reprimidos: 6 de los 20 son flagelares (FliG RS05635, FliM/FliN RS05620, MotA,
  FliQ, flagelina RS05045, FlgA). Aparecen también genes de fago (cápside mayor
  N4-gp56, proteína portal).
- Inducidos: 5 chaperonas o proteasas de choque térmico (ClpB, Hsp20, GroES, HtpX,
  GroEL); los dos primeros son un sistema toxina-antitoxina (AbrB/MazE, VapC).

### STRING (`10_descargar_string.sh`, `11_analisis_string.R`)

Versión **12.0** (`https://version-12-0.string-db.org`), descargada el
2026-09-26 y guardada en `datos/string/`. Identificadores: locus antiguos
(360095.BARBAKC583_xxxx). 1 147 locus consultados, 1 142 con identificador.
Universo con identificador: **1 115 de 1 138** (23 sin locus antiguo o sin
entrada). En ese universo: 275 reprimidos y 197 inducidos.

Enriquecimiento local (hipergeométrico, BH dentro de cada categoría, mismo
universo que la sección 6):

| Término (reprimidos) | Genes | Reprimidos | Esperados | q |
|---|---|---|---|---|
| KEGG bbk02040, Flagellar assembly | 25 | **25** | 6,2 | **1,9 × 10⁻¹⁴** |
| Palabra clave KW-0282, Flagellum | 26 | **26** | 6,4 | **3,3 × 10⁻¹⁵** |
| KEGG bbk00190, Oxidative phosphorylation | 34 | 20 | 8,4 | 6,1 × 10⁻⁴ |

Inducidos: ningún término con q < 0,05 (STRING no agrupa las chaperonas en un
término propio; nuestra categoría Chaperona sí es significativa, sección 6).

**Concordancia con KEGG.** Los 25 genes de bbk02040 están **todos** en el
conjunto flagelar de 31. El conjunto añade 6 con nombre flagelar que KEGG no
incluye en esa vía: FliO (RS04185), FliL (RS05555), FlaF, FlbT, FliK y MotC.

**Módulo en la red** (puntaje combinado ≥ 0,4): **400 aristas** entre los 31
genes flagelares; 10 000 conjuntos aleatorios del mismo tamaño dan en promedio
13,8 y como máximo 61 (p empírico = 1 / 10 001).

Red de los genes con cambio mayor a 2 veces: 195 genes con identificador, 697
aristas, 32 genes sin ninguna arista.

---

## 15. Figuras (Fase 3, 26 de septiembre de 2026)

`12_figuras.R` produce seis figuras en `figuras/` (PDF para imprenta y PNG a
300 ppp, 180 mm de ancho). Las cifras que muestran se recalcularon desde
`resultados/genes/tabla_de_genes.csv` y coinciden con `registros/12_figuras.log`.

| Figura | Contenido | Cifras de origen |
|---|---|---|
| 1 | Flujo de trabajo | — |
| 2 | PCA y lecturas asignadas a genes bacterianos, por muestra | secciones 3 y 4 |
| 3 | Volcanes de temperatura (37 frente a 25 °C) y de pH (6 frente a 8) | abajo |
| 4 | Genes flagelares y chaperonas: mapa de calor por réplica y curva térmica | abajo |
| 5 | Red STRING de los genes con cambio mayor a 2 veces | abajo y sección 14 |
| 6 | Ortólogos recíprocos de los 31 genes flagelares en el género; la cruz marca el acierto parálogo de FliI | sección 9 |

**Figura 3.** Se cuentan los transcritos analizables (1 186; sección 14), no el
universo deduplicado de la sección 6. En 1 184 hay valor q; los otros 2
(RS07210 y RS01255) tienen cero lecturas en las 12 muestras y DESeq2 no les
asigna valor q.

| Contraste | Reprimidos | Inducidos | Total (padj < 0,05) | Flagelares significativos |
|---|---|---|---|---|
| 37 °C frente a 25 °C | 283 | 205 | 488 | 31 de 31 |
| pH 6 frente a pH 8 | 187 | 188 | 375 | 6 de 31 |

**Figura 4.** 45 genes: los 31 flagelares y 14 chaperonas. Las chaperonas son la
categoría "Chaperona" de `scripts/categorias.R` (15 genes, definida por
palabras clave antes de ver los resultados) sin FlgA (RS05575), que es
flagelar. La categoría incluye chaperonas que no son de choque térmico (Hfq,
SecB, PCu(A)C, HemW, la chaperona del sistema de secreción tipo III, ATP12),
además de GroEL, GroES, DnaK, DnaJ, una proteína con dominio DnaJ, Hsp33, ClpB y
HtpX.
Hsp20 (RS02960, "Hsp20 family protein"), RpoH y Lon, que suben con la
temperatura (sección 8), no entran: su anotación no contiene las palabras clave.

- Mapa de calor: puntuación Z por gen de `vst(blind = FALSE)` del modelo de
  12 muestras, en las dos réplicas de Pl25, Pl30 y Pl37.
- Curva térmica: log2 del cambio respecto a 25 °C (contrastes 30 frente a 25 y
  37 frente a 25). Medianas:

| Grupo | 25 °C | 30 °C | 37 °C |
|---|---|---|---|
| Flagelares (n = 31) | 0 | −1,594 | −1,293 |
| Chaperonas (n = 14) | 0 | +0,394 | +0,839 |

A 30 °C, 30 de los 31 flagelares tienen log2 del cambio negativo, y 3 de las 14
chaperonas.

**Figura 5.** De los 195 genes con cambio mayor a 2 veces e identificador en
STRING (sección 14), se dibujan los **163 con al menos una arista**; **697
aristas** con puntaje combinado ≥ 0,4. Disposición de Fruchterman-Reingold con
semilla 42, idéntica en cada corrida. La posición de algunas etiquetas (ggrepel
0.9.8) puede variar entre corridas en la misma computadora: en 12 corridas
seguidas salieron dos versiones (5 y 7), con la misma disposición de la red.
Solo cambia el dibujo de las etiquetas, no los datos.
