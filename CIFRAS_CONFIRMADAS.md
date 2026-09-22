# Cifras confirmadas

Regeneración completa del análisis el 22 de septiembre de 2026, en carpeta
limpia (`v2/`), con comparación contra la versión anterior.
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
| CDS anotadas en el GFF | **1 217** |
| CDS con `protein_id` | 1 166 |
| Accesiones de proteína únicas | **1 139** (coincide con `protein.faa` de NCBI) |
| Accesiones duplicadas (región de ~28 kb) | 27 |
| Genes con anotación flagelar | **31** |

El conjunto flagelar se define en `scripts/01_definicion_conjunto_flagelar.sh`
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

Las muestras en sangre y en células son mayoritariamente ARN humano: solo
0,2–6 % de sus lecturas son bacterianas.

| Transcritos | Valor |
|---|---|
| En el índice de salmon | 1 189 |
| Analizables (`rowSums > 0`) | **1 186** |

---

## 4. Estructura de los datos (PCA)

`vst(blind=TRUE)`, `plotPCA(ntop=500)`.

- **PC1 = 53 %** de la varianza · **PC2 = 16 %**
- Correlación PC1 ~ tasa de mapeo: **0,817**
- Correlación PC2 ~ tasa de mapeo: 0,06

PC1 separa por **matriz biológica**, que determina simultáneamente qué
fracción de la librería es bacteriana. Ambos factores son inseparables en
este diseño.

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

## 8. Perfil térmico: dos programas independientes

De los 31 genes flagelares, el contraste **37 vs 30 °C no es significativo en
30 de 31** (único con p<0,05: FliF RS05510, p = 0,012, que no sobrevive la
corrección por 31 pruebas). **23 de 31 alcanzan su mínimo a 30 °C.**

Conteos normalizados (media por condición):

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

**La represión flagelar ocurre entre 25 y 30 °C** (flagelina −71 %), donde las
chaperonas apenas se mueven (GroEL +5 %). **La respuesta de choque térmico
ocurre entre 30 y 37 °C** (GroEL ×4,6, ClpB ×8,3), donde el flagelo ya no
cambia. Son dos programas con umbrales térmicos distintos; el flagelar
precede al de choque térmico.

Lectura correcta: la expresión flagelar es **máxima a 25 °C** y ya está
reprimida a 30 °C, sin descenso adicional a 37 °C.

---

## 9. Genómica comparada

### Panel de género (BLASTp, e<1e-5, de los 31 flagelares)

| Especie | Genes presentes |
|---|---|
| *B. bacilliformis* | 31 / 31 |
| *B. clarridgeiae* | 31 / 31 |
| *B. schoenbuchensis* | 31 / 31 |
| *B. ancashensis* | 31 / 31 |
| *B. tribocorum* | **2 / 31** |
| *B. quintana* | **2 / 31** |
| *B. henselae* | **2 / 31** |

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
| 1 216 CDS | **1 217** | Faltaba BARBAKC583_RS03925 |
| 1 176 analizables | **1 186** | Recuento anterior no reproducible |
| Universo 1 134 (o 1 175) | **1 138** | `lt2prot.tsv` anterior incluía 50 CDS sin `protein_id` |
| 29 de 31 flagelares | **31 de 31 significativos** (30 repr. + 1 ind.) | |
| 216 transcritos modificados | **207** (restringido, \|LFC\|>1) o **488** (padj<0,05) | Definición anterior no rastreable |
| q = 1,3 × 10⁻¹⁶ | **2,14 × 10⁻¹⁷** | |
| NES −2,46 · q 1,8 × 10⁻⁹ | **−2,44 · 2,43 × 10⁻⁹** | |
| 6 flagelares con umbral formal | **7** | |

---

## 12. Reproducibilidad

Todas las cifras de este documento se regeneran con los scripts de `v2/`:

| Script | Produce |
|---|---|
| `s01_anotacion.sh` | secciones 2 |
| `s02_deseq.R` | secciones 3, 4 |
| `s03_enriquecimiento.R` | sección 6 |
| `s04_robustez.R` | secciones 5, 7, 8 |
| `s05_pca_mapeo.R` | sección 4 |

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
```
