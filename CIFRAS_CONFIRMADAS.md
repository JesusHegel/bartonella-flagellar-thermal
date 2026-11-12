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

Las 23 corridas fueron depositadas por The University of Montana y
secuenciadas en Illumina HiSeq 2500, pareado de 2 × 150 pb.

**Corrección (12 de noviembre de 2026).** Este documento afirmaba antes que las
corridas provenían de **un solo centro** y que **no existían dos centros de
secuenciación**. Es falso. Wachter et al. 2020 indica en Materiales y métodos
que el ARN se envió al **Yale Center for Genomic Analysis** (Pl25, Pl30, Pl37,
pH06, pH07, pH08) o a **GENEWIZ** (PlBG, HUVE, HB37, HBBG). El campo de centro
de los metadatos de SRA registra quién deposita, no quién secuencia; de ahí
salió el error.

El corte entre los dos centros **coincide exactamente** con el corte entre las
12 muestras del modelo y las 11 excluidas. Dos consecuencias: el modelo
restringido queda libre de efecto de lote por centro, y la exclusión de las 11
gana una razón independiente de la matriz biológica y del PCA. Explica además
el hallazgo de adaptadores de la sección 3b.

### Diseño del estudio original

KC583, pases 4 a 7, crecida **cuatro días en placa HIBB a 30 °C** (agar
infusión corazón Bacto con 4 % de sangre desfibrinada y 2 % de suero, ambos de
oveja) y **trasladada 2 h** a cada condición. Tabla 1 del artículo.

| Condición | Medio | Traslado |
|---|---|---|
| Pl25, Pl30, Pl37 | placa HIBB | 2 h |
| PlBG | placa HIBB, atmósfera sanguínea | 2 h |
| pH06, pH07, pH08 | HIBB líquido | 2 h |
| HB37, HBBG | sangre humana | 2 h |
| HUVE | HUVEC en medio EGM-Plus | 24 h |

**Ninguna condición es crecimiento sostenido a esa temperatura.** Lo que mide
este reanálisis es la respuesta a un traslado de dos horas. Escribir "cultivos
crecidos a 25, 30 y 37 °C" sería incorrecto.

Wachter justifica las dos horas: con **cuatro horas a 37 °C** hubo muerte
bacteriana, que atribuye a la adaptación de KC583 al laboratorio.

### Réplicas descartadas por los autores originales

Se enviaron **tres réplicas biológicas por condición** (30 bibliotecas) y en el
SRA hay **23**. La tabla S2 del artículo lista las 30 y marca en amarillo las
eliminadas, con la leyenda "Replicates removed from global analysis". Son siete
y coinciden una a una con las que faltan en el depósito.

| Réplica | Alineadas | Ambiguas |
|---|---|---|
| pH08-1 | 12,4 % | 60,7 % |
| pH06-3 | 21,0 % | 44,4 % |
| HUVE-3 | 2,1 % | 81,5 % |
| Pl25-3 | 99,3 % | 56,6 % |
| Pl30-1 | 99,1 % | 31,1 % |
| Pl37-1 | 99,2 % | 33,6 % |
| pH07-3 | 91,0 % | 2,5 % |

Criterio declarado: descartaron las réplicas que **no correlacionaban** con las
otras de su condición, no las técnicamente defectuosas. Cuatro de las siete no
tienen problema técnico, y pH07-3 es la más limpia de sus tres hermanas.

**Consecuencia para el modelo.** Las 12 muestras son dos réplicas por
condición, y en cada una son las dos sobrevivientes de ese filtro. En Pl25 y
Pl30 se conservó la réplica con peor depleción ribosomal (80,3 % y 84,7 % de
ambiguas) y se descartó la mejor (56,6 % y 31,1 %); esas dos conservadas son
las mismas que la sección 3 marca por debajo de dos millones de lecturas
asignadas.

El filtro **no se puede revertir ni reproducir**, porque el umbral de
correlación no está publicado. Va declarado en limitaciones. El respaldo está
en la sección 7: los 31 genes conservan la dirección en las dos particiones
independientes, es decir con cualquiera de las dos réplicas de cada
condición.

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

## 3b. Control de calidad de las lecturas (2 y 3 de noviembre de 2026)

`fastqc` 0.12.1 y `multiqc` sobre las **lecturas sin recortar**, que son las que
entran a salmon. 46 informes, dos por corrida. Los produce
`scripts/00_descargar_y_cuantificar.sh` cuando `fastqc` está en el PATH; van en
un entorno aparte (`entorno/control_de_calidad.yml`) para no obligar a conda a
resolver de nuevo el entorno del análisis, que ya está verificado (12e).
**Ninguna cifra del análisis depende de este entorno**: fastqc solo describe las
lecturas. Salida en `control_de_calidad/`, fuera de `resultados/`, así que
`verificar.sh` no la compara.

### Hallazgo: el contenido de adaptadores separa exactamente las 12 muestras del modelo

| Módulo de fastqc | 12 del modelo (24 archivos) | 11 excluidas (22 archivos) |
|---|---|---|
| Contenido de adaptadores | **24 pasan** | **16 fallan, 6 advierten** |
| Estadística básica | 24 pasan | 22 pasan |
| Calidad por base | 24 pasan | 22 pasan |
| Calidad por secuencia | 24 pasan | 22 pasan |
| Contenido de N | 24 pasan | 22 pasan |
| Secuencias sobrerrepresentadas | 8 pasan, 16 advierten | 15 advierten, 7 fallan |
| Distribución de longitud | 24 advierten | 2 pasan, 20 advierten |
| Contenido por base | 24 fallan | 1 advierte, 21 fallan |
| GC por secuencia | 1 pasa, 8 advierten, 15 fallan | 1 pasa, 8 advierten, 13 fallan |
| Niveles de duplicación | 1 advierte, 23 fallan | 22 fallan |

El contenido de adaptadores es **el único módulo que reproduce la partición de
12 frente a 11 usada en el modelo**, y lo hace sin una sola excepción: los 22
fallos y advertencias caen en los 22 archivos de las 11 muestras excluidas.

**Explicación (12 de noviembre de 2026).** No es coincidencia. Las 12 del
modelo y las 11 excluidas se prepararon en **centros distintos**: Yale Center
for Genomic Analysis y GENEWIZ respectivamente (sección 1). La partición de
adaptadores es el efecto de dos preparaciones de biblioteca independientes. Las
excluidas difieren además en el tratamiento del ARN: HUVE, HB37 y HBBG se
enriquecieron con MICROBEnrich, y las de sangre se deplecionaron de globina con
GLOBINclear.

Las otras dos medidas que se usaron para excluir no dan esa partición:

| Medida | Qué separa |
|---|---|
| Asignación a CDS (sección 3) | 8 muestras (sangre 4,6 %, endotelio 0,3 %). PlBG tiene 45,5 %, **más alta** que los cultivos (38,3 %) |
| PC1 (sección 4) | PlBG queda intermedio (−3,4) entre cultivos (+14,7) y sangre (−23,2) |
| Contenido de adaptadores | Las 11, incluidas las tres PlBG |

PlBG se había excluido **solo por diseño**, por no pertenecer a las series de
temperatura ni de pH (sección 3). El contenido de adaptadores le da a esa
exclusión un respaldo medido que antes no tenía.

### Lo que este control no dice

**Las 12 muestras del modelo no pasan todos los módulos.** Fallan en contenido
por base (24 de 24), niveles de duplicación (23 de 24) y GC por secuencia (15 de
24). Los tres son los falsos positivos conocidos de RNA-seq: cebado no aleatorio
al inicio de la lectura, transcritos muy abundantes que fastqc cuenta como
duplicados, y una distribución de GC que no es unimodal porque el
transcriptoma no tiene el GC medio del genoma. **No puede escribirse que las
muestras del modelo "pasan el control de calidad"**, sino que pasan los módulos
de calidad de base, calidad de secuencia, contenido de N y contenido de
adaptadores.

Pendiente de verificar: los umbrales exactos con que fastqc 0.12.1 decide
`pass`, `warn` y `fail` en el módulo de adaptadores. Hay que citarlos de la
documentación de la herramienta antes de afirmar un porcentaje en el
manuscrito.

### Contenido de GC por condición

| Condición | Archivos | %GC |
|---|---|---|
| Pl25 | 4 | 39 – 41 |
| Pl30 | 4 | 41 – 42 |
| Pl37 | 4 | 42 – 43 |
| pH06 | 4 | 43 – 44 |
| pH07 | 4 | 45 – 46 |
| pH08 | 4 | 43 – 45 |
| PlBG | 6 | 41 – 42 |
| HB37 | 6 | 51 – 52 |
| HBBG | 6 | 51 – 52 |
| HUVE | 4 | 55 – 57 |

Genoma de referencia: **38,24 %** de GC (calculado sobre
`datos/referencia/GCF_000015445.1_ASM1544v1_genomic.fna`, 1 445 021 nt).

El GC separa **8 muestras**, las de sangre y endotelio, no las 11. Las tres de
PlBG tienen GC de cultivo (41 a 42 %) y aun así fallan en adaptadores: son dos
señales distintas y no deben presentarse como una sola. El exceso de GC en
sangre y endotelio es compatible con material del hospedero, que tiene GC más
alto, pero esto **no se ha verificado** asignando lecturas al genoma humano.

### Consecuencia para los métodos

La decisión de no recortar las lecturas queda respaldada por medición **para las
12 muestras del modelo**, que son las únicas que entran a los contrastes. Las 11
excluidas sí habrían necesitado recorte, lo que agrega una razón independiente a
su exclusión.

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
| un solo centro | **Yale n=12 vs GENEWIZ n=11** | Wachter et al. 2020, Materiales y métodos. El campo de centro de SRA registra quién **deposita** (University of Montana), no quién secuencia. Esta fila estaba invertida: retiraba el dato correcto. Restituida el 12-11-2026 |
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

## 12d. Réplica desde el SRA (nivel 2, 26 de septiembre de 2026)

`scripts/replica_nivel_2.sh` en la computadora de referencia: descarga de las
23 corridas del SRA, cuantificación con salmon (índice con las seis huellas del
original), análisis completo en carpetas `_desde_sra` y comparación con
`scripts/comparar_nivel_2.R`. Descarga y cuantificación: 40 minutos, sin fallos.

**Cuantificación.** Lecturas procesadas idénticas en las 23 corridas. Lecturas
asignadas: diferencia máxima de 0,054 %. Correlación mínima de NumReads entre
las dos cuantificaciones (log2): 0,9986.

**Cifras principales** (publicadas, desde el SRA):

| Cifra | Publicada | Desde el SRA |
|---|---|---|
| 37 vs 25: significativos | 488 | 491 |
| 37 vs 25: significativos y \|LFC\|>1 | 207 | 206 |
| 37 vs 25: flagelares significativos (reprimidos + inducidos) | 31 (30 + 1) | 31 (30 + 1) |
| Modelo de 23 muestras: significativos / flagelares | 538 / 30 | 532 / 30 |
| 30 vs 25: significativos / flagelares | 423 / 30 | 422 / 30 |
| 37 vs 30: flagelares no significativos | 30 | 30 |
| pH 6 vs 8: significativos / flagelares | 375 / 6 | 375 / 6 |
| Hipergeométrico Flagelar: reprimidos · q | 30 · 2,14 × 10⁻¹⁷ | 30 · 2,39 × 10⁻¹⁷ |
| Hipergeométrico Chaperona: q de inducidos | 0,0117 | 0,0125 |
| fgsea Flagelar: NES · q | −2,441 · 2,43 × 10⁻⁹ | −2,444 · 1,77 × 10⁻⁹ |
| Umbral formal lfcThreshold=1: significativos / flagelares | 49 / 7 | 50 / 7 |
| STRING KEGG bbk02040: reprimidos · q | 25 · 1,94 × 10⁻¹⁴ | 25 · 2,13 × 10⁻¹⁴ |

**Todas las cifras de los genes flagelares son idénticas.** Los conteos del
genoma completo varían en pocos genes y los valores q conservan su orden de
magnitud. Genes que cambian de significancia: 3 en 37 vs 25, 1 en 30 vs 25,
0 en 37 vs 30 y en pH, 8 en el modelo de 23 muestras; **ninguno flagelar**.
Todos tenían padj entre 0,043 y 0,055, salvo RS00085 (proteína hipotética) en
el modelo de 23 muestras (0,046 → 0,14): sus lecturas cambian en dos corridas
PlBG (650 → 1 069 y 541 → 685), que no están en el modelo de 12 muestras. Los
valores fraccionarios indican lecturas compartidas con otra secuencia, que
salmon reparte de forma distinta en cada corrida con varios hilos.

**Las cifras del manuscrito son las publicadas** (cuantificaciones incluidas en
`cuantificacion/`); desde el SRA se reproducen todas las conclusiones.

## 12e. Réplica independiente en una segunda computadora (17 y 24 de octubre de 2026)

Nivel 1 ejecutado por B. Castañeda en un equipo distinto al de referencia,
siguiendo únicamente `GUIA_DE_REPLICA.md`.

| | Computadora de referencia | Réplica independiente |
|---|---|---|
| Procesador | AMD Ryzen 5 PRO 4650G | Intel Core i7-13620H (13.ª gen.) |
| Sistema | Debian GNU/Linux 13 (trixie) | Ubuntu 26.04.1 LTS sobre WSL 2 |

Entorno creado en ambas desde `entorno/bartonella.yml` (salmon 2.4.1, R 4.5.3).
Son dos arquitecturas de fabricantes distintos, que es la condición en la que
esta comprobación tiene valor.

**Resultado: 28 de las 39 tablas idénticas byte a byte y 11 equivalentes dentro
de la tolerancia. Ninguna distinta, ninguna faltante.**

Los cinco contrastes de DESeq2 conservan exactamente los mismos conteos de genes
significativos que las tablas publicadas: 488 (37 vs 25), 538 (modelo de 23
muestras), 423 (30 vs 25), 128 (37 vs 30) y 375 (pH 6 vs 8). **Ningún gen cambia
de significancia ni de signo.** Las comprobaciones de significancia añadidas en
fd162ff confirman además 3 vías significativas en fgsea, 1 categoría inducida y
2 reprimidas en el hipergeométrico, y 0 y 46 términos de STRING en inducidos y
reprimidos.

Desviación máxima de toda la corrida: 1,0 × 10⁻³ absoluta y 2,1 × 10⁻³ relativa,
en `genes/tabla_de_genes.csv`. En los contrastes, como máximo 4,1 × 10⁻⁶
absoluta. `universo_genes.csv` resultó idéntico byte a byte, y los cuatro
archivos de `datos/string/` conservan su md5: la réplica usó exactamente los
mismos datos de STRING.

**Hallazgo: tres tablas mal clasificadas por el verificador.** En la primera
corrida (17-10-2026, antes de fd162ff), `enriquecimiento_hipergeometrico.csv` y
las dos tablas de enriquecimiento de STRING salieron `DISTINTO`. Medidas con
`comparar_con_tolerancia.R` resultaron equivalentes, con diferencias relativas
máximas de 3,3 × 10⁻¹⁵, 3,0 × 10⁻¹⁵ y 5,4 × 10⁻¹⁵, y ningún término cruzando el
umbral de 0,05. En el hipergeométrico los conteos son idénticos (31 genes
flagelares, 1 inducido, 30 reprimidos, 7,6 esperados) y solo varían los últimos
dígitos de `p_rep` y `q_rep`.

Esas tres tablas no estaban en la lista de comparación con tolerancia de
`verificar.sh`, pese a contener valores p y q sujetos a la misma variación entre
procesadores descrita en 12c. Y como no tienen columna `padj`, tampoco pasaban
por la comprobación de cambios de significancia: la tolerancia numérica por sí
sola habría aceptado un término cruzando el 0,05 sin avisar. Ambas cosas se
corrigieron en fd162ff, que añade las tres tablas a la lista y extiende la
comprobación de significancia a las columnas `q`, `q_ind` y `q_rep`. **El
verificador queda más estricto que antes, no menos.** La corrida del 24-10-2026
con el parche aplicado es la que da el resultado de arriba.

El hueco solo podía aparecer en una segunda computadora: en la de referencia
esas tres tablas salen idénticas a sí mismas y nunca llegan a compararse con
tolerancia.

## 12f. Segunda réplica desde el SRA (nivel 2, 2 y 3 de noviembre de 2026)

Descarga y cuantificación de las 23 corridas desde el SRA, con el informe de
calidad activado (3b), y análisis completo con `scripts/replica_nivel_2.sh`.
Registros en `registros/registro_calidad.txt`,
`registros/registro_calidad_2.txt` y `registros/registro_nivel2.txt`.

Equipo: AMD Ryzen 5 PRO 4650G with Radeon Graphics, Debian GNU/Linux 13 (trixie). Es **la misma computadora de referencia**
de 12c, 12d y 12e, con el mismo entorno conda. La única diferencia con 12d es
que esta corrida volvió a descargar y cuantificar desde cero, con el control de
calidad activado.

Las seis huellas del índice de salmon coinciden con el original. Descarga,
cuantificación y control de calidad: de 01:38 a 02:59, una hora y 21 minutos
para 22 corridas (la primera venía de un intento anterior). Análisis completo:
48 segundos.

**Cuantificación.** Lecturas procesadas idénticas en las 23 corridas. Diferencia
máxima en lecturas asignadas: **0,07 %** (SRR12653239). Correlación mínima de
NumReads entre las dos cuantificaciones (log2): **0,992541** (SRR12653251, que
es además la que más transcritos cambia, 150).

**Cifras principales:**

| Cifra | Publicada | Desde el SRA |
|---|---|---|
| 37 vs 25: significativos | 488 | 490 |
| 37 vs 25: significativos y \|LFC\|>1 | 207 | 206 |
| 37 vs 25: flagelares significativos (reprimidos + inducidos) | 31 (30 + 1) | 31 (30 + 1) |
| Modelo de 23 muestras: significativos / flagelares | 538 / 30 | 538 / 30 |
| 30 vs 25: significativos / flagelares | 423 / 30 | 424 / 30 |
| 37 vs 30: flagelares no significativos | 30 | 30 |
| pH 6 vs 8: significativos / flagelares | 375 / 6 | 373 / 6 |
| Universo: inducidos / reprimidos | 200 / 278 | 202 / 278 |
| Hipergeométrico Flagelar: reprimidos · q | 30 · 2,14 × 10⁻¹⁷ | 30 · 2,14 × 10⁻¹⁷ |
| Hipergeométrico Chaperona: q de inducidos | 0,0117 | 0,0125 |
| fgsea Flagelar: NES · q | −2,441 · 2,43 × 10⁻⁹ | −2,447 · 2,69 × 10⁻⁹ |
| Umbral formal lfcThreshold=1: significativos / flagelares | 49 / 7 | 50 / 7 |
| STRING KEGG bbk02040: reprimidos · q | 25 · 1,94 × 10⁻¹⁴ | 25 · 1,94 × 10⁻¹⁴ |

Conteos iguales en **12 de 18**; diferencia relativa máxima en valores q y NES:
**10,7 %** (q de fgsea Flagelar). **Las ocho cifras de genes flagelares son
idénticas**, y los dos valores q de las pruebas principales (hipergeométrico
2,14 × 10⁻¹⁷ y STRING 1,94 × 10⁻¹⁴) coinciden exactamente.

### Genes que cambian de significancia

**Ninguno flagelar, en ninguno de los cinco contrastes.**

| Contraste | Cambian | De borde | Fuera de borde |
|---|---|---|---|
| 37 vs 25 | 8 | 4 | 4 |
| 30 vs 25 | 11 | 9 | 2 |
| 37 vs 30 | 2 | 2 | 0 |
| pH 6 vs 8 | 4 | 4 | 0 |
| Modelo de 23 muestras | 12 | 6 | 6 |
| **Total** | **37** | **25** | **12** |

"De borde" significa padj entre 0,04 y 0,06 en ambas corridas: el gen no cambia
de valor, solo cruza el umbral de 0,05. Correlación de log2FC entre las dos
corridas: 0,99955 a 0,99999 según el contraste.

Los 12 cambios que no son de borde se concentran en **siete genes**, todos en la
mitad baja de expresión (percentiles 9 a 42 del TPM mediano):

| Locus | Producto | Contrastes | Factor máximo de cambio en padj |
|---|---|---|---|
| RS00030 | DNA polymerase III subunit epsilon | 3 | **1 200 ×** (0,12 → 1 × 10⁻⁴) |
| RS00015 | DNA polymerase I | 2 | 36 × (0,0014 → 0,051) |
| RS04140 | proteína hipotética (204 nt) | 1 | 27 × (0,21 → 0,0077) |
| RS01760 | ComF family protein | 1 | 10 × (0,043 → 0,41) |
| RS00040 | shikimate dehydrogenase | 2 | 4 × |
| RS00035 | dephospho-CoA kinase | 1 | 3 × |
| RS00065 | tRNA uridine-5-carboxymethylaminomethyl… | 1 | 2 × |

**RS00030 no es utilizable**: su padj se mueve tres órdenes de magnitud entre
corridas de los mismos datos, en tres contrastes distintos. Lo mismo vale, en
menor grado, para los otros seis.

Mecanismo: salmon con varios hilos no es determinista y reparte de forma
distinta las lecturas compartidas entre secuencias parecidas (12c, 12d). Se
comprobó que **no** es la región duplicada de ~28 kb (sección 2), que abarca
RS03130 a RS03445: ninguno de los siete genes tiene copias idénticas. Cuatro de
los siete (RS00015, RS00030, RS00035, RS00040) son consecutivos al inicio del
cromosoma, entre las posiciones 1 970 y 7 360, y RS00035 y RS00040 se solapan en
4 nt. La causa precisa **no se ha establecido**.

### Comparación con la primera réplica de nivel 2 (12d)

Las dos corridas de nivel 2 se hicieron **en la misma computadora, con el mismo
entorno y los mismos datos de origen**. Difieren solo en que son dos descargas y
dos cuantificaciones distintas. Aun así:

| | 12d (26 de septiembre) | 12f (3 de noviembre) |
|---|---|---|
| 37 vs 25: significativos | 491 | 490 |
| Modelo de 23 muestras: significativos | 532 | 538 |
| 30 vs 25: significativos | 422 | 424 |
| pH 6 vs 8: significativos | 375 | 373 |
| Genes que cambian de significancia | 3, 1, 0, 0, 8 | 8, 11, 2, 4, 12 |

**El número de genes que cruzan el umbral de 0,05 no es estable entre réplicas**,
ni siquiera en la misma máquina: 12 cambios en 12d frente a 37 en esta. Esto
aísla la causa, porque el equipo, el sistema, el entorno y los datos de entrada
son idénticos: lo único que varía es el reparto de lecturas que hace salmon con
varios hilos.

Lo que sí es idéntico en las tres corridas (publicada, 12d y 12f) son **todas
las cifras de genes flagelares**, y **ningún gen flagelar cambia de
significancia en ninguna**.

**Las cifras del manuscrito son las publicadas** (cuantificaciones incluidas en
`cuantificacion/`); desde el SRA se reproducen todas las conclusiones.

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

---

## 16. Candidatos publicados por otros estudios (3 de noviembre de 2026)

`scripts/13_candidatos_publicados.sh` ubica por BLASTp en el proteoma de KC583
proteínas identificadas en otros trabajos, como blanco vacunal o como objeto de
estudio funcional, y les añade los valores de expresión de este reanálisis.
Salida en `resultados/candidatos/`. Determinista: dos corridas seguidas dan md5
idéntico en las tres tablas, con `-num_threads 1`.

**67 consultas de siete estudios** (cerrado el 12 de noviembre de 2026):
5 de Gul et al. 2020, 18 de Rivera-Asencios et al. 2025 (Khan et al. 2020
propone FlhA, la misma proteína del constructo II de Rivera-Asencios), 7 de
Jimenez-Vasquez et al. 2023, 4 de Gomes et al. 2016, 22 de Dichter et al. 2021,
10 de Padilla Rojas et al. 2019 y la secuencia curada de IalB. Asignadas por
secuencia: 66; ambigua: 1 (A1USN1).

Henriquez-Camacho 2015, que la introducción también cita, es una revisión y no
propone candidatos propios: discute flagelina, IalB y Pap31, los tres ya
incluidos. **El cruce cubre todos los estudios que proponen blancos.**

### Deduplicación por locus

Siete consultas caen en loci que otro estudio ya había propuesto y se colapsan
antes del contraste, porque este cuenta genes distintos y no consultas: si tres
trabajos proponen Pap31, Pap31 pesa una vez. No flagelares 55 a 48.

### Cobertura del blanco

`qcovs` mide cobertura de la **consulta**. Cuando la consulta es más larga que
el blanco, la cobertura de consulta baja aunque el blanco esté cubierto entero.
El control original rechazaba ese caso, que es el seguro: lo peligroso es una
consulta **recortada**, que da 100 % de cobertura faltándole media proteína.
Ahora se calcula también la cobertura del blanco y se acepta la asignación si
el blanco está cubierto al 95 % o más con identidad suficiente.

Un caso lo motivó: `dic15_autoC`. Dichter cita WP_193741240.1 declarando 297 aa,
pero NCBI extendió esa entrada a 416 aa después de la publicación. El acierto en
KC583 es WP_052467648.1, de **297 aa exactos y 100 % de identidad**, con el
segundo acierto a 68,2 %. La discrepancia está en la base de datos, no en la
asignación. Se incluye, y su gen RS05525 está reprimido, de modo que incluirlo
**debilita** el contraste: es la decisión conservadora.

### Por qué el emparejamiento es por secuencia y no por nombre

La anotación de RefSeq repite **84 nombres de producto** en este genoma, entre
ellos seis `invasion associated locus B family protein`, seis
`autotransporter outer membrane beta-barrel domain-containing protein` y seis de
cada clase de transportador ABC. Ninguno de los 31 genes flagelares comparte su
nombre con otro gen, pero varios candidatos no flagelares sí, y otros caen
dentro de operones de transportadores donde genes vecinos llevan nombres
parecidos y distintos.

Dos de las 24 asignaciones tienen `n_genes_mismo_producto` igual a 6, es decir
que por nombre eran indeterminables: `ialB_publicado` y `A1UTU0`.

Columnas de control de `candidatos_publicados.tsv`: `asignacion` (solo
`por secuencia` es utilizable), `ambiguo`, `segundo_identidad_pct`,
`n_genes_mismo_producto` y las dos longitudes, porque BLASTp mide cobertura de
la consulta y una consulta recortada da 100 % aunque le falte media proteína.

### Procedencia de las secuencias

Los 23 accesos de UniProt **ya fueron retirados de UniProtKB**:
`rest.uniprot.org/uniprotkb/<acceso>.fasta` responde 200 con cuerpo vacío
(`content-length: 0`), que se confunde con un fallo de red. Solo quedan en
UniParc. Obtenidos el 3 de noviembre de 2026, versión de UniProt 2026_03, con
el `UPI` y la longitud de cada uno en
`datos/candidatos_publicados/procedencia_uniparc.txt`.

La secuencia de IalB es **P35641.1** (`IALB_BARBK`, 186 aa), entrada curada de
UniProtKB/Swiss-Prot obtenida del NCBI. El mnemónico `BARBK` corresponde a
*B. bacilliformis* KC583, la misma cepa del genoma de referencia. No debe
confundirse con **P35640.1** (`IALA_BARBK`, 170 aa), que es IalA.

### Resultado principal: los ocho blancos flagelares están reprimidos

Siete proteínas del constructo vacunal II de Rivera 2022, que aporta sus 12
epítopos de proteínas asociadas a virulencia, más FliP de Gul et al. 2020:

| Proteína | Locus | Estudio | log2FC 37/25 | padj |
|---|---|---|---|---|
| FliG | RS05635 | Rivera 2022 | −2,405 | 3,11 × 10⁻¹⁰ |
| MotA | RS05615 | Rivera 2022 | −2,111 | 2,45 × 10⁻⁹ |
| FlhB | RS05640 | Rivera 2022 | −1,596 | 8,74 × 10⁻¹¹ |
| FlgK | RS05465 | Rivera 2022 | −1,360 | 4,46 × 10⁻⁶ |
| Proteína del cuerpo basal | RS05605 | Rivera 2022 | −1,293 | 1,95 × 10⁻⁵ |
| FlhA | RS05690 | Rivera 2022; Khan 2019 | −1,276 | 1,27 × 10⁻⁷ |
| FliP | RS05550 | Gul 2020 | −1,265 | 1,70 × 10⁻³ |
| FliF, anillo M | RS05510 | Rivera 2022 | −0,902 | 5,25 × 10⁻⁴ |

**Ocho de ocho reprimidos, todos con padj menor que 0,002.**

El constructo II de Rivera se formó filtrando por "proteínas asociadas a
virulencia" y conservación entre cepas. Sus 12 epítopos finales salen de esas
siete proteínas flagelares: el filtro devolvió flagelo en su totalidad.

### Control: el constructo I, de proteínas de membrana externa

El mismo trabajo armó un segundo constructo con un filtro distinto, proteínas
de membrana externa, sobre las mismas cepas y con el mismo método. Sirve de
control para saber si la represión es propia del flagelo o un efecto general de
la temperatura.

| Proteína | Locus | log2FC 37/25 | padj | |
|---|---|---|---|---|
| LptD | RS02440 | −0,828 | 4,02 × 10⁻⁴ | reprimida |
| BamA | RS02830 | −0,138 | 0,471 | sin cambio |
| BamD | RS04545 | −0,114 | 0,467 | sin cambio |
| Lipoproteína putativa | RS01620 | −0,009 | 0,968 | sin cambio |
| OMP85 | RS06535 | +0,090 | 0,749 | sin cambio |
| TolC | RS03155 | +0,534 | 0,0126 | inducida |
| Autotransportador | RS05535 | +0,688 | 3,54 × 10⁻⁵ | inducida |

Los tres epítopos experimentales, presentes en los dos constructos: UvrA
(RS03570, −0,148, padj 0,578), 4-hidroxibenzoato octaprenyltransferasa
(RS01890, −0,032, padj 0,830) y citrato sintasa (RS02885, −0,584, padj 0,0227).

### Especificidad

`resultados/candidatos/especificidad.tsv`, calculada por el mismo script sobre
las 66 consultas asignadas por secuencia, deduplicadas por locus:

| Grupo | Reprimidos a 37 °C | Total |
|---|---|---|
| Flagelares | 11 | 11 |
| No flagelares | 15 | 48 |

Fisher exacto de una cola: **p = 2,761 × 10⁻⁵**.

Los once flagelares incluyen ahora flagelina (RS05045, −2,007), FlgE (RS05470,
−0,704) y FlgH (RS05560, −1,414), que aporta Dichter.

### Los candidatos publicados se reparten en tres grupos

Con el cruce ampliado aparece por primera vez un candidato **inducido** a
temperatura de hospedero, lo que cambia la lectura del conjunto:

| Comportamiento a 37 °C | Candidatos |
|---|---|
| Reprimidos | los 8 flagelares, más FecCD, LptD, citrato sintasa, Prot 492 y Pap31 |
| Inducidos | **GroEL** (RS05735, +2,057, padj 3,11 × 10⁻¹⁰) y TolC, autotransportador, Prot 288 |
| Sin cambio | SCS-α (RS00130, −0,090), SCS-β (RS00125, −0,207), IalB y el resto |

**GroEL es el hallazgo del cruce ampliado.** Gomes et al. 2016 lo propone como
antígeno (candidato A, secuenciación N-terminal) y aquí sube cuatro veces a
37 °C. El locus RS05735 está anotado como `chaperonin GroEL` y ya figuraba en
el grupo de chaperonas de la figura 4 y en el top 20 del contraste de
temperatura: el cruce lo reencuentra por una vía independiente, partiendo de la
secuencia publicada y no de la anotación.

Pap31 merece una nota: UniProt da Q3HM50 como **fragmento** y de otra cepa, de
modo que la consulta usa Omp25 de KC583 (A1UU13, vía UniParc), que Gomes
reporta como su homólogo al 100 %. Su producto comparte nombre con otros tres
genes, así que por nombre habría sido indeterminable.

### Los antígenos que funcionan en pacientes son los que se inducen a 37 °C

Dichter et al. 2021 montó un ELISA con **porina B y autotransportador E** que
alcanza 81 % de sensibilidad y 95 % de especificidad con suero de pacientes
peruanos. Los dos están **inducidos** en estos datos:

| Antígeno del ELISA | Locus | log2FC 37/25 | padj |
|---|---|---|---|
| Porina B | RS05915 | **+1,591** | 1,55 × 10⁻¹² |
| Autotransportador E | RS05535 | **+0,688** | 3,54 × 10⁻⁵ |

Es el respaldo más fuerte del criterio que propone este trabajo, y viene de
datos ajenos ya validados experimentalmente. Pap31, el antígeno clásico, es
vecino de porina B en el genoma (RS05920 frente a RS05915) y va al revés:
−1,615.

### FlgH: reprimida y reconocida por los pacientes

De las tres proteínas flagelares del panel de Dichter, flagelina y FlgE fueron
reactivas solo con **suero de conejo**; **FlgH** está entre las catorce que
reconoce el suero de pacientes. Y FlgH está reprimida: RS05560, −1,414, padj
9,93 × 10⁻⁶.

Ese caso solo impide cualquier lectura simple del tipo "reprimida, luego no
sirve como antígeno". Una proteína puede estar reprimida en cultivo a 37 °C y
aun así generar anticuerpos en el paciente.

### RS05650 es un candidato publicado

`pad09_ompHipC` cae en **RS05650** (−1,652, padj 1,26 × 10⁻⁸), el gen que está
dentro del clúster flagelar, entre FlhB (RS05640) y FliR (RS05685), anotado como
`hypothetical protein` y por eso fuera de los 31.

Resulta ser una de las proteínas de origen de los epítopos de M1 de Padilla
Rojas. Que el regulón sea más amplio que la anotación deja de ser una
observación interna y pasa a tener consecuencia sobre un candidato vacunal
publicado. Su nombre de producto lo comparten **98 genes**: por nombre era
imposible.

**Este es un contraste dentro del conjunto de candidatos publicados**, que es
pequeño y no fue elegido al azar: las proteínas provienen de tamizajes con
criterios propios y no son independientes entre sí. No sustituye a la prueba
genómica, que es el hipergeométrico de la sección 6
(q = 2,14 × 10⁻¹⁷). Lo que aporta es que la represión a 37 °C **no alcanza por
igual a todos los blancos publicados**: de los candidatos de membrana externa
del mismo autor, cuatro no cambian y dos suben.

### ialB corresponde a RS01570, y no cambia

P35641.1 empareja con **`BARBAKC583_RS01570`** al **100 % de identidad y 100 %
de cobertura**, sin segundo acierto por encima de e < 1 × 10⁻⁵. El genoma tiene
seis proteínas anotadas `invasion associated locus B family protein` y la
anotación **no asigna el símbolo `ialB` a ninguna**: la cadena `ialB` aparece
cero veces en el GFF.

| | log2FC 37/25 | padj | log2FC pH 6/8 | padj |
|---|---|---|---|---|
| RS01570 (ialB) | +0,391 | 0,214 | −0,259 | 0,453 |

**El transcrito de `ialB` no varía ni con la temperatura ni con el pH** en las
12 muestras del modelo.

El resultado no se debe a falta de potencia. `RS01570` es el transcrito
**número 17 de 1 189** por abundancia (TPM mediano 7 657,6; percentil 98,6):

| Condición | Pl25 | Pl30 | Pl37 | pH06 | pH07 | pH08 |
|---|---|---|---|---|---|---|
| TPM | 5 987,9 | 8 040,7 | 6 379,9 | 7 274,6 | 9 308,2 | 9 466,8 |

**Lo que este resultado no autoriza a decir.** No se ha leído la sección de
métodos de Coleman y Minnick 2003: se desconocen sus temperaturas, medio,
tiempos de muestreo y si midieron transcrito o proteína. Por tanto puede
escribirse que en este conjunto de datos la abundancia del transcrito no varía,
pero **no** que no se replicó ese trabajo. Esa lectura queda pendiente.

Los seis parálogos se parecen entre sí solo 20 a 27 % por BLASTp, salvo
`RS01500` y `RS01510`, que dan 94,9 % entre ellos y son una duplicación en
tándem separada por 5,5 kb en la misma hebra.

### Los cinco candidatos de Gul et al. 2020

| Consulta | Locus | Identidad | Producto | log2FC 37/25 | padj |
|---|---|---|---|---|---|
| A1UTU1 | RS05550 | 99,6 % | FliP, poro del sistema de secreción tipo III flagelar | −1,265 | 1,70 × 10⁻³ |
| A1US20 | RS02205 | 100 % | Permeasa ABC de la familia FecCD | −1,537 | 1,86 × 10⁻⁸ |
| A1USE3 | RS02820 | 100 % | Fosfatidato citidililtransferasa | +0,550 | 0,0257 |
| A1UU82 | RS06260 | 100 % | Permeasa ABC del hemo | +0,213 | 0,24 |
| A1URZ7 | RS02080 | 100 % | Hemo liasa, subunidad familia CcmF/NrfE | +0,174 | 0,378 |

Dos reprimidos, uno levemente inducido, dos sin cambio. Ninguno de los cinco
cambia con el pH: los cinco dan padj mayor a 0,25. La identidad de 99,6 % de
A1UTU1 corresponde a una reanotación del codón de inicio (247 aa en UniProt,
250 en RefSeq), no a dos proteínas distintas.

`RS02205` está en el agrupamiento de captación de hemo y hierro, que se reprime
entero (RS02200 −1,488, RS02205 −1,537, RS02210 −1,696, RS02215 −1,423; todos
con padj menor que 1,2 × 10⁻⁶). Al ser represión coordinada de un operón, la
observación no depende de haber acertado el locus exacto que eligió Gul.

### Asignación ambigua

`A1USN1` queda marcada `ambigua`: acierta al 100 % de identidad y cobertura,
pero el segundo acierto llega al 94,6 %. Es la pareja TolC, y **la tesis misma
la cita como "A1USK0/A1USN1"** porque tampoco las distingue. No debe elegirse
una de las dos: en el manuscrito se citan como un par.

### Cifras retiradas de este cruce

Cinco cifras obtenidas a mano, antes de existir el script, resultaron
equivocadas y **no deben usarse**:

| Cifra retirada | Cifra correcta |
|---|---|
| ialB en RS01500: −1,239 (padj 0,0032) en temperatura y +0,75 (padj 0,0317) en pH, presentado como concordante con Coleman y Minnick 2003 | ialB es RS01570: +0,391 (padj 0,214) y −0,259 (padj 0,453), **sin cambio en ninguno de los dos ejes** |
| A1UU82 en CcmA RS06250, inducido +1,063, padj 0,0112 | RS06260, +0,213, padj 0,24, sin cambio |
| A1URZ7 en el agrupamiento del hemo reprimido, −1,42 a −1,70 | RS02080, +0,174, padj 0,378, sin cambio |
| A1US20 en RS02200, RS02210 o RS02215 | RS02205, −1,537, padj 1,86 × 10⁻⁸ |
| Fisher de especificidad: p = 1,0 × 10⁻³ (2 de 10), p = 3,365 × 10⁻⁴ (3 de 15) y p = 7,088 × 10⁻⁵ (5 de 26) | **p = 2,761 × 10⁻⁵**, sobre 15 de 48 no flagelares deduplicados por locus | Cruce cerrado en siete estudios el 12-11-2026 |

`RS06250` (CcmA) sí está inducido, con +1,063 y padj 0,0112, pero **no es un
candidato de Gul**: era una atribución errónea.

La cifra retirada de `ialB` se había elegido por ser, de los seis parálogos, el
que reproducía el resultado esperado. **Cualquier elección del parálogo hecha
por nombre de producto, o por cuál de los candidatos da el resultado esperado,
es circular y no es admisible.** También se había dado 170 aa como longitud de
IalB; esa cifra corresponde a IalA (P35640.1), no a IalB.

---

## 17. Condiciones excluidas: no concluyente (12 de noviembre de 2026)

`scripts/14_exploratorio_sangre.py`. **Exploratorio, sin valores p, fuera de
`verificar.sh`.** Mide el percentil mediano del conjunto flagelar dentro de
cada muestra, a partir del TPM de `quant.sf`. El rango es insensible a la
profundidad y a la composición de la biblioteca, que es lo que difiere entre
matrices.

| Condición | Centro | Percentil mediano del conjunto flagelar |
|---|---|---|
| Pl25 | Yale | 49,1 · 62,6 |
| Pl30 | Yale | 18,0 · 19,7 |
| Pl37 | Yale | 25,0 · 16,2 |
| pH06 | Yale | 67,5 · 69,6 |
| pH07 | Yale | 50,3 · 53,2 |
| pH08 | Yale | 66,3 · 62,4 |
| PlBG | GENEWIZ | 64,0 · 75,1 · 64,3 |
| HB37 | GENEWIZ | 61,7 · 56,3 · 67,5 |
| HBBG | GENEWIZ | 74,6 · 59,5 · 60,5 |
| HUVE | GENEWIZ | 18,2 · 23,7 |

**El método es válido:** en placa reproduce el resultado conocido, alto a 25 °C
y bajo a 30 y 37 °C. No hay artefacto de ceros: entre 95 y 99 % de los
transcritos tienen señal en todas las muestras.

**Pero el resultado no es interpretable.** En sangre humana a 37 °C el conjunto
flagelar queda alto, lo contrario de lo que se esperaría. Eso no se puede leer
como biología: **PlBG es placa a 37 °C igual que Pl37** y marca 68 frente a 21.
Esas dos se diferencian en la atmósfera sanguínea y en el centro de preparación
de biblioteca. El centro está perfectamente confundido con la condición en las
11 muestras excluidas (sección 1).

**Qué hacer con esto.** Va en limitaciones, para explicar por qué no se usaron
las condiciones excluidas, y **no como apoyo del resultado principal**.

### Observación aparte, esta sí dentro del modelo

A 30 °C y pH 7,4, mismo centro y mismo tiempo de traslado, el conjunto flagelar
marca ~19 en placa (Pl30) y ~52 en líquido (pH07). Solo cambia sólido frente a
líquido. **Pendiente**: contraste formal con DESeq2, posible porque las dos
condiciones están entre las 12 del modelo.
