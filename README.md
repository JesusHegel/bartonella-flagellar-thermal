# Regulación térmica del regulón flagelar en *Bartonella bacilliformis*

Reanálisis de los transcritos codificantes de proteína del BioProject
**PRJNA647605** (Wachter *et al.*, 2020), generado originalmente para
caracterizar ARN pequeños no codificantes.

> **Todas las cifras válidas están en [`CIFRAS_CONFIRMADAS.md`](CIFRAS_CONFIRMADAS.md).**
> Ese archivo es la única fuente de verdad del proyecto. El análisis completo se
> regeneró desde cero el 22 de septiembre de 2026 en carpeta limpia, y cada cifra
> tiene un script que la produce.

---

## Resultado principal

El contraste 37 °C frente a 25 °C modifica **488 de 1 186 transcritos**
(207 con |log₂ cambio| > 1). El conjunto funcional más afectado es el aparato
flagelar: **los 31 genes flagelares anotados resultan significativos** — 30
reprimidos y uno inducido (FliO).

| Método | Estadístico | q |
|---|---|---|
| Hipergeométrico | 30 reprimidos frente a 7,6 esperados | **2,1 × 10⁻¹⁷** |
| fgsea, sobre el ranking completo sin umbrales | NES = −2,44 | **2,4 × 10⁻⁹** |

### El cambio se completa a 30 °C

El contraste 37 °C frente a 30 °C **no es significativo en 30 de los 31 genes**.
La represión flagelar ocurre entre 25 y 30 °C y no avanza más allá. La respuesta
de chaperonas, en cambio, sigue aumentando hasta 37 °C. Son dos programas
térmicos con formas distintas, de modo que la represión flagelar no es un efecto
secundario del choque térmico (`figures/Fig5_disociacion.png`).

### Especificidad

Frente a un contraste de pH de magnitud comparable en el genoma (375 transcritos
modificados, frente a 488 por temperatura), solo **6 de los 31** genes flagelares
responden. El conjunto no reacciona a cualquier estímulo.

### Distribución en el género

Los 31 genes están completos en *B. bacilliformis*, *B. clarridgeiae*,
*B. schoenbuchensis* y *B. ancashensis*, y se reducen a **2 de 31** en
*B. tribocorum*, *B. quintana* y *B. henselae*. Los 31 tienen ortólogo recíproco
en la cepa peruana USM-LMMB07.

---

## El conjunto flagelar se define antes del análisis

`scripts/01_anotacion.sh` extrae los 31 genes por patrón sobre el campo
`product` del archivo de anotación del genoma. **Ese script no consulta ningún
dato de expresión.** Es la evidencia de que el conjunto se estableció a priori y
no a partir de los resultados. Los 31 son la totalidad de los que tienen
anotación flagelar, no una selección.

---

## Estructura

    CIFRAS_CONFIRMADAS.md   Todas las cifras del proyecto, con su procedencia
    scripts/                Ocho scripts que regeneran cada resultado, en orden
    results/                Las tablas que producen
    figures/                Las cinco figuras (PNG para presentacion, PDF para imprenta)
    data/                   Tablas derivadas de la anotacion del genoma
    logs/                   Salida completa de cada script en la corrida del 22-09-2026

### Orden de ejecución

| Script | Produce |
|---|---|
| `01_anotacion.sh` | Anotación, conjunto flagelar a priori, duplicados |
| `02_deseq.R` | Objeto DESeq2, PCA, cinco contrastes |
| `03_enriquecimiento.R` | Enriquecimiento por dos métodos |
| `04_robustez.R` | Conteos, sensibilidad, particiones, curva térmica, especificidad |
| `05_pca_mapeo.R` | PC1 frente a tasa de mapeo |
| `06_figuras.R` | Las cinco figuras |
| `07_panel_genero_ortologia.sh` | Panel de género y ortología con la cepa peruana |
| `08_quimiotaxis_utr.sh` | Homología del sistema Che y estructura del 5'UTR |

---

## Datos de origen

| Recurso | Identificador |
|---|---|
| RNA-seq | BioProject PRJNA647605 (SRA) |
| Genoma de referencia | GCF_000015445.1 — *B. bacilliformis* KC583 |
| Cepa peruana | GCF_001624625.1 — USM-LMMB07 |
| Panel de género | GCF_000253015.1 · GCF_039555305.1 · GCF_001281405.1 · GCF_000196435.1 · GCF_009936175.1 · GCF_019930925.1 |
| Quimiotaxis | GCF_037023865.1 · GCF_000022005.1 · GCF_000092025.1 · GCF_003324715.1 |

Los datos crudos no se incluyen: se recuperan de SRA con los identificadores
anteriores. Las 23 corridas provienen de un solo centro de secuenciación
(University of Montana, HiSeq 2500).

---

## Alcance y limitaciones

El estudio es **enteramente computacional**. No se realizaron experimentos.
Los resultados corresponden a niveles de **transcrito, no de proteína**, y no
demuestran que el flagelo esté ausente ni que deje de funcionar.

- El diseño original tiene dos réplicas biológicas en siete de las diez condiciones.
- PC1 explica el 53 % de la varianza y correlaciona 0,82 con la fracción de
  lecturas bacterianas de cada librería; por eso los contrastes se restringen a
  las 12 muestras de matriz homogénea.
- Los 31 genes flagelares están en operones y no son independientes entre sí, de
  modo que el valor q nominal es más extremo que el real.
- Los conteos de DESeq2 no están normalizados por longitud: permiten comparar un
  mismo gen entre condiciones, no genes distintos entre sí.

---

## Versión anterior

El análisis previo a septiembre de 2026 está en el commit `c2bf641`. Sus cifras
fueron superadas; la sección 11 de `CIFRAS_CONFIRMADAS.md` detalla qué cambió.

---

## Licencia y cita

Código bajo licencia MIT. Los datos derivados de anotaciones de NCBI mantienen
las condiciones de su fuente original.

Repositorio asociado a un manuscrito en preparación.

Datos originales:
Wachter S, Hicks LD, Raghavan R, Minnick MF (2020). *Novel small RNAs expressed
by Bartonella bacilliformis under multiple conditions reveal potential mechanisms
for persistence in the sand fly vector and human host.*
PLoS Negl Trop Dis 14(11):e0008671.
