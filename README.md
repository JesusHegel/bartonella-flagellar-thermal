# Regulación térmica del regulón flagelar en *Bartonella bacilliformis*

Reanálisis de los transcritos codificantes de proteína del BioProject
**PRJNA647605** (Wachter *et al.*, 2020), conjunto generado originalmente para
caracterizar ARN pequeños no codificantes.

Este repositorio contiene el código, los conjuntos de datos derivados y las
figuras del manuscrito asociado.

---

## Resultado principal

El contraste entre 25 °C y 37 °C modifica 216 de 1 176 transcritos. El conjunto
funcional más afectado es el aparato flagelar: 29 de los 31 genes flagelares
anotados presentan expresión reducida a temperatura de hospedero.

Dos métodos de enriquecimiento independientes coinciden en identificarlo como el
conjunto más desplazado: test hipergeométrico (q = 1,3 × 10⁻¹⁶) y análisis sobre
el ranking completo sin umbrales (NES = −2,46; q = 1,8 × 10⁻⁹).

---

## Estructura
---

## Archivo clave: `data/flag_gff.tsv`

Contiene los 31 genes flagelares empleados en todo el análisis.

El conjunto se define en `scripts/01_definicion_conjunto_flagelar.sh` mediante
búsqueda por patrón sobre el campo `product` del archivo de anotación GFF del
genoma de referencia. **Ese script no consulta ningún dato de expresión.**

Se publica de forma explícita porque constituye la evidencia de que el conjunto
se estableció antes del análisis diferencial y no a partir de sus resultados.
Los 31 genes son la totalidad de los que tienen anotación flagelar en el genoma,
no una selección.

---

## Datos de origen

| Recurso | Identificador |
|---|---|
| Datos de RNA-seq | BioProject PRJNA647605 (SRA) |
| Genoma de referencia | GCF_000015445.1 (*B. bacilliformis* KC583) |
| Cepa para ortología | GCF_001624625.1 (USM-LMMB07) |
| Genomas del género | GCF_000253015.1, GCF_039555305.1, GCF_001281405.1, GCF_000196435.1, GCF_009936175.1, GCF_019930925.1 |
| Organismos para quimiotaxis | GCF_037023865.1, GCF_000022005.1, GCF_000092025.1, GCF_003324715.1 |

Los datos crudos no se incluyen en este repositorio: están disponibles en SRA y
se recuperan mediante los identificadores anteriores.

---

## Programas utilizados

salmon v2.4.1 · DESeq2 v1.50.2 · ashr · fgsea · BLAST+ · skani · CheckM2 ·
ViennaRNA 2.7.0 · NCBI Datasets · R 4.5.3

---

## Orden de ejecución

1. `01_definicion_conjunto_flagelar.sh` — extrae el conjunto flagelar de la anotación
2. `02_expresion_diferencial.R` — modelo global y modelo restringido a un lote técnico
3. `03_enriquecimiento.R` — enriquecimiento por dos métodos independientes
4. `04_panel_genero.sh` — distribución del regulón en siete especies del género
5. `05_quimiotaxis.sh` — búsqueda por homología de componentes del sistema Che
6. `06_figuras.R` — generación de las figuras

---

## Alcance y limitaciones

El estudio es enteramente computacional. No se realizaron experimentos de
laboratorio. Los resultados corresponden a niveles de transcrito, no de proteína.
El diseño original dispone de dos réplicas biológicas para siete de las diez
condiciones.

Las limitaciones completas se detallan en el manuscrito asociado.

---

## Licencia

Código bajo licencia MIT. Los datos derivados de anotaciones públicas de NCBI
mantienen las condiciones de uso de su fuente original.

---

## Cita

Repositorio asociado a un manuscrito en preparación. La referencia se
incorporará tras su publicación.

Los datos originales corresponden a:
Wachter S, Hicks LD, Raghavan R, Minnick MF (2020). *Novel small RNAs expressed
by Bartonella bacilliformis under multiple conditions reveal potential
mechanisms for persistence in the sand fly vector and human host.*
PLoS Negl Trop Dis 14(11):e0008671.
