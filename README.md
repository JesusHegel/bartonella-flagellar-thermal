# Regulación térmica del regulón flagelar en *Bartonella bacilliformis*

Reanálisis de los transcritos codificantes de proteína del BioProject
**PRJNA647605** (Wachter *et al.*, 2020), generado originalmente para
caracterizar ARN pequeños no codificantes.

> **Todas las cifras válidas están en [`CIFRAS_CONFIRMADAS.md`](CIFRAS_CONFIRMADAS.md).**
> Es la única fuente de verdad del proyecto. Cada cifra tiene un script que la produce.
>
> **Para replicar el análisis en otra computadora, sigue [`GUIA_DE_REPLICA.md`](GUIA_DE_REPLICA.md).**

---

## Resultado principal

El contraste 37 °C frente a 25 °C modifica **488 de 1 186 transcritos**
(207 con |log₂ cambio| > 1). El conjunto funcional más afectado es el aparato
flagelar: **los 31 genes flagelares anotados resultan significativos**, 30
reprimidos y uno inducido (FliO).

| Método | Estadístico | q |
|---|---|---|
| Hipergeométrico | 30 reprimidos frente a 7,6 esperados | **2,1 × 10⁻¹⁷** |
| fgsea, sobre el ranking completo sin umbrales | NES = −2,44 | **2,4 × 10⁻⁹** |

**El cambio se completa a 30 °C.** El contraste 37 °C frente a 30 °C no es
significativo en 30 de los 31 genes. Las chaperonas principales, en cambio,
siguen aumentando entre 30 y 37 °C. Los dos programas tienen perfiles térmicos
distintos.

**Especificidad.** Frente a un contraste de pH de magnitud comparable en el
genoma (375 transcritos modificados, frente a 488 por temperatura), solo
**6 de los 31** genes flagelares responden.

**Distribución en el género.** Los 31 genes están completos en
*B. bacilliformis*, *B. clarridgeiae*, *B. schoenbuchensis* y *B. ancashensis*,
y se reducen a **2 de 31** en *B. tribocorum*, *B. quintana* y *B. henselae*.
Los 31 tienen ortólogo recíproco en la cepa peruana USM-LMMB07.

El conjunto flagelar se define en `scripts/01_anotacion.sh` por palabras clave
sobre la anotación del genoma. **Ese script no consulta ningún dato de
expresión**: el conjunto se fijó antes de mirar los resultados, y los 31 son la
totalidad de los genes con anotación flagelar.

---

## Estructura

    CIFRAS_CONFIRMADAS.md   Todas las cifras del proyecto, con su procedencia
    GUIA_DE_REPLICA.md      Paso a paso para replicar el análisis en otra computadora
    entorno/                Programas y versiones exactas (entorno de conda)
    datos/                  Tabla de muestras, genoma de referencia KC583 del NCBI y copia de STRING
    cuantificacion/         Cuantificación de salmon de las 23 corridas (punto de partida)
    scripts/                Scripts numerados en el orden en que se ejecutan
    resultados/             Tablas que producen los scripts
    figuras/                Figuras (PNG para pantalla, PDF para imprenta)
    registros/              Salida de pantalla de cada script en la última corrida

### Scripts

| Script | Produce |
|---|---|
| `00_descargar_y_cuantificar.sh` | Solo para la réplica completa: descarga las 23 corridas del SRA y las cuantifica con salmon |
| `01_anotacion.sh` | Anotación de las CDS y conjunto flagelar a priori (`resultados/anotacion/`) |
| `02_expresion_diferencial.R` | Modelos de DESeq2 y los cinco contrastes (`resultados/expresion_diferencial/`) |
| `03_enriquecimiento.R` | Enriquecimiento hipergeométrico y fgsea (`resultados/enriquecimiento/`) |
| `04_robustez.R` | Umbrales, modelos, curva térmica, especificidad, particiones (`resultados/robustez/`) |
| `05_pca_y_tasa_mapeo.R` | PCA frente a tasa de asignación por muestra (`resultados/muestras/`) |
| `06_panel_genero_y_ortologia.sh` | BLASTp en siete especies y ortología con USM-LMMB07 (`resultados/genomica_comparada/`) |
| `07_quimiotaxis_y_utr.sh` | Homólogos del sistema Che y estructura del 5'UTR de la flagelina |
| `08_tabla_de_genes.R` | Qué gen es cada transcrito, TPM por condición y los 20 más inducidos y reprimidos (`resultados/genes/`) |
| `09_tabla_de_muestras.R` | Lecturas de cada corrida y motivo de inclusión o exclusión (`resultados/muestras/`) |
| `10_descargar_string.sh` | Descarga STRING para KC583 solo si no está guardado en `datos/string/` |
| `11_analisis_string.R` | Enriquecimiento con la anotación de STRING y módulo flagelar en la red (`resultados/string/`) |
| `12_figuras.R` | Las seis figuras del manuscrito (`figuras/`): flujo, muestras, volcanes, flagelares y chaperonas, red STRING y panel del género |
| `ejecutar_todo.sh` | Corre del 01 al 12 y verifica el resultado |
| `verificar.sh` | Compara cada tabla regenerada con la versión publicada |

### Uso rápido

```bash
conda env create -f entorno/bartonella.yml     # solo la primera vez
conda activate bartonella
bash scripts/ejecutar_todo.sh
```

La corrida completa tarda unos minutos y termina comparando cada tabla con la
versión publicada. Casi todas deben ser idénticas byte a byte. Los contrastes de
DESeq2, fgsea y la PCA pueden variar en el último decimal según el procesador;
en ese caso se comparan con tolerancia y solo se aceptan si no cambia ninguna
decisión de significancia (ver `scripts/comparar_con_tolerancia.R`).

---

## Datos de origen

| Recurso | Identificador |
|---|---|
| RNA-seq | BioProject PRJNA647605 (SRA), 23 corridas pareadas |
| Genoma de referencia | GCF_000015445.1, *B. bacilliformis* KC583 (incluido en `datos/referencia/`) |
| Cepa peruana | GCF_001624625.1, USM-LMMB07 |
| Panel de género | GCF_000253015.1 · GCF_039555305.1 · GCF_001281405.1 · GCF_000196435.1 · GCF_009936175.1 · GCF_019930925.1 |
| Quimiotaxis | GCF_037023865.1 · GCF_000022005.1 · GCF_000092025.1 · GCF_003324715.1 |

Las lecturas crudas no se incluyen: se recuperan del SRA con el script 00. Las
23 corridas provienen de un solo centro de secuenciación (University of Montana,
HiSeq 2500). El genoma de referencia sí se incluye, porque el NCBI puede
reanotar un genoma con el tiempo y eso cambiaría los resultados sin aviso.

---

## Alcance y limitaciones

El estudio es **enteramente computacional**. Los resultados corresponden a
niveles de **transcrito, no de proteína**, y no demuestran que el flagelo esté
ausente ni que deje de funcionar.

- El diseño original tiene dos réplicas biológicas en siete de las diez condiciones.
- PC1 explica el 53 % de la varianza y correlaciona 0,82 con la tasa de
  asignación de lecturas a genes bacterianos. Por eso los contrastes se
  restringen a las 12 muestras de placa y caldo.
- Los 31 genes flagelares están en operones y no son independientes entre sí,
  de modo que el valor q nominal es más extremo que el real.
- Los conteos de DESeq2 no están normalizados por longitud: permiten comparar
  un mismo gen entre condiciones, no genes distintos entre sí.

---

## Licencia y cita

Código bajo licencia MIT. Los datos derivados de anotaciones del NCBI mantienen
las condiciones de su fuente original.

Repositorio asociado a un manuscrito en preparación.

Datos originales:
Wachter S, Hicks LD, Raghavan R, Minnick MF (2020). *Novel small RNAs expressed
by Bartonella bacilliformis under multiple conditions reveal potential mechanisms
for persistence in the sand fly vector and human host.*
PLoS Negl Trop Dis 14(11):e0008671.
