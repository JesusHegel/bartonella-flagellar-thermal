# Guía de réplica

Esta guía permite repetir el análisis completo en otra computadora y comprobar
que los resultados son los mismos. Está pensada para Windows con WSL (Linux
dentro de Windows), pero los pasos 3 a 6 son iguales en cualquier Linux o Mac.

Hay dos niveles de réplica:

| Nivel | Parte de | Resultado esperado | Tiempo aproximado |
|---|---|---|---|
| **1** (obligatorio) | Las cuantificaciones de salmon incluidas en `cuantificacion/` | Tablas **idénticas byte a byte**, salvo las numéricas de DESeq2 y derivadas, que pueden variar en el último decimal (se verifican con tolerancia) | 10–20 min, más la instalación |
| **2** (si el equipo lo permite) | Las lecturas crudas del SRA | Mismas conclusiones, con diferencias mínimas en decimales | 1–2 h, necesita internet y 25 GB libres |

El nivel 2 no puede dar archivos idénticos: salmon con varios hilos reparte de
forma ligeramente distinta unas pocas lecturas ambiguas en cada corrida (en una
prueba con la muestra SRR12653239, 214 de 1 040 336 lecturas; ver
`CIFRAS_CONFIRMADAS.md`, sección 12b).

---

## 1. Instalar Linux dentro de Windows (WSL)

Solo la primera vez. Abre **PowerShell como administrador** y ejecuta:

```powershell
wsl --install -d Debian
```

Reinicia la computadora cuando lo pida. Al volver, se abre una ventana de
Debian que pide crear un usuario y una contraseña de Linux. Desde ese momento,
todo se hace en esa ventana ("terminal de Debian").

## 2. Instalar herramientas básicas y conda

En la terminal de Debian:

```bash
sudo apt update && sudo apt install -y git wget unzip
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b
~/miniforge3/bin/conda init bash
exec bash
```

Al final, el inicio de la línea debe mostrar `(base)`.

## 3. Descargar el repositorio

```bash
git clone https://github.com/JesusHegel/bartonella-flagellar-thermal.git ~/bartonella
cd ~/bartonella
```

## 4. Crear el entorno con las versiones exactas

```bash
conda env create -f entorno/bartonella.yml
conda activate bartonella
```

Tarda entre 15 y 40 minutos. Instala R, DESeq2, salmon, BLAST y el resto de
programas en las mismas versiones que se usaron en el análisis. Usar otras
versiones puede cambiar algunos resultados (por ejemplo, el umbral formal de
DESeq2 da números distintos entre versiones).

## 5. Nivel 1: correr todo y verificar

```bash
cd ~/bartonella
conda activate bartonella
bash scripts/ejecutar_todo.sh
```

El script borra `resultados/`, `figuras/` e `intermedios/`, los regenera desde
cero y al final compara cada tabla con la versión publicada. La última línea
debe decir una de estas dos cosas:

```
RESULTADO: todo identico.
RESULTADO: correcto. N tablas identicas y M equivalentes dentro de la tolerancia.
```

La segunda es normal en otra computadora: los contrastes de DESeq2, fgsea, la
PCA y las tablas que se derivan de ellos pueden variar en el último decimal
según el procesador. Se aceptan solo si
ningún gen cambia de significativo a no significativo (ver
`CIFRAS_CONFIRMADAS.md`, sección 12c).

Los scripts 06 y 07 descargan proteomas del NCBI, así que necesitan internet.
Las figuras se revisan a la vista y no se comparan por md5: en otra computadora
cambian las fuentes y el dibujo varía un poco aunque los datos sean los mismos.
La figura 5 puede variar incluso entre dos corridas en la misma computadora:
ggrepel coloca algunas etiquetas en otra posición. Los nodos, las aristas y la
disposición de la red son siempre los mismos.

**Qué enviar a Hegel:** la salida completa de la pantalla y la carpeta
`registros/`.

Si aparece algún `DISTINTO` o `FALTA`, no sigas: envía la salida tal cual.
Para volver a los archivos originales: `git checkout -- resultados figuras registros`.

## 6. Nivel 2: desde las lecturas crudas del SRA

Necesita internet y unos 25 GB libres en el disco de Windows (cada corrida se
descarga, se cuantifica y se borra antes de la siguiente).

```bash
cd ~/bartonella
conda activate bartonella
bash scripts/replica_nivel_2.sh
```

Hace tres cosas, en orden:

1. Reconstruye el índice de salmon, comprueba que sus seis huellas coinciden con
   el original, y descarga y cuantifica las 23 corridas en
   `cuantificacion_desde_sra/` (script 00). Si se corta, vuelve a correr el mismo
   comando: retoma donde quedó.
2. Corre el análisis completo sobre esas cuantificaciones, en carpetas aparte que
   terminan en `_desde_sra`. Nunca toca `cuantificacion/` ni `resultados/`.
3. Compara las cifras principales con las publicadas (`comparar_nivel_2.R`).

No corras `ejecutar_todo.sh` mientras trabaja: borra `intermedios/`, donde está
el índice de salmon.

Desde el SRA **no se espera identidad byte a byte**: salmon con varios hilos no
es determinista. La tabla A muestra las lecturas de cada corrida y la tabla B
las cifras principales publicadas junto a las obtenidas desde el SRA. Lo que
importa es que no cambie ninguna conclusión.

**Qué enviar a Hegel:** `registros_desde_sra/comparar_nivel_2.log`.

---

## Si algo falla

- `ERROR: ejecuta los scripts desde la carpeta raiz del repositorio`: haz
  `cd ~/bartonella` y vuelve a intentarlo.
- `FALTA Rscript` (u otro programa): falta `conda activate bartonella`.
- Error de un script concreto: su salida completa está en
  `registros/<nombre del script>.log`.
