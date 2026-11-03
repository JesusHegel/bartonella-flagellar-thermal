# Candidatos publicados por otros estudios

Secuencias de las proteinas que otros trabajos propusieron como blanco vacunal
o de diagnostico en *Bartonella bacilliformis*. El script
`scripts/13_candidatos_publicados.sh` las busca por blastp en el proteoma de
KC583 y les pega al lado los valores de expresion de este reanalisis.

## Por que las secuencias se guardan aqui y no se descargan al correr

Son pocas y pesan muy poco, asi que se versionan en el repositorio. Asi la
replica no depende de que UniProt o el NCBI sigan en linea ni de que no cambien
la entrada. El origen de cada una queda anotado en `fuentes.tsv` y abajo.

## Archivos

`fuentes.tsv`, una fila por consulta, separado por tabuladores:

| Columna | Contenido |
|---|---|
| consulta | Identificador corto. Debe ser igual a la cabecera en `secuencias.faa` |
| acceso | Numero de acceso de donde se bajo la secuencia |
| estudio | Quien la propuso |
| nombre_publicado | Nombre que le da el estudio, o `NA` si no le da ninguno |
| lo_que_reporta | En una linea, que afirma ese estudio sobre esta proteina |

`secuencias.faa`, las secuencias en FASTA. La cabecera de cada una tiene que
ser exactamente el valor de `consulta`, sin nada mas:

```
>A1UTU1
MKR...
```

El script se detiene si las cabeceras y la columna `consulta` no coinciden.

## Como obtener las secuencias

### Candidatos de UniProt (Gul et al. 2020)

Son proteinas de esta misma cepa con otro identificador, asi que se espera
identidad cercana al 100% contra el proteoma de referencia.

**Los cinco accesos ya no estan en UniProtKB.** Pedirlos a
`rest.uniprot.org/uniprotkb/<acceso>.fasta` devuelve respuesta 200 con el
cuerpo vacio (`content-length: 0`), que es facil de confundir con un error de
red. Fueron retirados de la base activa y hoy viven solo en UniParc, el
archivo historico. Se bajan de ahi, buscando por el acceso retirado:

```bash
cd datos/candidatos_publicados
: > secuencias.faa
: > procedencia_uniparc.txt
for A in A1UTU1 A1UU82 A1URZ7 A1US20 A1USE3; do
  curl -sL --max-time 30 \
    "https://rest.uniprot.org/uniparc/search?query=$A&format=fasta&size=1" \
    -o /tmp/up_$A.faa
  N=$(grep -c '^>' /tmp/up_$A.faa)
  if [ "$N" -ne 1 ]; then echo "FALLO $A (devolvio $N registros)"; continue; fi
  printf '%s\t%s\t%s aa\n' "$A" \
    "$(head -1 /tmp/up_$A.faa | sed 's/^>//')" \
    "$(grep -v '^>' /tmp/up_$A.faa | tr -d '\n' | wc -c)" >> procedencia_uniparc.txt
  awk -v a=$A 'NR==1{print ">"a; next}{print}' /tmp/up_$A.faa >> secuencias.faa
done
grep -c '^>' secuencias.faa      # deben salir 5
cat procedencia_uniparc.txt
```

El `awk` cambia la cabecera de UniParc por el identificador corto, que es lo
que el script espera. `procedencia_uniparc.txt` anota a que `UPI` resolvio cada
acceso y de cuantos aminoacidos es, y se versiona junto con las secuencias:
es el registro de donde salio cada una.

Lo obtenido el 3 de noviembre de 2026, version de UniProt 2026_03:

```
A1UTU1  UPI0000ED0A9F status=active     247 aa
A1UU82  UPI0000674951 status=active     259 aa
A1URZ7  UPI000067478A status=active     659 aa
A1US20  UPI000067479D status=active     366 aa
A1USE3  UPI0000ED09FE status=active     272 aa
```

Los 247 aa de A1UTU1 coinciden con la ficha de UniParc de FliP, que es la
unica de las cinco cuyo nombre da Gul en el articulo.

No sirve `rest.uniprot.org/unisave/<acceso>?format=fasta`: responde, pero
devuelve todas las versiones historicas de la entrada concatenadas, no una
secuencia.

### IalB (Mitchell y Minnick 1995)

Falta el numero de acceso. Buscalo en el NCBI y pega la secuencia de proteina:

<https://www.ncbi.nlm.nih.gov/protein/?term=ialB+Bartonella+bacilliformis>

Cuando lo tengas, anotalo en la columna `acceso` de `fuentes.tsv` (hoy dice
`PENDIENTE`) y agrega la secuencia:

```bash
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=ACCESO&rettype=fasta&retmode=text" \
  | awk 'NR==1{print ">ialB_publicado"; next}{print}' >> secuencias.faa
```

Esta es la consulta que importa mas. El genoma tiene **seis** proteinas
anotadas como `invasion associated locus B family protein` y la anotacion no le
asigna el simbolo `ialB` a ninguna. Sin esta secuencia no se puede decir cual
de las seis es el gen del que hablan Coleman y Minnick, y cualquier eleccion
hecha por nombre, o peor, hecha porque su expresion coincide con lo esperado,
es circular.

### Candidatos flagelares de otros estudios

Khan et al. 2019 (FlhA) y Rivera 2022 (FlgK, FliF, FliG, FlhB, FlhA, MotA). Si
el trabajo da numero de acceso, se agregan igual que los de arriba. Si solo da
el nombre de la proteina, el emparejamiento por nombre **si** es valido en este
caso particular, porque ninguno de los 31 genes flagelares comparte su nombre
de producto con otro gen del genoma. Eso se comprueba con la columna
`n_genes_mismo_producto` de la tabla de salida, que vale 1 para todos ellos.

## Salida

`resultados/candidatos/candidatos_publicados.tsv`, una fila por consulta.
Columnas de control que hay que leer antes de usar cualquier fila:

| Columna | Que mirar |
|---|---|
| asignacion | Solo `por secuencia` sirve para el manuscrito. `ambigua`, `acierto debil` y `sin acierto` se revisan a mano |
| ambiguo | `TRUE` si el segundo acierto casi empata con el primero, es decir si la secuencia tampoco distingue entre dos genes |
| n_genes_mismo_producto | Cuantos genes del genoma comparten ese nombre de producto. Si es mayor que 1, emparejar por nombre nunca habria bastado |

`resultados/candidatos/candidatos_blastp.tsv` trae todos los aciertos de cada
consulta, en orden, para revisar los casos marcados.
