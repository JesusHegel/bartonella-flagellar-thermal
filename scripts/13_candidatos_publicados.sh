#!/bin/bash
# 13. Candidatos publicados por otros estudios, cruzados con la expresion.
#
# Localiza en el proteoma de KC583 las proteinas que otros trabajos propusieron
# como blanco vacunal o de diagnostico, y les pega al lado los valores de
# expresion medidos en este reanalisis.
#
# El emparejamiento se hace por SECUENCIA (blastp), nunca por nombre de
# producto. Razon: la anotacion de RefSeq repite nombres. En este genoma hay
# seis proteinas llamadas "invasion associated locus B family protein" y seis
# "ABC transporter permease", entre 84 nombres repetidos en total. Elegir por
# nombre deja la asignacion al azar, y elegir el parálogo cuya expresion
# coincide con lo esperado es razonamiento circular: con seis candidatos
# siempre hay uno que da el resultado que uno queria.
#
# Por eso la salida trae dos columnas de control:
#   n_genes_mismo_producto  cuantos genes del genoma comparten ese nombre
#                           (si es > 1, el nombre nunca habria bastado)
#   ambiguo                 TRUE si el segundo acierto casi empata con el
#                           primero, es decir si la secuencia tampoco decide
#
# Entrada (se versiona en el repositorio):
#   datos/candidatos_publicados/secuencias.faa  secuencias a buscar
#   datos/candidatos_publicados/fuentes.tsv     de donde viene cada una
#   Ver datos/candidatos_publicados/LEEME.md para como obtenerlas.
#
# Salida:
#   resultados/candidatos/candidatos_publicados.tsv  un renglon por consulta
#   resultados/candidatos/candidatos_blastp.tsv      todos los aciertos
#
# Se ejecuta despues del script 08 (necesita la tabla de genes).
set -uo pipefail
source scripts/rutas.sh

ENT=datos/candidatos_publicados
SAL=$DIR_RES/candidatos
mkdir -p "$SAL"

# Si todavia no estan las secuencias, el script avisa y sale bien, para que
# ejecutar_todo.sh no se detenga por esto: ningun otro resultado depende de
# este script.
for f in "$ENT/secuencias.faa" "$ENT/fuentes.tsv"; do
  if [ ! -s "$f" ]; then
    echo "AVISO: falta $f (o esta vacio); no se genera la tabla de candidatos."
    echo "       Como obtener las secuencias: $ENT/LEEME.md"
    exit 0
  fi
done
[ -s "$DIR_RES/genes/tabla_de_genes.csv" ] || { echo "ERROR: falta la tabla de genes; corre antes el script 08."; exit 1; }

N_CONS=$(grep -c '^>' "$ENT/secuencias.faa")
echo "Consultas en secuencias.faa: $N_CONS"

# Toda secuencia de secuencias.faa tiene que estar declarada en fuentes.tsv:
# una secuencia sin procedencia no puede entrar al analisis.
grep '^>' "$ENT/secuencias.faa" | sed 's/^>//; s/[[:space:]].*//' | sort -u > "$DIR_INT/cand_faa.txt"
tail -n +2 "$ENT/fuentes.tsv" | cut -f1 | sed '/^$/d' | sort -u > "$DIR_INT/cand_tsv.txt"
HUERFANAS=$(comm -23 "$DIR_INT/cand_faa.txt" "$DIR_INT/cand_tsv.txt")
if [ -n "$HUERFANAS" ]; then
  echo "ERROR: hay secuencias en secuencias.faa que no figuran en fuentes.tsv:"
  echo "$HUERFANAS" | sed 's/^/    /'
  echo "  Agregalas a fuentes.tsv o quitalas del FASTA."
  exit 1
fi

# Lo contrario si se permite: una fila de fuentes.tsv sin secuencia queda en la
# tabla marcada como "sin secuencia". Asi un candidato que todavia no se pudo
# conseguir sigue documentado en los resultados en vez de desaparecer.
FALTANTES=$(comm -13 "$DIR_INT/cand_faa.txt" "$DIR_INT/cand_tsv.txt")
if [ -n "$FALTANTES" ]; then
  echo "AVISO: estos candidatos de fuentes.tsv todavia no tienen secuencia:"
  echo "$FALTANTES" | sed 's/^/    /'
  echo "       Quedan en la tabla como 'sin secuencia'. Ver $ENT/LEEME.md"
else
  echo "Identificadores: todos los candidatos de fuentes.tsv tienen secuencia"
fi

# ---- 1. Base de datos del proteoma de referencia ----
DB=$DIR_INT/db_kc583
mkdir -p "$(dirname "$DB")"
makeblastdb -in "$REF/protein.faa" -dbtype prot -out "$DB" > "$DIR_INT/makeblastdb_kc583.log" 2>&1 \
  || { echo "ERROR: fallo makeblastdb"; exit 1; }
echo "Proteoma de referencia: $(grep -c '^>' "$REF/protein.faa") proteinas"

# ---- 2. BLASTp ----
# -num_threads 1 para que el orden de los aciertos empatados sea siempre el
# mismo: la tabla de salida se compara byte a byte con verificar.sh.
# Se guardan hasta 10 aciertos por consulta para poder medir la ambiguedad.
blastp -query "$ENT/secuencias.faa" -db "$DB" \
  -evalue 1e-5 -max_hsps 1 -max_target_seqs 10 -num_threads 1 \
  -outfmt '6 qseqid sseqid pident qcovs qlen slen bitscore' \
  > "$DIR_INT/blastp_candidatos.tsv" 2> "$DIR_INT/blastp_candidatos.log" \
  || { echo "ERROR: fallo blastp"; exit 1; }
echo "Aciertos de blastp: $(wc -l < "$DIR_INT/blastp_candidatos.tsv")"

# ---- 3. Cruce con la anotacion y la expresion ----
python3 - "$ENT/fuentes.tsv" "$DIR_INT/blastp_candidatos.tsv" \
         "$DIR_RES/anotacion/locus_a_proteina.tsv" "$DIR_RES/genes/tabla_de_genes.csv" \
         "$SAL/candidatos_publicados.tsv" "$SAL/candidatos_blastp.tsv" \
         "$DIR_INT/cand_faa.txt" <<'PY'
import csv, sys, os, collections
from math import comb

f_fuentes, f_blast, f_lp, f_genes, f_sal, f_sal_todo, f_faa = sys.argv[1:8]
con_secuencia = {l.strip() for l in open(f_faa) if l.strip()}

# Umbrales de asignacion. Los candidatos publicados suelen ser proteinas de
# esta misma cepa con otro identificador, asi que se espera identidad cercana
# al 100%. Se baja a 90/80 para tolerar reanotaciones del codon de inicio
# (por ejemplo FliP figura con 247 aa en UniProt y 250 en RefSeq).
MIN_IDENT, MIN_COBERTURA = 90.0, 80.0
# Si el segundo acierto llega al 95% del puntaje del primero, la secuencia no
# distingue entre los dos y la asignacion queda marcada como ambigua.
MARGEN = 0.95
# blastp mide la cobertura DE LA CONSULTA, asi que una consulta recortada da
# 100% de cobertura aunque le falte media proteina. Por eso se compara ademas
# el largo de las dos secuencias: deben parecerse.
LARGO_MIN, LARGO_MAX = 0.80, 1.25

fuentes = list(csv.DictReader(open(f_fuentes), delimiter='\t'))
orden = [r['consulta'] for r in fuentes]
fuente = {r['consulta']: r for r in fuentes}

# protein_id -> locus
locus_de = {}
for ln in open(f_lp):
    p = ln.rstrip('\n').split('\t')
    if len(p) >= 2:
        locus_de[p[1]] = p[0]

genes = {r['locus']: r for r in csv.DictReader(open(f_genes))}
n_producto = collections.Counter(r['producto'] for r in genes.values())

# aciertos por consulta, ya ordenados por blastp (mejor primero)
aciertos = collections.defaultdict(list)
for ln in open(f_blast):
    q, s, pid, cov, ql, sl, bits = ln.rstrip('\n').split('\t')
    aciertos[q].append({'proteina': s, 'identidad': float(pid),
                        'cobertura': float(cov), 'largo_q': int(ql),
                        'largo_s': int(sl), 'puntaje': float(bits)})

COLS = ['consulta', 'acceso', 'estudio', 'nombre_publicado',
        'proteina_kc583', 'identidad_pct', 'cobertura_pct',
        'largo_publicado_aa', 'largo_kc583_aa',
        'segundo_acierto', 'segundo_identidad_pct',
        'locus', 'producto', 'n_genes_mismo_producto', 'conjunto_flagelar',
        'log2FC_37vs25', 'padj_37vs25', 'log2FC_pH6vspH8', 'padj_pH6vspH8',
        'asignacion', 'ambiguo']

vacio = {c: 'NA' for c in COLS}
filas, filas_todo = [], []

for q in orden:
    r = dict(vacio)
    r['consulta'] = q
    for c in ('acceso', 'estudio', 'nombre_publicado'):
        r[c] = fuente[q].get(c, 'NA') or 'NA'
    hs = aciertos.get(q, [])

    # tabla larga: todos los aciertos, para poder revisar a mano
    for i, h in enumerate(hs, 1):
        L = locus_de.get(h['proteina'], 'NA')
        filas_todo.append({'consulta': q, 'puesto': i, 'proteina_kc583': h['proteina'],
                           'locus': L, 'identidad_pct': '%.1f' % h['identidad'],
                           'cobertura_pct': '%.0f' % h['cobertura'],
                           'largo_publicado_aa': h['largo_q'],
                           'largo_kc583_aa': h['largo_s'],
                           'puntaje': '%.0f' % h['puntaje'],
                           'producto': genes.get(L, {}).get('producto', 'NA')})

    if q not in con_secuencia:
        # Declarado en fuentes.tsv pero sin secuencia todavia. Se deja en la
        # tabla para que el pendiente quede a la vista en los resultados.
        r['asignacion'] = 'sin secuencia'
        r['ambiguo'] = 'NA'
        filas.append(r)
        continue

    if not hs:
        r['asignacion'] = 'sin acierto'
        r['ambiguo'] = 'NA'
        filas.append(r)
        continue

    m = hs[0]
    r['proteina_kc583'] = m['proteina']
    r['identidad_pct'] = '%.1f' % m['identidad']
    r['cobertura_pct'] = '%.0f' % m['cobertura']
    r['largo_publicado_aa'] = str(m['largo_q'])
    r['largo_kc583_aa'] = str(m['largo_s'])

    # El segundo acierto se reporta siempre, para poder juzgar a ojo cuanto
    # margen hay. No alcanza con la bandera: entre RS01500 y RS01510, que son
    # una duplicacion en tandem, el segundo llega al 94.9% de identidad.
    if len(hs) > 1:
        r['segundo_acierto'] = hs[1]['proteina']
        r['segundo_identidad_pct'] = '%.1f' % hs[1]['identidad']
    else:
        r['segundo_acierto'] = 'ninguno'
        r['segundo_identidad_pct'] = 'NA'

    ambiguo = len(hs) > 1 and hs[1]['puntaje'] >= MARGEN * m['puntaje']
    r['ambiguo'] = 'TRUE' if ambiguo else 'FALSE'
    razon = m['largo_q'] / m['largo_s'] if m['largo_s'] else 0.0

    # qcovs mide cobertura de la CONSULTA. Si la consulta es mas larga que el
    # blanco, por ejemplo porque la accesion citada por el articulo fue
    # extendida en la base de datos despues de publicarse, la cobertura de
    # consulta baja aunque el blanco este cubierto entero. Ese es el caso
    # seguro, no el peligroso: lo que el control debe cazar es una consulta
    # RECORTADA, que da 100 % de cobertura faltandole media proteina.
    # Se calcula tambien la cobertura del BLANCO, con los campos que ya hay.
    cob_s = min(100.0, m['cobertura'] / 100.0 * m['largo_q'] / m['largo_s'] * 100) \
        if m['largo_s'] else 0.0
    cob_efectiva = max(m['cobertura'], cob_s)
    blanco_entero = cob_s >= 95.0 and m['identidad'] >= MIN_IDENT
    rescatados = globals().setdefault('rescatados', [])

    if m['identidad'] < MIN_IDENT or cob_efectiva < MIN_COBERTURA:
        r['asignacion'] = 'acierto debil'
    elif not (LARGO_MIN <= razon <= LARGO_MAX) and not blanco_entero:
        r['asignacion'] = 'largo dispar'
    elif ambiguo:
        r['asignacion'] = 'ambigua'
    else:
        r['asignacion'] = 'por secuencia'

    if r['asignacion'] == 'por secuencia' and (
            m['cobertura'] < MIN_COBERTURA or not (LARGO_MIN <= razon <= LARGO_MAX)):
        rescatados.append((r['consulta'], m['proteina'], m['cobertura'], cob_s,
                           m['largo_q'], m['largo_s']))

    L = locus_de.get(m['proteina'], 'NA')
    r['locus'] = L
    g = genes.get(L)
    if g:
        r['producto'] = g['producto']
        r['n_genes_mismo_producto'] = str(n_producto[g['producto']])
        r['conjunto_flagelar'] = g['conjunto_flagelar']
        if g['analizable'] == 'TRUE':
            for c in ('log2FC_37vs25', 'padj_37vs25', 'log2FC_pH6vspH8', 'padj_pH6vspH8'):
                r[c] = g[c]
        else:
            for c in ('log2FC_37vs25', 'padj_37vs25', 'log2FC_pH6vspH8', 'padj_pH6vspH8'):
                r[c] = 'no analizable'
    filas.append(r)

with open(f_sal, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=COLS, delimiter='\t', lineterminator='\n')
    w.writeheader()
    w.writerows(filas)

CT = ['consulta', 'puesto', 'proteina_kc583', 'locus', 'identidad_pct',
      'cobertura_pct', 'largo_publicado_aa', 'largo_kc583_aa',
      'puntaje', 'producto']
with open(f_sal_todo, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=CT, delimiter='\t', lineterminator='\n')
    w.writeheader()
    w.writerows(filas_todo)

# ---- Resumen en pantalla ----
print()
if globals().get('rescatados'):
    print()
    print('=== ASIGNADAS POR COBERTURA DEL BLANCO (consulta mas larga que el blanco) ===')
    for c, pr, cq, cs, lq, ls in rescatados:
        print('  %-18s %-18s cobertura consulta %.0f%%, del blanco %.1f%% '
              '(%d aa consultados frente a %d en KC583)' % (c, pr, cq, cs, lq, ls))
print()
print('=== ASIGNACION ===')
for k, n in collections.Counter(r['asignacion'] for r in filas).most_common():
    print('  %-16s %d' % (k, n))

print()
print('=== CANDIDATOS ASIGNADOS POR SECUENCIA ===')
print('%-18s %-20s %9s %8s %7s %5s %9s %10s' %
      ('consulta', 'locus', 'identidad', '2o mejor', 'nombres',
       'flag', 'log2FC', 'padj'))
for r in filas:
    if r['asignacion'] != 'por secuencia':
        continue
    print('%-18s %-20s %8s%% %7s%% %7s %5s %9s %10s' %
          (r['consulta'], r['locus'], r['identidad_pct'],
           r['segundo_identidad_pct'], r['n_genes_mismo_producto'],
           'si' if r['conjunto_flagelar'] == 'TRUE' else 'no',
           r['log2FC_37vs25'], r['padj_37vs25']))

duda = [r for r in filas if r['asignacion'] != 'por secuencia']
if duda:
    print()
    print('=== REVISAR A MANO (no entran al manuscrito asi) ===')
    for r in duda:
        print('  %-18s %-14s mejor acierto: %s (%s%%, cobertura %s%%, '
              '%s aa publicados frente a %s aa en KC583)' %
              (r['consulta'], r['asignacion'], r['proteina_kc583'],
               r['identidad_pct'], r['cobertura_pct'],
               r['largo_publicado_aa'], r['largo_kc583_aa']))

riesgo = [r for r in filas if r['n_genes_mismo_producto'] not in ('NA', '1')]
if riesgo:
    print()
    print('=== ASIGNACIONES QUE EL NOMBRE NO HABRIA RESUELTO ===')
    for r in riesgo:
        print('  %-10s %-20s %s genes comparten "%s"' %
              (r['consulta'], r['locus'], r['n_genes_mismo_producto'], r['producto']))

# ---- Especificidad: flagelares frente a no flagelares ----
# Comprueba si la represion a 37 C es propia del aparato flagelar o un efecto
# general de la temperatura. El universo son los candidatos publicados que se
# pudieron asignar por secuencia; es una comparacion DENTRO de ese conjunto,
# pequeno y no elegido al azar, no una prueba sobre el genoma (esa es el
# hipergeometrico de la seccion 6).

def reprimido(r):
    if r['padj_37vs25'] in ('NA', 'no analizable'):
        return False
    return float(r['padj_37vs25']) < 0.05 and float(r['log2FC_37vs25']) < 0

def fisher_una_cola(a, na, b, nb):
    """Cola superior de la hipergeometrica para la tabla [a, na-a; b, nb-b]."""
    tot, k = na + nb, a + b
    return sum(comb(na, i) * comb(nb, k - i)
               for i in range(a, min(na, k) + 1)) / comb(tot, k)

utiles = [r for r in filas if r['asignacion'] == 'por secuencia']
fl = [r for r in utiles if r['conjunto_flagelar'] == 'TRUE']
nf = [r for r in utiles if r['conjunto_flagelar'] != 'TRUE']
# Deduplicacion por locus. Varios estudios proponen la misma proteina y el
# contraste debe contar genes distintos, no consultas: si tres trabajos
# proponen Pap31, Pap31 pesa una vez.
def un_locus(lista):
    vistos, salida, colapsados = set(), [], []
    for r in sorted(lista, key=lambda x: x['consulta']):
        if r['locus'] in vistos:
            colapsados.append((r['consulta'], r['locus']))
            continue
        vistos.add(r['locus']); salida.append(r)
    return salida, colapsados

fl_n, nf_n = len(fl), len(nf)
fl, col_fl = un_locus(fl)
nf, col_nf = un_locus(nf)
if col_fl or col_nf:
    print()
    print('=== CONSULTAS COLAPSADAS POR LOCUS (misma proteina en varios estudios) ===')
    for c, l in col_fl + col_nf:
        print('  %-22s %s' % (c, l))
    print('  flagelares %d -> %d, no flagelares %d -> %d' % (fl_n, len(fl), nf_n, len(nf)))

a, b = sum(map(reprimido, fl)), sum(map(reprimido, nf))

if fl and nf:
    p = fisher_una_cola(a, len(fl), b, len(nf))
    print()
    print('=== ESPECIFICIDAD: REPRESION A 37 C ===')
    print('  flagelares    : %d reprimidos de %d' % (a, len(fl)))
    print('  no flagelares : %d reprimidos de %d' % (b, len(nf)))
    for r in nf:
        if reprimido(r):
            print('      no flagelar reprimido: %-8s %-20s log2FC %s  padj %s'
                  % (r['consulta'], r['locus'], r['log2FC_37vs25'], r['padj_37vs25']))
    print('  Fisher exacto, una cola: p = %.3e' % p)

    f_esp = os.path.join(os.path.dirname(f_sal), 'especificidad.tsv')
    with open(f_esp, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['grupo', 'reprimidos_37vs25', 'total'])
        w.writerow(['flagelar', a, len(fl)])
        w.writerow(['no_flagelar', b, len(nf)])
        w.writerow(['fisher_una_cola_p', '%.6e' % p, ''])
    print('  Tabla : %s' % f_esp)
PY

echo
echo "Tabla  : $SAL/candidatos_publicados.tsv"
echo "Aciertos: $SAL/candidatos_blastp.tsv"
