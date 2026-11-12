# -*- coding: utf-8 -*-
"""
Identifica el ialB publicado en el genoma de KC583 y reporta su expresion.

KC583 tiene seis genes anotados como "invasion associated locus B family
protein". Elegir uno por abundancia seria una conjetura, asi que se empareja
por secuencia contra la entrada curada de Swiss-Prot P35641, que es el gen
clonado por Mitchell y Minnick 1995 y el que Coleman y Minnick 2003 usaron
como sonda en su ensayo de proteccion con nucleasa S1.

El resultado sirve de control de especificidad: ialB es un factor de invasion
muy expresado y, si no se mueve en ninguno de los dos ejes donde el regulon
flagelar si se mueve, la senal flagelar no es fisiologia global.

    python3 scripts/17_ialb.py

Salida: resultados/candidatos/ialB_identificacion.tsv
"""
import os, csv, sys, subprocess

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def r(*p): return os.path.join(RAIZ, *p)

CONSULTA = r("datos/referencia/ialB_publicado.faa")
BLANCO   = r("datos/referencia/protein.faa")
for a in (CONSULTA, BLANCO):
    if not os.path.exists(a):
        print("falta", a); sys.exit(1)

# evalue permisivo a proposito: interesa cuanto se separa el primero del segundo
cmd = ["blastp", "-query", CONSULTA, "-subject", BLANCO, "-evalue", "10",
       "-outfmt", "6 sseqid pident length qlen slen qcovs bitscore evalue"]
sal = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
hits = [l.split("\t") for l in sal.strip().split("\n") if l.strip()]
if not hits:
    print("BLASTp no devolvio ningun hit"); sys.exit(1)
hits.sort(key=lambda c: -float(c[6]))

l2p, p2l = {}, {}
for f in csv.reader(open(r("resultados/anotacion/locus_a_proteina.tsv")), delimiter="\t"):
    if len(f) >= 2:
        l2p[f[0]] = f[1]; p2l[f[1]] = f[0]

print("=== EMPAREJAMIENTO POR SECUENCIA, P35641 CONTRA EL PROTEOMA ===")
print("%-18s %-22s %8s %6s %6s %7s %9s" %
      ("proteina", "locus", "ident", "largo", "qcovs", "bits", "evalue"))
for c in hits[:6]:
    p = c[0].split("|")[-1] if "|" in c[0] else c[0]
    print("%-18s %-22s %8s %6s %6s %7s %9s" %
          (p, p2l.get(p, "?"), c[1], c[2], c[5], c[6], c[7]))

mejor = hits[0]
prot = mejor[0].split("|")[-1] if "|" in mejor[0] else mejor[0]
locus = p2l.get(prot)
ident, cob = float(mejor[1]), float(mejor[5])
if ident < 95.0 or cob < 95.0:
    print("\nEl mejor hit no llega a 95 por ciento de identidad y cobertura.")
    print("La identificacion no se puede sostener por secuencia."); sys.exit(1)
if len(hits) > 1:
    sep = float(mejor[6]) / float(hits[1][6])
    print("\nSeparacion contra el segundo hit: %.1f veces en bitscore" % sep)
else:
    print("\nUnico hit en todo el proteoma: no hay ambiguedad de paralogo.")
print("ialB publicado = %s (%s), %.1f%% identidad sobre %s aa" %
      (locus, prot, ident, mejor[2]))

# ---- expresion de los seis parologos de la familia ----------------------
sl = {}
csl = r("resultados/condiciones/solido_vs_liquido.csv")
if os.path.exists(csl):
    sl = {f["locus"]: f for f in csv.DictReader(open(csl))}

def n(v):
    try: return float(v)
    except (TypeError, ValueError): return float("nan")

familia = []
for f in csv.DictReader(open(r("resultados/genes/tabla_de_genes.csv"))):
    if "invasion associated locus B" not in f["producto"]:
        continue
    m = sl.get(f["locus"], {})
    familia.append({
        "locus": f["locus"], "proteina": f["proteina"],
        "es_ialB_publicado": "si" if f["locus"] == locus else "no",
        "TPM_Pl25": n(f["TPM_Pl25"]), "TPM_Pl30": n(f["TPM_Pl30"]),
        "TPM_Pl37": n(f["TPM_Pl37"]), "TPM_pH07": n(f["TPM_pH07"]),
        "log2FC_37vs30": n(f["log2FC_37vs30"]), "padj_37vs30": n(f["padj_37vs30"]),
        "log2FC_37vs25": n(f["log2FC_37vs25"]), "padj_37vs25": n(f["padj_37vs25"]),
        "log2FC_liq_vs_plac": n(m.get("log2FoldChange")),
        "padj_liq_vs_plac": n(m.get("padj")),
    })

print("\n=== LOS SEIS PARALOGOS DE LA FAMILIA ialB ===")
print("%-22s %-4s %9s %9s | %8s %8s | %8s %8s | %8s %8s" % (
      "locus", "ialB", "TPM_Pl30", "TPM_pH07", "37vs30", "padj",
      "37vs25", "padj", "liq/plac", "padj"))
for d in familia:
    print("%-22s %-4s %9.1f %9.1f | %+8.3f %8.3g | %+8.3f %8.3g | %+8.3f %8.3g" % (
        d["locus"], d["es_ialB_publicado"], d["TPM_Pl30"], d["TPM_pH07"],
        d["log2FC_37vs30"], d["padj_37vs30"], d["log2FC_37vs25"],
        d["padj_37vs25"], d["log2FC_liq_vs_plac"], d["padj_liq_vs_plac"]))

PADJ = 0.05
mov = [d["locus"] for d in familia if d["padj_37vs30"] < PADJ]
print("\nParalogos con cambio significativo de 30 a 37 C: %s" % (mov or "ninguno"))
print("Coleman y Minnick 2003 reportan caida de ialB de 30 a 37 C en caldo de")
print("recuperacion, cepa JB584, a los 30 minutos, por proteccion con S1.")
print("Estos datos son 2 h en agar, cepa KC583, por RNA-seq, y no la reproducen.")

cols = ["locus", "proteina", "es_ialB_publicado", "TPM_Pl25", "TPM_Pl30",
        "TPM_Pl37", "TPM_pH07", "log2FC_37vs30", "padj_37vs30",
        "log2FC_37vs25", "padj_37vs25", "log2FC_liq_vs_plac", "padj_liq_vs_plac"]
os.makedirs(r("resultados/candidatos"), exist_ok=True)
dest = r("resultados/candidatos/ialB_identificacion.tsv")
with open(dest, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for d in sorted(familia, key=lambda x: x["locus"]):
        w.writerow(d)
print("\nTabla:", dest)
