# -*- coding: utf-8 -*-
"""
Asociacion entre represion a 37 C y reconocimiento por suero de pacientes.

Usa las 22 proteinas de Dichter et al. 2021, que es el unico trabajo con
serologia de pacientes reales. La condicion de cada una esta en fuentes.tsv
como seroreactiva=si o seroreactiva=no, tomada de su texto: catorce reactivas
con suero de pacientes peruanos y el resto no.

La pregunta: una proteina reprimida a 37 C, se reconoce menos?

    python3 scripts/15_seroreactividad.py

Salida: resultados/candidatos/seroreactividad.tsv
"""
import os, csv, sys
from math import comb

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def r(*p): return os.path.join(RAIZ, *p)

PADJ, LFC = 0.05, 0.0          # mismo criterio que el script 13

def fisher_una_cola(a, n1, b, n2):
    """P(X >= a) hipergeometrica, cola superior."""
    tot, exitos = n1 + n2, a + b
    p = 0.0
    for k in range(a, min(n1, exitos) + 1):
        if exitos - k > n2:
            continue
        p += comb(n1, k) * comb(n2, exitos - k) / comb(tot, exitos)
    return p

sero = {}
with open(r("datos/candidatos_publicados/fuentes.tsv")) as f:
    for fila in csv.DictReader(f, delimiter="\t"):
        rep = fila.get("lo_que_reporta", "") or ""
        if "seroreactiva=si" in rep:
            sero[fila["consulta"]] = True
        elif "seroreactiva=no" in rep:
            sero[fila["consulta"]] = False
if not sero:
    print("No hay filas con seroreactiva= en fuentes.tsv"); sys.exit(1)

filas = []
with open(r("resultados/candidatos/candidatos_publicados.tsv")) as f:
    for x in csv.DictReader(f, delimiter="\t"):
        if x["consulta"] not in sero or x["asignacion"] != "por secuencia":
            continue
        try:
            lfc, padj = float(x["log2FC_37vs25"]), float(x["padj_37vs25"])
        except (ValueError, KeyError):
            continue
        filas.append({
            "consulta": x["consulta"], "locus": x["locus"],
            "flagelar": x["conjunto_flagelar"] == "TRUE",
            "log2FC": lfc, "padj": padj,
            "reprimida": lfc < LFC and padj < PADJ,
            "seroreactiva": sero[x["consulta"]],
        })

if len(filas) < 10:
    print("Solo %d proteinas de Dichter asignadas; revisar" % len(filas)); sys.exit(1)

sr = [x for x in filas if x["seroreactiva"]]
ns = [x for x in filas if not x["seroreactiva"]]
a, b = sum(x["reprimida"] for x in sr), sum(x["reprimida"] for x in ns)
p = fisher_una_cola(a, len(sr), b, len(ns))

print("Proteinas de Dichter asignadas por secuencia: %d" % len(filas))
print()
print("%-18s %-20s %6s %10s %10s %s" %
      ("consulta", "locus", "log2FC", "padj", "reprimida", "seroreactiva"))
for x in sorted(filas, key=lambda d: d["log2FC"]):
    print("%-18s %-20s %6.3f %10.3g %10s %s" %
          (x["consulta"], x["locus"], x["log2FC"], x["padj"],
           "si" if x["reprimida"] else "no",
           "si" if x["seroreactiva"] else "no"))

print()
print("=== REPRESION FRENTE A SEROREACTIVIDAD ===")
print("  reactivas con suero de pacientes : %d reprimidas de %d" % (a, len(sr)))
print("  no reactivas                     : %d reprimidas de %d" % (b, len(ns)))
print("  Fisher exacto, una cola: p = %.4g" % p)

sal = r("resultados/candidatos/seroreactividad.tsv")
with open(sal, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["grupo", "reprimidas_37vs25", "total"])
    w.writerow(["seroreactiva", a, len(sr)])
    w.writerow(["no_seroreactiva", b, len(ns)])
    w.writerow(["fisher_una_cola_p", "%.6e" % p, ""])
print("  Tabla : %s" % sal)
