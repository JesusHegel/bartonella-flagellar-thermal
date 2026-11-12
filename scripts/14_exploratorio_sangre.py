# -*- coding: utf-8 -*-
"""
EXPLORATORIO. Posicion del conjunto flagelar en sangre humana (HB37, HBBG).

Estas condiciones estan FUERA del modelo: solo el 4,6 % de sus lecturas asignan
a CDS, se prepararon en otro centro (GENEWIZ), se enriquecieron con MICROBEnrich
y se deplecionaron de globina. No se calculan valores p ni cambios de expresion.

Lo unico que se mide es el RANGO de cada gen flagelar dentro de su propia
muestra. El rango es insensible a la profundidad y a la composicion de la
biblioteca, que es justamente lo que difiere entre matrices.

    python3 scripts/14_exploratorio_sangre.py

Salida: resultados/exploratorio/rangos_flagelares.tsv

RESULTADO (12 de noviembre de 2026): NO CONCLUYENTE, y por una razon concreta.

El metodo funciona donde los datos son buenos: en placa reproduce el hallazgo
conocido, con el conjunto flagelar alto a 25 C (percentil ~56) y bajo a 30 y
37 C (~19 y ~21).

En las condiciones excluidas el conjunto flagelar queda ALTO: HB37 ~62, HBBG
~65, PlBG ~68. Pero eso no se puede leer como biologia, porque PlBG es placa a
37 C igual que Pl37, y marca 68 frente a 21. Esas dos condiciones se diferencian
en la atmosfera sanguinea y en el CENTRO de preparacion de biblioteca: Pl37 fue
a Yale y PlBG a GENEWIZ (seccion 1 de CIFRAS). El centro esta perfectamente
confundido con la condicion en las 11 muestras excluidas.

No es artefacto de ceros: entre el 95 y el 99 % de los transcritos tienen senal
en todas las muestras, asi que los rangos estan bien definidos.

CONSECUENCIA: estas condiciones no sirven como evidencia a favor ni en contra.
Van en limitaciones, no como apoyo.

OBSERVACION APARTE, esta si dentro del modelo y sin confusion de centro: a 30 C
y pH 7,4, el conjunto flagelar marca ~19 en placa (Pl30) y ~52 en liquido
(pH07). Solo cambia solido frente a liquido. Merece un contraste formal con
DESeq2, posible porque las dos condiciones estan entre las 12 del modelo.
"""
import os, sys, csv, statistics

from rutas import r   # redirige resultados/, intermedios/ y figuras/ a $DIR_*

ORDEN = ["Pl25", "Pl30", "Pl37", "pH06", "pH07", "pH08", "PlBG", "HB37", "HBBG", "HUVE"]
EN_MODELO = {"Pl25", "Pl30", "Pl37", "pH06", "pH07", "pH08"}

locus2prot = {}
with open(r("resultados/anotacion/locus_a_proteina.tsv")) as f:
    for fila in csv.reader(f, delimiter="\t"):
        if len(fila) >= 2:
            locus2prot[fila[0]] = fila[1]

flagelares = set()
with open(r("resultados/anotacion/conjunto_flagelar.tsv")) as f:
    for fila in csv.reader(f, delimiter="\t"):
        if fila:
            flagelares.add(fila[0])
prot_flagelar = {locus2prot[l] for l in flagelares if l in locus2prot}
if len(prot_flagelar) != len(flagelares):
    print("AVISO: %d de %d genes flagelares sin proteina asociada"
          % (len(flagelares) - len(prot_flagelar), len(flagelares)))

muestras = []
with open(r("datos/muestras.tsv")) as f:
    next(f)
    for fila in csv.reader(f, delimiter="\t"):
        if fila:
            muestras.append((fila[0], fila[1]))

def prot_de(nombre):
    partes = nombre.split("_cds_")
    if len(partes) != 2:
        return None
    return "_".join(partes[1].split("_")[:2])

filas = []
for run, cond in muestras:
    ruta = r("cuantificacion", run, "quant.sf")
    if not os.path.exists(ruta):
        print("falta", ruta); sys.exit(1)
    tpm = {}
    with open(ruta) as f:
        next(f)
        for ln in f:
            c = ln.rstrip("\n").split("\t")
            p = prot_de(c[0])
            if p:
                tpm[p] = tpm.get(p, 0.0) + float(c[3])
    n = len(tpm)
    no_cero = sum(1 for v in tpm.values() if v > 0)
    orden = sorted(tpm, key=lambda k: -tpm[k])
    rango = {p: i + 1 for i, p in enumerate(orden)}
    rf = [rango[p] for p in prot_flagelar if p in rango]
    tf = [tpm[p] for p in prot_flagelar if p in tpm]
    filas.append({
        "run": run, "cond": cond, "en_modelo": "si" if cond in EN_MODELO else "no",
        "transcritos": n, "no_cero": no_cero,
        "pct_no_cero": round(100 * no_cero / n, 1), "n_flagelares": len(rf),
        "rango_mediano": statistics.median(rf),
        "percentil_mediano": round(100 * (1 - statistics.median(rf) / n), 1),
        "tpm_mediano_flagelar": round(statistics.median(tf), 1),
        "tpm_mediano_global": round(statistics.median(tpm.values()), 1),
    })

filas.sort(key=lambda d: (ORDEN.index(d["cond"]), d["run"]))
cols = ["run", "cond", "en_modelo", "transcritos", "no_cero", "pct_no_cero",
        "n_flagelares", "rango_mediano", "percentil_mediano",
        "tpm_mediano_flagelar", "tpm_mediano_global"]

os.makedirs(r("resultados/exploratorio"), exist_ok=True)
sal = r("resultados/exploratorio/rangos_flagelares.tsv")
with open(sal, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for d in filas:
        w.writerow(d)

print("%-13s %-6s %-4s %8s %9s %12s %9s" %
      ("run", "cond", "mod", "rangoMed", "pctilMed", "TPMmedFlag", "%noCero"))
for d in filas:
    print("%-13s %-6s %-4s %8.0f %8.1f%% %12.1f %8.1f%%"
          % (d["run"], d["cond"], d["en_modelo"], d["rango_mediano"],
             d["percentil_mediano"], d["tpm_mediano_flagelar"], d["pct_no_cero"]))

print("\nResumen por condicion (percentil mediano del conjunto flagelar):")
for c in ORDEN:
    v = [d["percentil_mediano"] for d in filas if d["cond"] == c]
    if v:
        print("  %-6s %-4s %s"
              % (c, "" if c in EN_MODELO else "(fuera)",
                 "  ".join("%.1f%%" % x for x in v)))
print("\nSalida:", sal)
print("EXPLORATORIO: sin valores p, pero la salida si entra a verificar.sh.")
