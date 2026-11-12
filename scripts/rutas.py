# -*- coding: utf-8 -*-
# Carpetas del proyecto. Todos los scripts de python cargan este archivo.
# Equivalente de rutas.sh y rutas.R, con las mismas variables de entorno.
# Para la replica desde el SRA se cambian con variables de entorno (ver README).
import os

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if not os.path.exists(os.path.join(RAIZ, "CIFRAS_CONFIRMADAS.md")):
    raise SystemExit("ERROR: no encuentro la raiz del repositorio.")


def _v(nombre, por_defecto):
    return os.environ.get(nombre) or por_defecto


DIR_CUANT = _v("DIR_CUANT", "cuantificacion")
DIR_RES = _v("DIR_RES", "resultados")
DIR_INT = _v("DIR_INT", "intermedios")
DIR_FIG = _v("DIR_FIG", "figuras")

SUBCARPETAS = ["anotacion", "expresion_diferencial", "enriquecimiento", "robustez",
               "muestras", "genomica_comparada", "genes", "string", "candidatos",
               "condiciones", "exploratorio"]
for _d in [os.path.join(RAIZ, DIR_RES, _s) for _s in SUBCARPETAS] + \
          [os.path.join(RAIZ, DIR_INT), os.path.join(RAIZ, DIR_FIG)]:
    os.makedirs(_d, exist_ok=True)

_REDIR = (("resultados", DIR_RES), ("intermedios", DIR_INT), ("figuras", DIR_FIG))


def r(*p):
    """Ruta absoluta desde la raiz del repositorio.

    Las que empiezan en resultados/, intermedios/ o figuras/ se redirigen a
    DIR_RES, DIR_INT y DIR_FIG, para que la replica desde el SRA escriba en su
    propio arbol en vez de encima de los resultados publicados.
    """
    ruta = "/".join(p)
    for nombre, destino in _REDIR:
        if ruta == nombre or ruta.startswith(nombre + "/"):
            return os.path.join(RAIZ, destino, ruta[len(nombre):].lstrip("/"))
    return os.path.join(RAIZ, ruta)


def res(*p): return os.path.join(RAIZ, DIR_RES, *p)
def inter(*p): return os.path.join(RAIZ, DIR_INT, *p)
def figu(*p): return os.path.join(RAIZ, DIR_FIG, *p)
