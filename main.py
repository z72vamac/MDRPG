# Incluimos primero los paquetes
import gurobipy as gp
import pdb
from gurobipy import GRB
import numpy as np
from itertools import product
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.lines as mlines
from data import *
from entorno import *
import copy
import estimacion_M as eM
import networkx as nx
from PDST import PDST
from PDMTZ import PDMTZ
from PDSEC import PDSEC
from TDST import TDST
from TDMTZ import TDMTZ
from TDSEC import TDSEC
from NDMTZ import NDMTZ
from NDST import NDST
from PDMTZ_heuristic2 import PDMTZ_heuristic
from NDMTZ_heuristic2 import NDMTZ_heuristic
# from PDMTZ2 import *


# Definicion de los datos
""" P: conjunto de poligonales a agrupar
    E: conjunto de entornos
    T: sucesion de etapas
    C(e): centro del entorno e
    R(e): radio del entorno e
    p: indice de las poligonales
    e: indice de los entornos
    t: indice de las etapas
    n: dimension del problema
"""

# for i in range(8, 3200):
# i = 7
    # print(i)

# for i in range(0, 100):
np.random.seed(4)

lista = list(4*np.ones(6, np.int))
nG = len(lista)
datos = Data([], m=nG, tmax = 600, grid = False, alpha = False,
             init=True,
             show=True,
             seed=2)

datos.generar_grid()
datos.generar_grafos(lista)

grafo = Data([], m = 1, grid = True, tmax=600, alpha = True,
                init = True,
                show = True,
                seed = 2)
grafo.generar_grafo_personalizado(mode = 3)

# grafo = grafos_T.mostrar_datos()[0]

#
# ciclos = Data([], m = 1, r = 6, grid = False, tmax=120, alpha = True,
#                 init = True,
#                 show = True,
#                 seed = 2)
# ciclos.generar_ciclo()
#
# ciclo = ciclos.mostrar_datos()[0]
NDMTZ(datos, grafo)
# NDMTZ(datos, grafo)
# NDST(datos, grafo)
# PDMTZ(datos)

# PDMTZ_heuristic(datos)
# NDMTZ_heuristic(datos, grafo)
# TDST(datos, ciclo)
# TDMTZ(datos, ciclo)
# TDSEC(datos, ciclo)

# xL_dict = {}
# xR_dict = {}
#
# for g in range(nG+1):
#     xL_dict[g] = [xL[(g, 0)].X, xL[(g, 1)].X]
#     xR_dict[g] = [xR[(g, 0)].X, xR[(g, 1)].X]


# PDST(datos) #, xL_dict, xR_dict)
