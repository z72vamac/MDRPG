"""Tenemos un conjunto E de entornos y un conjunto de poligonales P de las que queremos recorrer un porcentaje alfa p . Buscamos
   un tour de mínima distancia que alterne poligonal-entorno y que visite todas las poligonales"""


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

# error en la semilla 11101
np.random.seed(1)
# 2102.2
# semilla 1120: las apariencias engañan

nG = 2
datos = Data([], m=nG, r=1, modo=4, tmax=120, alpha = True,
             init=True,
             show=True,
             seed=2)
datos.generar_grafos()

# grafo1
# A = [90, 10]
# B = [90, 90]
# C = [10, 90]
# D = [10, 10]
#
# V = np.array([A, B, C, D])
# Ar = np.zeros((len(V), len(V)))
#
# Ar[0, 1] = 0.7
# Ar[1, 2] = 0.6
# Ar[2, 3] = 0.3
# Ar[0, 3] = .5
#
# grafos = [Grafo(V, Ar, 0.5)]

grafos = datos.mostrar_datos()


T_index = range(datos.m + 2)
T_index_prima = range(1, datos.m+1)
T_index_primaprima = range(datos.m+1)

grafos_T = Data([], m = 1, r = 6, modo = 4, tmax=120, alpha = True,
                init = True,
                show = True,
                seed = 2)
grafos_T.generar_grafo()

grafo_T = grafos_T.mostrar_datos()[0]

vD = 1.5

vC = 1

# Creamos el modelo8
MODEL = gp.Model("TD-Graph")

# Variables que modelan las distancias
# Variable binaria ugit = 1 si en la etapa t entramos por el segmento sgi
ugit_index = []

for g in T_index_prima:
    for i in grafos[g-1].aristas:
        for t in T_index_prima:
            ugit_index.append((g, i, t))


ugit = MODEL.addVars(ugit_index, vtype=GRB.BINARY, name='ugit')

# Variable continua no negativa dgLit que indica la distancia desde el punto de lanzamiento hasta el segmento
# sgi.
dgLit_index = ugit_index

dgLit = MODEL.addVars(dgLit_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dgLit')
auxgLit = MODEL.addVars(
    dgLit_index, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='dgRit')

# Variable continua no negativa pgLit = ugit * dgLit
pgLit_index = ugit_index

pgLit = MODEL.addVars(pgLit_index, vtype=GRB.CONTINUOUS, lb=0.0, name='pgLit')


# Variable binaria vgit = 1 si en la etapa t salimos por el segmento sgi
vgit_index = ugit_index

vgit = MODEL.addVars(vgit_index, vtype=GRB.BINARY, name='vgit')

# Variable continua no negativa dgRit que indica la distancia desde el punto de salida del segmento sgi hasta el
# punto de recogida del camion
dgRit_index = ugit_index

dgRit = MODEL.addVars(dgRit_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dgRit')
auxgRit = MODEL.addVars(
    dgRit_index, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='dgRit')


# Variable continua no negativa pgRit = vgit * dgRit
pgRit_index = ugit_index

pgRit = MODEL.addVars(pgRit_index, vtype=GRB.CONTINUOUS, lb=0.0, name='pgRit')


# Variable binaria zgij = 1 si voy del segmento i al segmento j del grafo g.
zgij_index = []
sgi_index = []

for g in T_index_prima:
    for i in grafos[g-1].aristas:
        sgi_index.append((g, i))
        for j in grafos[g-1].aristas:
            if i != j:
                zgij_index.append((g, i, j))

zgij = MODEL.addVars(zgij_index, vtype=GRB.BINARY, name='zgij')
sgi = MODEL.addVars(sgi_index, vtype=GRB.CONTINUOUS, lb=0, name='sgi')

# Variable continua no negativa dgij que indica la distancia entre los segmentos i j en el grafo g.
dgij_index = zgij_index

dgij = MODEL.addVars(dgij_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dgij')
auxgij = MODEL.addVars(
    dgij_index, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='dgij')

# Variable continua no negativa pgij = zgij * dgij
pgij_index = zgij_index

pgij = MODEL.addVars(pgij_index, vtype=GRB.CONTINUOUS, lb=0.0, name='pgij')

# Variable continua no negativa dRLt que indica la distancia entre el punto de recogida en la etapa t y el punto de
# salida para la etapa t+1
# dRLt_index = T_index_primaprima
#
# dRLt = MODEL.addVars(dRLt_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dRLt')
# auxRLt = MODEL.addVars(
#     dRLt_index, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='dRLt')

# Parametrizacion del grafo
muxLt = MODEL.addVars(grafo_T.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'muxLt')
minet = MODEL.addVars(grafo_T.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'minet')
maxet = MODEL.addVars(grafo_T.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'maxet')
aet = MODEL.addVars(grafo_T.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0,  name = 'aet')

#aet = MODEL.addVars(grafo_T.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'aet')
entryet = MODEL.addVars(grafo_T.aristas, T_index, vtype = GRB.BINARY, name = 'entryet')

# Variables binaria que es 1 si en la etapa t el punto xL se encuentra en el segmento e
met = MODEL.addVars(grafo_T.aristas, T_index, vtype = GRB.BINARY, name = 'met')

mmuet = MODEL.addVars(grafo_T.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'mmuet')

# Variable binaria que toma el valor 1 cuando el camion entra al vertice i
bit = MODEL.addVars(grafo_T.num_puntos, T_index, vtype = GRB.BINARY, name = 'bit')

# Variable que modela el producto de muxLt con muLjt
pLeit = MODEL.addVars(grafo_T.aristas, grafo_T.num_puntos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pLeit')

# Variable binaria que toma el valor 1 cuando el camion atraviesa la arista e
qijt = MODEL.addVars(grafo_T.num_puntos, grafo_T.num_puntos, T_index, vtype = GRB.BINARY, name = 'qet')

muxRt = MODEL.addVars(grafo_T.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'muxRt')

# Variables binaria que es 1 si en la etapa t el punto xR se encuentra en el segmento e
uet = MODEL.addVars(grafo_T.aristas, T_index, vtype = GRB.BINARY, name = 'uet')

umuet = MODEL.addVars(grafo_T.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'umuet')

# Variable binaria que toma el valor 1 cuando el camion entra al vertice i
cit = MODEL.addVars(grafo_T.num_puntos, T_index, vtype = GRB.BINARY, name = 'cit')

# Variable que modela el producto de muxLt con muLjt
pReit = MODEL.addVars(grafo_T.aristas, grafo_T.num_puntos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pReit')


zeet = MODEL.addVars(grafo_T.aristas, grafo_T.aristas, T_index, vtype = GRB.BINARY, name = 'zeet')



deet = MODEL.addVars(grafo_T.aristas, grafo_T.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'deet')

peet = MODEL.addVars(grafo_T.aristas, grafo_T.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'peet')

qeet = MODEL.addVars(grafo_T.aristas, grafo_T.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'qeet')

# Variable continua no negativa dLRt que indica la distancia que recorre el camión en la etapa t mientras el dron se mueve
dLRt_index = T_index

dLRt = MODEL.addVars(dLRt_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dLRt')
auxLRt = MODEL.addVars(
    dLRt_index, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='dRLt')

# Variables que modelan los puntos de entrada o recogida
# xLt: punto de salida del dron del camion en la etapa t
xLt_index = []

for t in T_index:
    for dim in range(2):
        xLt_index.append((t, dim))

xLt = MODEL.addVars(xLt_index, vtype=GRB.CONTINUOUS, name='xLt')

# xRt: punto de recogida del dron del camion en la etapa t
xRt_index = []

for t in T_index:
    for dim in range(2):
        xRt_index.append((t, dim))

xRt = MODEL.addVars(xRt_index, vtype=GRB.CONTINUOUS, name='xRt')

# Rgi: punto de recogida del dron para el segmento sgi
Rgi_index = []
rhogi_index = []

for g in T_index_prima:
    for i in grafos[g-1].aristas:
        rhogi_index.append((g, i))
        for dim in range(2):
            Rgi_index.append((g, i, dim))

Rgi = MODEL.addVars(Rgi_index, vtype=GRB.CONTINUOUS, name='Rgi')
rhogi = MODEL.addVars(rhogi_index, vtype=GRB.CONTINUOUS,
                      lb=0.0, ub=1.0, name='rhogi')

# Lgi: punto de lanzamiento del dron del segmento sgi
Lgi_index = Rgi_index
landagi_index = rhogi_index

Lgi = MODEL.addVars(Lgi_index, vtype=GRB.CONTINUOUS, name='Rgi')
landagi = MODEL.addVars(
    landagi_index, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='rhogi')

# Variables auxiliares para modelar el valor absoluto
mingi = MODEL.addVars(rhogi_index, vtype=GRB.CONTINUOUS, lb=0.0, ub = 1.0, name='mingi')
maxgi = MODEL.addVars(rhogi_index, vtype=GRB.CONTINUOUS, lb=0.0, ub = 1.0, name='maxgi')
entrygi = MODEL.addVars(rhogi_index, vtype=GRB.BINARY, name='entrygi')
mugi = MODEL.addVars(rhogi_index, vtype = GRB.BINARY, name = 'mugi')
pgi = MODEL.addVars(rhogi_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pgi')
alphagi = MODEL.addVars(rhogi_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'alphagi')


MODEL.update()

# En cada etapa hay que visitar/salir un segmento de un grafo
MODEL.addConstrs(ugit.sum('*', '*', t) == 1 for t in T_index_prima)
MODEL.addConstrs(vgit.sum('*', '*', t) == 1 for t in T_index_prima)

# # Para cada grafo g, existe un segmento i y una etapa t donde hay que recoger al dron
MODEL.addConstrs(ugit.sum(g, '*', '*') == 1 for g in T_index_prima)
MODEL.addConstrs(vgit.sum(g, '*', '*') == 1 for g in T_index_prima)

# MODEL.addConstrs(ugit.sum('*', i, '*') == 1 for i in range(nG))
# MODEL.addConstrs(vgit.sum('*', i, '*') == 1 for g in range(nG))

# De todos los segmentos hay que salir y entrar menos de aquel que se toma como entrada al grafo y como salida del grafo
MODEL.addConstrs(mugi[g, i] - ugit.sum(g, i, '*') == zgij.sum(g, '*', i) for g, i, j in zgij.keys())
MODEL.addConstrs(mugi[g, i] - vgit.sum(g, i, '*') == zgij.sum(g, i, '*') for g, i, j in zgij.keys())

# MODEL.addConstrs(ugit.sum(g, i, '*') <= mugi[g, i] for g, i in rhogi.keys())
# MODEL.addConstrs(vgit.sum(g, i, '*') <= mugi[g, i] for g, i in rhogi.keys())

# MODEL.addConstrs(mugi[g, i] <= zgij.sum(g, '*', i) for g, i, j in zgij.keys())
# MODEL.addConstrs(mugi[g, j] == zgij.sum(g, '*', j) for g, i, j in zgij.keys())

MODEL.addConstrs(ugit.sum(g, '*', t) - vgit.sum(g, '*', t) == 0 for t in T_index_prima for g in T_index_prima)

# MODEL.addConstrs(  (ugit.sum(g, i, '*') + zgij.sum(g, '*', i))
#                  - (vgit.sum(g, i, '*') + zgij.sum(g, i, '*')) == 0 for g in T_index_prima for i in grafos[g-1].aristas)

#MODEL.addConstrs(maxgi[g, j] + mingi[g, j] <= mugi[g, j] for g, j in rhogi.keys())

#MODEL.addConstrs(mugi[g, i] >= zgij[g, i, j] for g, i, j in zgij.keys())
#MODEL.addConstrs(mugi[g, i] == zgij.sum(g, '*', i) for g, i, j in zgij.keys())

MODEL.addConstrs(pgi[g, i] >= mugi[g, i] + alphagi[g, i] - 1 for g, i in rhogi.keys())
MODEL.addConstrs(pgi[g, i] <= mugi[g, i] for g, i in rhogi.keys())
MODEL.addConstrs(pgi[g, i] <= alphagi[g, i] for g, i in rhogi.keys())

# MODEL.addConstrs(ugit.sum(g, '*', t) - vgit.sum(g, '*', t) == 0 for t in T_index_prima for g in T_index_prima)

# MODEL.addConstrs(  (ugit.sum(g, i, '*') + zgij.sum(g, '*', i))
#                  - (vgit.sum(g, i, '*') + zgij.sum(g, i, '*')) == 0 for g in T_index_prima for i in grafos[g-1].aristas)
# MODEL.addConstrs(1 - ugit[g, i, t] == zgij.sum(g, '*', i)
#                  for g, i, t in ugit.keys())
# MODEL.addConstrs(1 - ugit[g, i, t] == zgij.sum(g, i, '*')
#                  for g, i, t in vgit.keys())

# MODEL.addConstr(ugit[0, 101, 0] == 0)
# MODEL.addConstr(ugit[0, 101, 1] == 0)
# MODEL.addConstr(ugit.sum(4, '*', 1) == 1)

# Eliminación de subtours
for g in T_index_prima:
    for i in grafos[g-1].aristas[0:]:
        for j in grafos[g-1].aristas[0:]:
            if i != j:
                MODEL.addConstr(grafos[g-1].num_aristas - 1 >= (sgi[g, i] - sgi[g, j]) + grafos[g-1].num_aristas * zgij[g, i, j])

# for g in range(nG):
#     MODEL.addConstr(sgi[g, grafos[g].aristas[0]] == 0)

for g in T_index_prima:
    for i in grafos[g-1].aristas[0:]:
        MODEL.addConstr(sgi[g, i] >= 0)
        MODEL.addConstr(sgi[g, i] <= grafos[g-1].num_aristas - 1)


# Restricciones de distancias y producto
MODEL.addConstrs((auxgLit[g, i, t, dim] >=   xLt[t, dim] - Rgi[g, i, dim]) for g, i, t, dim in auxgLit.keys())
MODEL.addConstrs((auxgLit[g, i, t, dim] >= - xLt[t, dim] + Rgi[g, i, dim]) for g, i, t, dim in auxgLit.keys())

MODEL.addConstrs(auxgLit[g, i, t, 0]*auxgLit[g, i, t, 0] + auxgLit[g, i, t, 1] * auxgLit[g, i, t, 1] <= dgLit[g, i, t] * dgLit[g, i, t] for g, i, t in ugit.keys())

# SmallM = 10000
#
#
# BigM = 0
# for g in T_index_prima:
#     for v in grafos[g-1].V:
#         for w in grafo_T.V:
#             BigM = max([np.linalg.norm(w - v), BigM])
#             SmallM = min([np.linalg.norm(w - v), SmallM])
SmallM = 0
BigM = 10000

#BigM += 5
#BigM = max([np.linalg.norm(origin-grafos[g].V) for g in range(nG)])
MODEL.addConstrs((pgLit[g, i, t] >= SmallM * ugit[g, i, t]) for g, i, t in ugit.keys())
MODEL.addConstrs((pgLit[g, i, t] >= dgLit[g, i, t] - BigM * (1 - ugit[g, i, t])) for g, i, t in ugit.keys())

MODEL.addConstrs((auxgij[g, i, j, dim] >=   Lgi[g, i, dim] - Rgi[g, j, dim]) for g, i, j, dim in auxgij.keys())
MODEL.addConstrs((auxgij[g, i, j, dim] >= - Lgi[g, i, dim] + Rgi[g, j, dim]) for g, i, j, dim in auxgij.keys())

MODEL.addConstrs(auxgij[g, i, j, 0]*auxgij[g, i, j, 0] + auxgij[g, i, j, 1] * auxgij[g, i, j, 1] <= dgij[g, i, j] * dgij[g, i, j] for g, i, j in dgij.keys())


for g, i, j in zgij.keys():
    first_i = i // 100 - 1
    second_i = i % 100
    first_j = j // 100 - 1
    second_j = j % 100

    segm_i = Poligonal(np.array([[grafos[g-1].V[first_i, 0], grafos[g-1].V[first_i, 1]], [
                       grafos[g-1].V[second_i, 0], grafos[g-1].V[second_i, 1]]]), grafos[g-1].A[first_i, second_i])
    segm_j = Poligonal(np.array([[grafos[g-1].V[first_j, 0], grafos[g-1].V[first_j, 1]], [
                       grafos[g-1].V[second_j, 0], grafos[g-1].V[second_j, 1]]]), grafos[g-1].A[first_j, second_j])

    BigM_local = eM.estima_BigM_local(segm_i, segm_j)
    SmallM_local = eM.estima_SmallM_local(segm_i, segm_j)
    MODEL.addConstr((pgij[g, i, j] >= SmallM_local * zgij[g, i, j]))
    MODEL.addConstr((pgij[g, i, j] >= dgij[g, i, j] - BigM_local * (1 - zgij[g, i, j])))

MODEL.addConstrs((auxgRit[g, i, t, dim] >=   Lgi[g, i, dim] - xRt[t, dim]) for g, i, t, dim in auxgRit.keys())
MODEL.addConstrs((auxgRit[g, i, t, dim] >= - Lgi[g, i, dim] + xRt[t, dim]) for g, i, t, dim in auxgRit.keys())

MODEL.addConstrs(auxgRit[g, i, t, 0]*auxgRit[g, i, t, 0] + auxgRit[g, i, t, 1] * auxgRit[g, i, t, 1] <= dgRit[g, i, t] * dgRit[g, i, t] for g, i, t in vgit.keys())


# SmallM = 0
#BigM = 10000
MODEL.addConstrs((pgRit[g, i, t] >= SmallM * vgit[g, i, t]) for g, i, t in vgit.keys())
MODEL.addConstrs((pgRit[g, i, t] >= dgRit[g, i, t] - BigM * (1 - vgit[g, i, t])) for g, i, t in vgit.keys())

# MODEL.addConstrs((auxRLt[t, dim] >= xRt[t, dim] - xLt[t + 1, dim])
#                  for t in range(nG+1) for dim in range(2))
# MODEL.addConstrs((auxRLt[t, dim] >= - xRt[t, dim] + xLt[t + 1, dim])
#                  for t in range(nG+1) for dim in range(2))
# MODEL.addConstrs(auxRLt[t, 0]*auxRLt[t, 0] + auxRLt[t, 1]
#                  * auxRLt[t, 1] <= dRLt[t] * dRLt[t] for t in range(nG+1))
#
# MODEL.addConstrs((auxLRt[t, dim] >= xLt[t, dim] - xRt[t, dim])
#                  for t, dim in auxLRt.keys())
# MODEL.addConstrs((auxRLt[t, dim] >= - xLt[t, dim] + xRt[t, dim])
#                  for t, dim in auxLRt.keys())
# MODEL.addConstrs(auxLRt[t, 0]*auxLRt[t, 0] + auxLRt[t, 1]
#                  * auxLRt[t, 1] <= dLRt[t] * dLRt[t] for t in dLRt.keys())

MODEL.addConstrs((pgLit[g, i, t]
                  + pgij.sum(g, '*', '*') + (alphagi[g, i])*grafos[g-1].longaristas[i // 100 - 1, i % 100]
                  + pgRit[g, i, t])/vD <= dLRt[t]/vC for g, i, t in pgLit.keys())

longitudes = []
for g in T_index_prima:
    longitudes.append(sum([grafos[g-1].A[i // 100 - 1, i % 100]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas]))

# MODEL.addConstrs((pgLit.sum('*', '*', t) +
#                   pgij.sum(g, '*', '*') +
#                   ugit.sum(g, '*', '*')*longitudes[g-1] +
#                   pgRit.sum('*', '*', t))/vD <= dLRt[t]/vC for t in T_index_prima for g in T_index_prima)
# MODEL.addConstrs((dLRt[t]/vD <= 50) for t in T_index_prima)

for g, i in rhogi.keys():
    first = i // 100 - 1
    second = i % 100
    MODEL.addConstr(rhogi[g, i] - landagi[g, i] == maxgi[g, i] - mingi[g, i])
    MODEL.addConstr(maxgi[g, i] + mingi[g, i] == alphagi[g, i])
    if datos.alpha:
        MODEL.addConstr(pgi[g, i] == grafos[g-1].A[first, second])
    MODEL.addConstr(maxgi[g, i] <= 1 - entrygi[g, i])
    MODEL.addConstr(mingi[g, i] <= entrygi[g, i])
    MODEL.addConstr(Rgi[g, i, 0] == rhogi[g, i] * grafos[g-1].V[first, 0] + (1 - rhogi[g, i]) * grafos[g-1].V[second, 0])
    MODEL.addConstr(Rgi[g, i, 1] == rhogi[g, i] * grafos[g-1].V[first, 1] + (1 - rhogi[g, i]) * grafos[g-1].V[second, 1])
    MODEL.addConstr(Lgi[g, i, 0] == landagi[g, i] * grafos[g-1].V[first, 0] + (1 - landagi[g, i]) * grafos[g-1].V[second, 0])
    MODEL.addConstr(Lgi[g, i, 1] == landagi[g, i] * grafos[g-1].V[first, 1] + (1 - landagi[g, i]) * grafos[g-1].V[second, 1])

if not(datos.alpha):
    for g in T_index_prima:
        MODEL.addConstr(gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas) >= grafos[g-1].alpha*grafos[g-1].longitud)


for e, t in uet.keys():
    MODEL.addConstr(maxet[e, t] - minet[e, t] == muxRt[e, t] - muxLt[e, t])
    MODEL.addConstr(maxet[e, t] <= entryet[e, t])
    MODEL.addConstr(minet[e, t] <= 1 - entryet[e, t])
    # MODEL.addConstr(maxet[e, t] <= 1 - entryet[e, t])
    # MODEL.addConstr(minet[e, t] <= entryet[e, t])


    first = e // 100 - 1
    second = e % 100

    # MODEL.addConstr(gp.quicksum(met[e1, t] for e1 in grafo.T_aristas >= bit[first, t] + bit[second, t])
    # MODEL.addConstr(uet[e, t] >= cit[first, t] + cit[second, t])

    MODEL.addConstr(bit[first, t] + bit[second, t] >= gp.quicksum(zeet[e, e1, t] for e1 in grafo_T.aristas if e1 != e))
    # MODEL.addConstr(cit[first, t] + cit[second, t] >= gp.quicksum(zeet[e1, e, t] for e1 in grafo_T.aristas if e1 != e))


    for i in range(grafo_T.num_puntos):
        MODEL.addConstr(pLeit[e, i, t] >= muxLt[e, t] + bit[i, t] - 1)
        MODEL.addConstr(pLeit[e, i, t] <= muxLt[e, t])
        MODEL.addConstr(pLeit[e, i, t] <= bit[i, t])

        # MODEL.addConstr(pLeit[e, first, t] >= muxLt[e, t] + bit[first, t] - 1)
        # MODEL.addConstr(pLeit[e, first, t] <= muxLt[e, t])
        # MODEL.addConstr(pLeit[e, first, t] <= bit[first, t])
        # MODEL.addConstr(pLeit[e, second, t] >= muxLt[e, t] + bit[second, t] - 1)
        # MODEL.addConstr(pLeit[e, second, t] <= muxLt[e, t])
        # MODEL.addConstr(pLeit[e, second, t] <= bit[second, t])
        #
        MODEL.addConstr(pReit[e, i, t] >= muxRt[e, t] + cit[i, t] - 1)
        MODEL.addConstr(pReit[e, i, t] <= muxRt[e, t])
        MODEL.addConstr(pReit[e, i, t] <= cit[i, t])

        MODEL.addConstr(gp.quicksum(met[e, t] for e in grafo_T.aristas if e // 100 - 1 == i or e % 100 == i) >= bit[i, t])
        MODEL.addConstr(gp.quicksum(uet[e, t] for e in grafo_T.aristas if e // 100 - 1 == i or e % 100 == i) >= cit[i, t])

        # MODEL.addConstr(uet.sum('*', t) >= cit[i, t])
        #
        # MODEL.addConstr(pReit[e, first, t] >= muxRt[e, t] + cit[first, t] - 1)
        # MODEL.addConstr(pReit[e, first, t] <= muxRt[e, t])
        # MODEL.addConstr(pReit[e, first, t] <= cit[first, t])
        # MODEL.addConstr(pReit[e, second, t] >= muxRt[e, t] + cit[second, t] - 1)
        # MODEL.addConstr(pReit[e, second, t] <= muxRt[e, t])
        # MODEL.addConstr(pReit[e, second, t] <= cit[second, t])

    MODEL.addConstr(mmuet[e, t] >= met[e, t] + muxLt[e, t] - 1)
    MODEL.addConstr(mmuet[e, t] <= met[e, t])
    MODEL.addConstr(mmuet[e, t] <= muxLt[e, t])

    MODEL.addConstr(umuet[e, t] >= uet[e, t] + muxRt[e, t] - 1)
    MODEL.addConstr(umuet[e, t] <= uet[e, t])
    MODEL.addConstr(umuet[e, t] <= muxRt[e, t])

# MODEL.addConstr(cit[2, 1] == 1)
# for t in T_index:
#     suma_L1 = 0
#     suma_L2 = 0
#     suma_R1 = 0
#     suma_R2 = 0
#     for e in grafo_T.aristas:
#         first = e // 100 - 1
#         second = e % 100
#
#         suma_L1 += met[e, t]*grafo_T.V[first][0] + mmuet[e, t]*(grafo_T.V[second][0] - grafo_T.V[first][0])
#         suma_L2 += met[e, t]*grafo_T.V[first][1] + mmuet[e, t]*(grafo_T.V[second][1] - grafo_T.V[first][1])
#
#         suma_R1 += uet[e, t]*grafo_T.V[first][0] + umuet[e, t]*(grafo_T.V[second][0] - grafo_T.V[first][0])
#         suma_R2 += uet[e, t]*grafo_T.V[first][1] + umuet[e, t]*(grafo_T.V[second][1] - grafo_T.V[first][1])
#
#     MODEL.addConstr(xLt[t, 0] == suma_L1)
#     MODEL.addConstr(xLt[t, 1] == suma_L2)
#     MODEL.addConstr(xRt[t, 0] == suma_R1)
#     MODEL.addConstr(xRt[t, 1] == suma_R2)

for t in T_index:
    for dim in range(2):
        MODEL.addConstr(xLt[t, dim] == gp.quicksum(met[e, t]*grafo_T.V[e // 100 - 1][dim] + mmuet[e, t]*(grafo_T.V[e % 100][dim] - grafo_T.V[e // 100 - 1][dim]) for e in grafo_T.aristas))
        MODEL.addConstr(xRt[t, dim] == gp.quicksum(uet[e, t]*grafo_T.V[e // 100 - 1][dim] + umuet[e, t]*(grafo_T.V[e % 100][dim] - grafo_T.V[e // 100 - 1][dim]) for e in grafo_T.aristas))


for e1, e2, t in deet.keys():
    if e1 == e2:
        first_e1 = e1 // 100 - 1
        second_e1 = e1 % 100

        # MODEL.addConstr(deet[e1, e2, t] == (aet[e2, t])*grafo_T.longaristas[first_e1, second_e1])
        MODEL.addConstr(deet[e1, e2, t] == (maxet[e1, t] + minet[e1, t])*grafo_T.longaristas[first_e1, second_e1])

    else:
        first_e1 = e1 // 100 - 1
        second_e1 = e1 % 100

        first_e2 = e2 // 100 - 1
        second_e2 = e2 % 100

        MODEL.addConstr(deet[e1, e2, t] == (pLeit[e1, first_e1, t] + bit[second_e1, t] - pLeit[e1, second_e1, t])*grafo_T.longaristas[first_e1, second_e1] + gp.quicksum(qijt[grafo_T.aristas[ind]//100-1, grafo_T.aristas[ind]%100, t]*grafo_T.longitudes[ind] for ind in range(grafo_T.num_aristas)) + gp.quicksum(qijt[grafo_T.aristas[ind]%100, grafo_T.aristas[ind]//100-1, t]*grafo_T.longitudes[ind] for ind in range(grafo_T.num_aristas)) + (pReit[e2, first_e2, t] + cit[second_e2, t] - pReit[e2, second_e2, t])*grafo_T.longaristas[first_e2, second_e2])

        # MODEL.addConstr(deet[e1, e2, t] == (pLeit[e1, first_e1, t] + bit[second_e1, t] - pLeit[e1, second_e1, t])*grafo_T.longaristas[first_e1, second_e1] + gp.quicksum(qijt[grafo_T.aristas[ind]//100-1, grafo_T.aristas[ind]%100, t]*grafo_T.longitudes[ind] for ind in range(grafo_T.num_aristas)) + (pReit[e2, first_e2, t] + cit[second_e2, t] - pReit[e2, second_e2, t])*grafo_T.longaristas[first_e2, second_e2])

    MODEL.addConstr(peet[e1, e2, t] <= grafo_T.longitud* zeet[e1, e2, t])
    MODEL.addConstr(peet[e1, e2, t] <= deet[e1, e2, t])
    MODEL.addConstr(peet[e1, e2, t] >= 0)
    MODEL.addConstr(peet[e1, e2, t] >= deet[e1, e2, t] - grafo_T.longitud*(1 - zeet[e1, e2, t]))

    MODEL.addConstr(zeet[e1, e2, t] >= met[e1, t] + uet[e2, t] - 1)
    MODEL.addConstr(zeet[e1, e2, t] <= met[e1, t])
    MODEL.addConstr(zeet[e1, e2, t] <= uet[e2, t])


#    MODEL.addConstr(peet[e1, e2, t] >= 0)


MODEL.addConstrs(met.sum('*', t) == 1 for t in T_index)
MODEL.addConstrs(uet.sum('*', t) == 1 for t in T_index)


#MODEL.addConstrs(xLt[t+1, dim] == xRt[t, dim] for t in T_index_primaprima for dim in range(2))

MODEL.addConstrs(dLRt[t] == gp.quicksum(peet[e1, e2, t] for e1 in grafo_T.aristas for e2 in grafo_T.aristas) for t in T_index)

for t in T_index:
    for i in range(grafo_T.num_puntos):
        MODEL.addConstr(bit[i, t] + gp.quicksum(qijt[j, i, t] for j in range(grafo_T.num_puntos) if ((j+1) * 100 + i) in grafo_T.aristas or ((i+1) * 100 + j) in grafo_T.aristas) == gp.quicksum(qijt[i, j, t] for j in range(grafo_T.num_puntos) if ((j+1) * 100 + i) in grafo_T.aristas or ((i+1) * 100 + j) in grafo_T.aristas) + cit[i, t])
    # MODEL.addConstr(gp.quicksum((pRjt[j, t] - pLjt[j, t])*ciclo.longitudes[j] for j in range(ciclo.num_segmentos)) == dLRt[t])

    # MODEL.addConstr(muxRt[0, t] <= muRjt[0, t])
    # MODEL.addConstrs(muxRt[k, t] <= muRjt[k-1, t] + muRjt[k, t] for k in range(1, ciclo.num_segmentos - 1))
    # MODEL.addConstr(muxRt[ciclo.num_segmentos-1, t] <= muRjt[ciclo.num_segmentos-2, t])
#
# Distancias
# MODEL.addConstrs(gp.quicksum((muRjt[j, t] -pRjt[j, t])*ciclo.longitudes[j] for j in range(ciclo.num_segmentos))
# + gp.quicksum(vjt[j, t]*ciclo.longitudes[j] for j in range(ciclo.num_segmentos))
# + gp.quicksum(pLjt[j, t + 1]*ciclo.longitudes[j] for j in range(ciclo.num_segmentos)) == dRLt[t] for t in range(nG+1))

# Distancias
# MODEL.addConstrs(gp.quicksum((muLjt[j, t] - pLjt[j, t])*ciclo.longitudes[j] for j in range(ciclo.num_segmentos))
# + gp.quicksum(ujt[j, t]*ciclo.longitudes[j] for j in range(ciclo.num_segmentos))
# + gp.quicksum(pRjt[j, t]*ciclo.longitudes[j] for j in range(ciclo.num_segmentos)) == dLRt[t] for t in T_index_prima)





# Distancias
# MODEL.addConstrs(gp.quicksum(landaLt[t] - j*(muLjt[j, t] - pLjt[j, t])*ciclo.longitudes[j] for j in range(ciclo.num_segmentos))
# + gp.quicksum(ujt[j, t]*ciclo.longitudes[j] for j in range(ciclo.num_segmentos))
# + gp.quicksum((landaRt[t] - j*pRjt[j, t])*ciclo.longitudes[j] for j in range(ciclo.num_segmentos)) == dLRt[t] for t in T_index_prima)

#MODEL.addConstrs(gp.quicksum((landaLt[t] - j)*mu + ujt[j, t] + landaRt[t] - j)*ciclo.longitudes[j] for j in range(ciclo.num_segmentos)) == dLRt[t] for t in T_index_prima)

# MODEL.addConstrs(ciclo.longitud*(landaRt[t] - landaLt[t])/ciclo.num_segmentos == dLRt[t] for t in T_index)

# Origen y destino
MODEL.addConstrs(xLt[0, dim] == grafo_T.V[0][dim] for dim in range(2))
MODEL.addConstrs(xRt[0, dim] == grafo_T.V[0][dim] for dim in range(2))

MODEL.addConstrs(xRt[nG+1, dim] == grafo_T.V[0][dim] for dim in range(2))
MODEL.addConstrs(xLt[nG+1, dim] == grafo_T.V[0][dim] for dim in range(2))

MODEL.update()

objective = gp.quicksum(pgLit[g, i, t] + pgRit[g, i, t] for g, i, t in pgRit.keys()) + gp.quicksum(pgij[g, i, j] for g, i, j in pgij.keys()) + gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for g in T_index_prima for i in grafos[g-1].aristas) + gp.quicksum(dLRt[t] for t in T_index)

# objective = gp.quicksum(dLRt[t] for t in T_index)


MODEL.setObjective(objective, GRB.MINIMIZE)
MODEL.Params.Threads = 8
# MODEL.Params.NonConvex = 2
MODEL.Params.timeLimit = 150

MODEL.update()

MODEL.write('modelo.lp')
MODEL.optimize()

MODEL.write('solucion.json')
MODEL.update()

if MODEL.Status == 3:
    MODEL.computeIIS()
    MODEL.write('casa.ilp')

fig, ax = plt.subplots()
plt.axis([0, 100, 0, 100])
ax.set_aspect('equal')

vals_u = MODEL.getAttr('x', ugit)
selected_u = gp.tuplelist((g, i, t) for g, i, t in vals_u.keys() if vals_u[g, i, t] > 0.5)
print('ugit = ' + str(selected_u))

vals_z = MODEL.getAttr('x', zgij)
selected_z = gp.tuplelist((g, i, j) for g, i, j in vals_z.keys() if vals_z[g, i, j] > 0.5)
print('zgij = ' + str(selected_z))

# vals_p = MODEL.getAttr('x', pgij)
# selected_p = gp.tuplelist(vals_p[(g, i, j)] for g, i, j in vals_p.keys() if vals_z[g, i, j] > 0.5)
# print('pgij = ' + str(selected_p))


vals_v = MODEL.getAttr('x', vgit)
selected_v = gp.tuplelist((g, i, t) for g, i, t in vals_v.keys() if vals_v[g, i, t] > 0.5)
print('vgit = ' + str(selected_v))

print(dLRt)

vals_met = MODEL.getAttr('x', met)
selected_met = gp.tuplelist((e, t) for e, t in vals_met.keys() if vals_met[e, t] > 0.5)
print('met = ' + str(selected_met))

vals_uet = MODEL.getAttr('x', uet)
selected_uet = gp.tuplelist((e, t) for e, t in vals_uet.keys() if vals_uet[e, t] > 0.5)
print('uet = ' + str(selected_uet))

vals_bit = MODEL.getAttr('x', bit)
selected_bit = gp.tuplelist((i, t) for i, t in vals_bit.keys() if vals_bit[i, t] > 0.5)
print('bit = ' + str(selected_bit))

vals_qijt = MODEL.getAttr('x', qijt)
selected_qijt = gp.tuplelist((i, j, t) for i, j, t in vals_qijt.keys() if vals_qijt[i, j, t] > 0.5)
print('qijt = ' + str(selected_qijt))

vals_cit = MODEL.getAttr('x', cit)
selected_cit = gp.tuplelist((i, t) for i, t in vals_cit.keys() if vals_cit[i, t] > 0.5)
print('cit = ' + str(selected_cit))

print(deet)

vals_zeet = MODEL.getAttr('x', zeet)
selected_zeet = gp.tuplelist((e1, e2, t) for e1, e2, t in vals_zeet.keys() if vals_zeet[e1, e2, t] > 0.5)
print('zeet = ' + str(selected_zeet))

#print(zeet)
# print(muxLt)
# print(muxRt)
# print(met)
# print(uet)
# print(umuet)
# print(mmuet)
# print(aet)
# print(bit)
# print(qijt)
# print(cit)

# print(maxet)
# print(minet)
# print(peet)
# print(muxLt)
# print(muxRt)

# print([(i, j, k) for (i, j, k) in vals_zeet.keys() if vals_peet[i, j, k] > 0.5])
# print('Rgi = ' + str(Rgi))
# print('Lgi = ' + str(Lgi))

#
#
# vals_muLjt = MODEL.getAttr('x', muLjt)
# selected_muLjt = gp.tuplelist((j, t) for j, t in vals_muLjt.keys() if vals_muLjt[j, t] > 0.5)
# print('muLjt = ' + str(selected_muLjt))
#
# vals_muRjt = MODEL.getAttr('x', muRjt)
# selected_muRjt = gp.tuplelist((j, t) for j, t in vals_muRjt.keys() if vals_muRjt[j, t] > 0.5)
# print('muRjt = ' + str(selected_muRjt))
#
# print('landaLt' + str(MODEL.getAttr('x', landaLt)))
# print('landaRt' + str(MODEL.getAttr('x', landaRt)))


# vals_dRLt = MODEL.getAttr('x', dRLt)
# print('dRLt = ' + str(vals_dRLt))

# vals_dLRt = MODEL.getAttr('x', dLRt)
# print('dLRt = ' + str(vals_dLRt))
# muxLtprint = np.zeros((ciclo.num_puntos, nG+1))
# muxRtprint = np.zeros((ciclo.num_puntos, nG+1))

# print('entrygi = ' + str(MODEL.getAttr('x', entrygi)))

# with np.printoptions(precision=3, suppress=True):
#     print('muxLt = ' + str(np.array(MODEL.getAttr('x', muxLt).values()).reshape((ciclo.num_segmentos, nG+2))).replace('\n', '\n\t'))
#     print('muLjt = ' + str(np.array(MODEL.getAttr('x', muLjt).values()).reshape((ciclo.num_segmentos, nG+2))).replace('\n', '\n\t'))
#     print('pLjt = ' + str(np.array(MODEL.getAttr('x', pLjt).values()).reshape((ciclo.num_segmentos, nG+2))).replace('\n', '\n\t'))
#     print('muxRt = ' + str(np.array(MODEL.getAttr('x', muxRt).values()).reshape((ciclo.num_segmentos, nG+2))).replace('\n', '\n\t'))
#     print('muRjt = ' + str(np.array(MODEL.getAttr('x', muRjt).values()).reshape((ciclo.num_segmentos, nG+2))).replace('\n', '\n\t'))
#     print('pRjt = ' + str(np.array(MODEL.getAttr('x', pRjt).values()).reshape((ciclo.num_segmentos, nG+2))).replace('\n', '\n\t'))
# for t, v in pLjt.keys():
#     muxLtprint[t, v] = np.round(pLjt[t, v].X, 2)
#     muxRtprint[t, v] = np.round(pRjt[t, v].X, 2)
#
#
# print(str(muxLtprint.T))
# print(str(muxRtprint.T))
# print('muxLt = ' + str(MODEL.getAttr('x', muxLt)))
# print('muxRt = ' + str(MODEL.getAttr('x', muxRt)))

# muxLtprint = np.zeros((ciclo.num_puntos, nG+2))
# muxRtprint = np.zeros((ciclo.num_puntos, nG+2))
#
#
# for t, v in pLjt.keys():
#     muxLtprint[t, v] = muxLt[t, v].X
#     muxRtprint[t, v] = muxRt[t, v].X
#
#
# print(str(muxLtprint.T))
# print(str(muxRtprint.T))

# muxLtprint = np.zeros((ciclo.num_puntos, nG+1))
# muxRtprint = np.zeros((ciclo.num_puntos, nG+1))
#
#
# for t, v in pLjt.keys():
#     muxLtprint[t, v] = muLjt[t, v].X
#     muxRtprint[t, v] = muRjt[t, v].X
#
#
# print(str(muxLtprint.T))
# print(str(muxRtprint.T))
# print('muLjt = ' + str(MODEL.getAttr('x', muLjt)))
# print('muRjt = ' + str(MODEL.getAttr('x', muRjt)))


# muLjtprint = np.zeros((ciclo.num_segmentos, nG+2))
# muRjtprint = np.zeros((ciclo.num_segmentos, nG+2))
#
# for v, t in muLjt.keys():
#     muLjtprint[v, t] = np.round(muLjt[v, t].X, 2)
#     muRjtprint[v, t] = np.round(muRjt[v, t].X, 2)
#
#
# print(str(muLjtprint.T))
# print(str(muRjtprint.T))
#
# vals_dRLt = MODEL.getAttr('x', ujt)
# print(vals_dRLt)
#print(pgij)
ind = 0
path_C = []
paths_D = []

#path_C.append(origin)
path_C.append([xLt[0, 0].X, xLt[0, 1].X])
for t in T_index_prima:
    #    if ind < nG:
    path_C.append([xLt[t, 0].X, xLt[t, 1].X])
    if ind < nG:
        path_D = []
        path_D.append([xLt[t, 0].X, xLt[t, 1].X])
        index_g = 0
        index_i = 0
        for g, i, ti in selected_u:
            if ti == t:
                index_g = g
                index_i = i

        count = 0
        path_D.append([Rgi[index_g, index_i, 0].X, Rgi[index_g, index_i, 1].X])
        path_D.append([Lgi[index_g, index_i, 0].X, Lgi[index_g, index_i, 1].X])
        limite = sum([1 for g, i, j in selected_z if g == index_g])
        while count < limite:
            for g, i, j in selected_z:
                if index_g == g and index_i == i:
                    count += 1
                    index_i = j
                    path_D.append([Rgi[index_g, index_i, 0].X, Rgi[index_g, index_i, 1].X])
                    path_D.append([Lgi[index_g, index_i, 0].X, Lgi[index_g, index_i, 1].X])

        ind += 1
        path_D.append([xRt[t, 0].X, xRt[t, 1].X])
    paths_D.append(path_D)
    path_C.append([xRt[t, 0].X, xRt[t, 1].X])

path_C.append([xLt[nG+1, 0].X, xLt[nG+1, 1].X])

for t in range(len(path_C) - 1):
    if t % 2 == 1 and t > 0:
        print('xL[' + str(int((t + 1)/2)) + '] = ' + str(path_C[t]))
    else:
        print('xR[' + str(int(t / 2)) + '] = ' + str(path_C[t]))

for g, i in rhogi.keys():
    plt.plot(Rgi[g, i, 0].X, Rgi[g, i, 1].X, 'ko', markersize=0.5, color='cyan')
    plt.plot(Lgi[g, i, 0].X, Lgi[g, i, 1].X, 'ko', markersize=0.5, color='cyan')
#
# path_C = []
for t in T_index_prima:
    # path_C.append([xLt[t, 0].X, xLt[t, 1].X])
    # path_C.append([xRt[t, 0].X, xRt[t, 1].X])
    plt.plot(xLt[t, 0].X, xLt[t, 1].X, 'ko', alpha = 0.3, markersize=10, color='green')
    ax.annotate("L" + str(t), xy = (xLt[t, 0].X, xLt[t, 1].X))
    plt.plot(xRt[t, 0].X, xRt[t, 1].X, 'ko', markersize=5, color='blue')
    ax.annotate("R" + str(t), xy = (xRt[t, 0].X, xRt[t, 1].X))
# ax.add_artist(Polygon(path_C, fill=False, animated=False,
#               linestyle='-', alpha=1, color='blue'))

for path in paths_D:
    ax.add_artist(Polygon(path, fill=False, closed=False,
                  animated=False, alpha=1, color='red'))

# ax.add_artist(Polygon(path_D, fill=False, animated=False,
#               linestyle='dotted', alpha=1, color='red'))

nx.draw(grafo_T.G, grafo_T.pos, node_size=100, node_color='black', alpha=0.3, edge_color='gray')
nx.draw_networkx_labels(grafo_T.G, grafo_T.pos)

for g in range(nG):
    grafo = grafos[g]
    nx.draw(grafo.G, grafo.pos, node_size=10,
            node_color='black', alpha=0.3, edge_color='gray')
    # nx.draw_networkx_labels(grafo.G, grafo.pos)

plt.savefig('TD-alphaE.png')

plt.show()
