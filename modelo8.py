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

# np.random.seed(111111)
# error en la semilla 11101
# np.random.seed(1)
np.random.seed(1)

nG = 2
datos = Data([], m=nG, r=2, modo=4, tmax=120,
             init=True,
             show=True,
             seed=2)
datos.generar_muestra()

grafos = datos.mostrar_datos()
T_index = range(datos.m + 2)
T_index_prima = range(datos.m)

origin = [0, 0]
dest = [0, 0]

vD = 3
vC = 1

# Creamos el modelo8
MODEL = gp.Model("TD-Graph")

# Variables que modelan las distancias
# Variable binaria ugit = 1 si en la etapa t entramos por el segmento sgi
ugit_index = []

for g in range(datos.m):
    for i in grafos[g].aristas:
        for t in T_index_prima:
            ugit_index.append((g, i, t))


ugit = MODEL.addVars(ugit_index, vtype=GRB.BINARY, name='ugit')

# Variable continua no negativa dgLit que indica la distancia desde el punto de lanzamiento hasta el segmento
# sgi.
dgLit_index = ugit_index

dgLit = MODEL.addVars(dgLit_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dgLit')

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

# Variable continua no negativa pgRit = vgit * dgRit
pgRit_index = ugit_index

pgRit = MODEL.addVars(pgRit_index, vtype=GRB.CONTINUOUS, lb=0.0, name='pgRit')


# Variable binaria zgij = 1 si voy del segmento i al segmento j del grafo g.
zgij_index = []
sgi_index = []

for g in range(datos.m):
    for i in grafos[g].aristas:
        sgi_index.append((g, i))
        for j in grafos[g].aristas:
            if i != j:
                zgij_index.append((g, i, j))

zgij = MODEL.addVars(zgij_index, vtype=GRB.BINARY, name='zgij')
sgi = MODEL.addVars(sgi_index, vtype=GRB.CONTINUOUS, lb=0, name='sgi')

# Variable continua no negativa dgij que indica la distancia entre los segmentos i j en el grafo g.
dgij_index = zgij_index

dgij = MODEL.addVars(dgij_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dgij')

# Variable continua no negativa pgij = zgij * dgij
pgij_index = zgij_index

pgij = MODEL.addVars(pgij_index, vtype=GRB.CONTINUOUS, lb=0.0, name='pgij')

# Variable continua no negativa dRLt que indica la distancia entre el punto de recogida en la etapa t y el punto de
# salida para la etapa t+1
dRLt_index = T_index

dRLt = MODEL.addVars(dRLt_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dRLt')

# Variable continua no negativa dLRt que indica la distancia que recorre el camión en la etapa t mientras el dron se mueve
dLRt_index = T_index

dLRt = MODEL.addVars(dLRt_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dLRt')

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

for g in range(datos.m):
    for i in grafos[g].aristas:
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
mingi = MODEL.addVars(rhogi_index, vtype=GRB.CONTINUOUS, lb=0.0, name='mingi')
maxgi = MODEL.addVars(rhogi_index, vtype=GRB.CONTINUOUS, lb=0.0, name='maxgi')
entrygi = MODEL.addVars(rhogi_index, vtype=GRB.BINARY, name='entrygi')

MODEL.update()

# En cada etapa hay que visitar/salir un segmento de un grafo
MODEL.addConstrs(ugit.sum('*', '*', t) == 1 for t in T_index_prima)
MODEL.addConstrs(vgit.sum('*', '*', t) == 1 for t in T_index_prima)

# # Para cada grafo g, existe un segmento i y una etapa t donde hay que recoger al dron
MODEL.addConstrs(ugit.sum(g, '*', '*') == 1 for g in range(nG))
MODEL.addConstrs(vgit.sum(g, '*', '*') == 1 for g in range(nG))

# MODEL.addConstrs(ugit.sum('*', i, '*') == 1 for i in range(nG))
# MODEL.addConstrs(vgit.sum('*', i, '*') == 1 for g in range(nG))

# De todos los segmentos hay que salir y entrar menos de aquel que se toma como entrada al grafo y como salida del grafo
MODEL.addConstrs(1 - ugit.sum(g, i, '*') == zgij.sum(g, '*', i)
                 for g, i, j in zgij.keys())
MODEL.addConstrs(1 - vgit.sum(g, i, '*') == zgij.sum(g, i, '*')
                 for g, i, j in zgij.keys())

MODEL.addConstrs(ugit.sum(g, '*', t) - vgit.sum(g, '*', t)
                 == 0 for t in T_index_prima for g in T_index_prima)

MODEL.addConstrs((ugit.sum(g, i, '*') + zgij.sum(g, '*', i))
                 - (vgit.sum(g, i, '*') + zgij.sum(g, i, '*')) == 0 for g in range(nG) for i in grafos[g].aristas)
# MODEL.addConstrs(1 - ugit[g, i, t] == zgij.sum(g, '*', i)
#                  for g, i, t in ugit.keys())
# MODEL.addConstrs(1 - ugit[g, i, t] == zgij.sum(g, i, '*')
#                  for g, i, t in vgit.keys())

# MODEL.addConstr(ugit[0, 101, 0] == 0)
# MODEL.addConstr(ugit[0, 101, 1] == 0)


# Eliminación de subtours
for g in range(nG):
    for i in grafos[g].aristas[0:]:
        for j in grafos[g].aristas[0:]:
            if i != j:
                MODEL.addConstr(grafos[g].num_aristas - 1 >= (sgi[g, i]
                                - sgi[g, j]) + grafos[g].num_aristas * zgij[g, i, j])

# for g in range(nG):
#     MODEL.addConstr(sgi[g, grafos[g].aristas[0]] == 0)

for g in range(nG):
    for i in grafos[g].aristas[0:]:
        MODEL.addConstr(sgi[g, i] >= 0)
        MODEL.addConstr(sgi[g, i] <= grafos[g].num_aristas - 1)


# Restricciones de distancias y producto
MODEL.addConstrs((xLt[t+1, 0] - Rgi[g, i, 0]) * (xLt[t+1, 0] - Rgi[g, i, 0])
                 + (xLt[t+1, 1] - Rgi[g, i, 1]) * (xLt[t+1, 1] - Rgi[g, i, 1]) <= dgLit[g, i, t] * dgLit[g, i, t] for g, i, t in ugit.keys())

SmallM = 0
BigM = 10000

# BigM = 0
# for g in range(nG):
#     for v in grafos[g].V:
#         BigM = max([np.linalg.norm(origin - v), BigM])

#BigM = max([np.linalg.norm(origin-grafos[g].V) for g in range(nG)])
MODEL.addConstrs((pgLit[g, i, t] >= SmallM * ugit[g, i, t])
                 for g, i, t in ugit.keys())
MODEL.addConstrs((pgLit[g, i, t] >= dgLit[g, i, t]
                 - BigM * (1 - ugit[g, i, t])) for g, i, t in ugit.keys())

MODEL.addConstrs((Lgi[g, i, 0] - Rgi[g, j, 0]) * (Lgi[g, i, 0] - Rgi[g, j, 0])
                 + (Lgi[g, i, 1] - Rgi[g, j, 1]) * (Lgi[g, i, 1] - Rgi[g, j, 1]) <= dgij[g, i, j] * dgij[g, i, j] for g, i, j in dgij.keys())

for g, i, j in zgij.keys():
    first_i = i // 100 - 1
    second_i = i % 100
    first_j = j // 100 - 1
    second_j = j % 100

    segm_i = Poligonal(np.array([[grafos[g].V[first_i, 0], grafos[g].V[first_i, 1]], [
                       grafos[g].V[second_i, 0], grafos[g].V[second_i, 1]]]), grafos[g].A[first_i, second_i])
    segm_j = Poligonal(np.array([[grafos[g].V[first_j, 0], grafos[g].V[first_j, 1]], [
                       grafos[g].V[second_j, 0], grafos[g].V[second_j, 1]]]), grafos[g].A[first_j, second_j])

    BigM = eM.estima_BigM_local(segm_i, segm_j)
    SmallM = eM.estima_SmallM_local(segm_i, segm_j)
    MODEL.addConstr((pgij[g, i, j] >= SmallM * zgij[g, i, j]))
    MODEL.addConstr((pgij[g, i, j] >= dgij[g, i, j] - BigM
                    * (1 - zgij[g, i, j])))

MODEL.addConstrs((Lgi[g, i, 0] - xRt[t+1, 0]) * (Lgi[g, i, 0] - xRt[t+1, 0])
                 + (Lgi[g, i, 1] - xRt[t+1, 1]) * (Lgi[g, i, 1] - xRt[t+1, 1]) <= dgRit[g, i, t] * dgRit[g, i, t] for g, i, t in vgit.keys())

SmallM = 0
#BigM = 10000
MODEL.addConstrs((pgRit[g, i, t] >= SmallM * vgit[g, i, t])
                 for g, i, t in vgit.keys())
MODEL.addConstrs((pgRit[g, i, t] >= dgRit[g, i, t]
                 - BigM * (1 - vgit[g, i, t])) for g, i, t in vgit.keys())


MODEL.addConstrs((xRt[t, 0] - xLt[t + 1, 0]) * (xRt[t, 0] - xLt[t + 1, 0])
                 + (xRt[t, 1] - xLt[t + 1, 1]) * (xRt[t, 1] - xLt[t + 1, 1]) <= dRLt[t] * dRLt[t] for t in range(nG+1))

MODEL.addConstrs((xLt[t, 0] - xRt[t, 0]) * (xLt[t, 0] - xRt[t, 0])
                 + (xLt[t, 1] - xRt[t, 1]) * (xLt[t, 1] - xRt[t, 1]) <= dLRt[t] * dLRt[t] for t in dLRt.keys())

MODEL.addConstrs((pgLit[g, i, t]
                  + pgij.sum(g, '*', '*') + grafos[g].A[i // 100 - 1, i % 100]*grafos[g].longaristas[i // 100 - 1, i % 100]
                  + pgRit[g, i, t])/vD <= dLRt[t]/vC for g, i, t in pgLit.keys())
# MODEL.addConstrs((dLRt[t]/vD <= 50) for t in T_index_prima)

for g, i in rhogi.keys():
    first = i // 100 - 1
    second = i % 100
    MODEL.addConstr(rhogi[g, i] - landagi[g, i] == maxgi[g, i] - mingi[g, i])
    MODEL.addConstr(maxgi[g, i] + mingi[g, i] == grafos[g].A[first, second])
    MODEL.addConstr(maxgi[g, i] <= 1 - entrygi[g, i])
    MODEL.addConstr(mingi[g, i] <= entrygi[g, i])
    MODEL.addConstr(Rgi[g, i, 0] == rhogi[g, i] * grafos[g].V[first, 0] + (1 - rhogi[g, i]) * grafos[g].V[second, 0])
    MODEL.addConstr(Rgi[g, i, 1] == rhogi[g, i] * grafos[g].V[first, 1] + (1 - rhogi[g, i]) * grafos[g].V[second, 1])
    MODEL.addConstr(Lgi[g, i, 0] == landagi[g, i] * grafos[g].V[first, 0] + (1 - landagi[g, i]) * grafos[g].V[second, 0])
    MODEL.addConstr(Lgi[g, i, 1] == landagi[g, i] * grafos[g].V[first, 1] + (1 - landagi[g, i]) * grafos[g].V[second, 1])


# Origen y destino
MODEL.addConstrs(xLt[0, dim] == origin[dim] for dim in range(2))
MODEL.addConstrs(xRt[0, dim] == origin[dim] for dim in range(2))

MODEL.addConstrs(xRt[nG+1, dim] == dest[dim] for dim in range(2))
MODEL.addConstrs(xLt[nG+1, dim] == dest[dim] for dim in range(2))

MODEL.update()

objective = gp.quicksum(pgLit[g, i, t] + pgRit[g, i, t] for g, i, t in pgRit.keys()) + gp.quicksum(pgij[g, i, j] for g, i, j in pgij.keys()) + gp.quicksum(1*dLRt[t] + 1*dRLt[t] for t in T_index)

# objective = gp.quicksum(dRLt[t] + dLRt[t] for t in T_index)

MODEL.setObjective(objective, GRB.MINIMIZE)
MODEL.Params.Threads = 8
MODEL.Params.NonConvex = 2
MODEL.Params.timeLimit = 450

MODEL.update()

MODEL.write('modelo.lp')
MODEL.optimize()


MODEL.update()

if MODEL.Status == 3:
    MODEL.computeIIS()
    MODEL.write('casa.ilp')

fig, ax = plt.subplots()
plt.axis([0, 100, 0, 100])

vals_u = MODEL.getAttr('x', ugit)
selected_u = gp.tuplelist((g, i, t)
                          for g, i, t in vals_u.keys() if vals_u[g, i, t] > 0.5)
print(selected_u)

vals_z = MODEL.getAttr('x', zgij)
selected_z = gp.tuplelist((g, i, j)
                          for g, i, j in vals_z.keys() if vals_z[g, i, j] > 0.5)
print(selected_z)

vals_v = MODEL.getAttr('x', vgit)
selected_v = gp.tuplelist((g, i, t)
                          for g, i, t in vals_v.keys() if vals_v[g, i, t] > 0.5)
print(selected_v)


ind = 0
path_C = []
paths_D = []

#path_C.append(origin)
path_C.append([xLt[0, 0].X, xLt[0, 1].X])
for t in T_index_prima:
    #    if ind < nG:
    path_C.append([xLt[t+1, 0].X, xLt[t+1, 1].X])
    if ind < nG:
        path_D = []
        path_D.append([xLt[t+1, 0].X, xLt[t+1, 1].X])
        index_g = 0
        index_i = 0
        for g, i, ti in selected_u:
            if ti == t:
                index_g = g
                index_i = i
        count = 0
        while count < grafos[index_g].num_aristas-1:
            for g, i, j in selected_z:
                if index_g == g and index_i == i:
                    path_D.append([Rgi[g, i, 0].X, Rgi[g, i, 1].X])
                    path_D.append([Lgi[g, i, 0].X, Lgi[g, i, 1].X])
                    count += 1
                    index_i = j
        path_D.append([Rgi[index_g, index_i, 0].X, Rgi[index_g, index_i, 1].X])
        path_D.append([Lgi[index_g, index_i, 0].X, Lgi[index_g, index_i, 1].X])
        ind += 1
        path_D.append([xRt[t+1, 0].X, xRt[t+1, 1].X])
    paths_D.append(path_D)
    path_C.append([xRt[t+1, 0].X, xRt[t+1, 1].X])

path_C.append([xLt[nG+1, 0].X, xLt[nG+1, 1].X])

print(path_C)

for g, i in rhogi.keys():
    plt.plot(Rgi[g, i, 0].X, Rgi[g, i, 1].X, 'ko', markersize=2)
    plt.plot(Lgi[g, i, 0].X, Lgi[g, i, 1].X, 'ko', markersize=2)
#
# path_C = []
for t in T_index:
    # path_C.append([xLt[t, 0].X, xLt[t, 1].X])
    # path_C.append([xRt[t, 0].X, xRt[t, 1].X])
    plt.plot(xLt[t, 0].X, xLt[t, 1].X, 'ko', markersize=5, color='green')
    plt.plot(xRt[t, 0].X, xRt[t, 1].X, 'ko', markersize=5, color='blue')

ax.add_artist(Polygon(path_C, fill=False, animated=False,
              linestyle='-', alpha=1, color='black'))

for path in paths_D:
    ax.add_artist(Polygon(path, fill=False, animated=False,
                  linestyle='dotted', alpha=1, color='red'))

ax.add_artist(Polygon(path_D, fill=False, animated=False,
              linestyle='dotted', alpha=1, color='red'))

for g in range(nG):
    grafo = grafos[g]
    nx.draw(grafo.G, grafo.pos, node_size=20,
            node_color='r', alpha=0.5, edge_color='gray')

plt.savefig('modelo8a.png')

plt.show()
