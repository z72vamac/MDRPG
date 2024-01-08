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
# np.random.seed(2)
# 2102.2
# semilla 1120: las apariencias engañan

# datos.m = 5
# datos = Data([], m=datos.m, r=2, modo=4, tmax=120, alpha = True,
#              init=True,
#              show=True,
#              seed=2)
# datos.generar_grafos()
#
# grafos = datos.mostrar_datos()

def TDST(datos, ciclo):

    grafos = datos.mostrar_datos()
    T_index = range(datos.m + 2)
    T_index_prima = range(1, datos.m+1)
    T_index_primaprima = range(datos.m+1)

    result = []

    vD = 2
    vC = 1

    # Creamos el modelo8
    MODEL = gp.Model("TD-Stages")

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
    difgLit = MODEL.addVars(dgLit_index, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='difgLit')

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
    difgRit = MODEL.addVars(dgRit_index, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='difgRit')


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
    difgij = MODEL.addVars(dgij_index, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='difgij')

    # Variable continua no negativa pgij = zgij * dgij
    pgij_index = zgij_index

    pgij = MODEL.addVars(pgij_index, vtype=GRB.CONTINUOUS, lb=0.0, name='pgij')

    # Variable continua no negativa dRLt que indica la distancia entre el punto de recogida en la etapa t y el punto de
    # salida para la etapa t+1
    # dRLt_index = T_index_primaprima
    #
    # dRLt = MODEL.addVars(dRLt_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dRLt')
    # difRLt = MODEL.addVars(
    #     dRLt_index, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='dRLt')

    # Parametrizacion del ciclo
    muxLt = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'muxLt')
    muxRt = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'muxRt')

    # Variables binaria que es 1 si en la etapa t el punto xL se encuentra en el segmento j
    muLit = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.BINARY, name = 'muLit')

    # Variables binaria que es 1 si en la etapa t el punto xR se encuentra en el segmento j
    muRit = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.BINARY, name = 'muRit')

    minit = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'minit')
    maxit = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'maxit')
    entryit = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.BINARY, name = 'entryit')

    minit1 = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'minit1')
    maxit1 = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'maxit1')
    entryit1 = MODEL.addVars(ciclo.num_segmentos, T_index_primaprima, vtype = GRB.BINARY, name = 'entryit1')

    # Variable continua no negativa dLRt que indica la distancia que recorre el camión en la etapa t mientras el dron se mueve
    dLRt_index = T_index

    dLRt = MODEL.addVars(dLRt_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dLRt')

    zijLRt = MODEL.addVars(ciclo.num_segmentos, ciclo.num_segmentos, dLRt_index, vtype = GRB.BINARY, name = 'zijLRt')

    dijLRt = MODEL.addVars(ciclo.num_segmentos, ciclo.num_segmentos, dLRt_index, vtype = GRB.CONTINUOUS, name = 'dijLRt')

    pijLRt = MODEL.addVars(ciclo.num_segmentos, ciclo.num_segmentos, dLRt_index, vtype = GRB.CONTINUOUS, name = 'pijLRt')

    # Variable continua no negativa dRLt que indica la distancia que recorre el camión junto al dron en la etapa t
    dRLt_index = T_index_primaprima

    dRLt = MODEL.addVars(dRLt_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dRLt')

    zijRLt = MODEL.addVars(ciclo.num_segmentos, ciclo.num_segmentos, dRLt_index, vtype = GRB.BINARY, name = 'zijRLt')

    dijRLt = MODEL.addVars(ciclo.num_segmentos, ciclo.num_segmentos, dRLt_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dijRLt')

    pijRLt = MODEL.addVars(ciclo.num_segmentos, ciclo.num_segmentos, dRLt_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'pijRLt')

    # Variable que modela el producto de muxLt con muLit
    pLit = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pLit')

    # Variable que modela el producto de muxRt con muRit
    pRit = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pRit')

    # Variables que indica la posicion absoluta en la poligonal
    landaLt = MODEL.addVars(T_index, vtype = GRB.CONTINUOUS, lb = 0, name = 'landaLt')

    # Variables que indica la posicion absoluta en la poligonal
    landaRt = MODEL.addVars(T_index, vtype = GRB.CONTINUOUS, lb = 0, name = 'landaRt')

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

    Lgi = MODEL.addVars(Lgi_index, vtype=GRB.CONTINUOUS, name='Lgi')
    landagi = MODEL.addVars(
        landagi_index, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='landagi')

    # Variables difiliares para modelar el valor absoluto
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

    # MODEL.addConstrs(ugit.sum('*', i, '*') == 1 for i in range(datos.m))
    # MODEL.addConstrs(vgit.sum('*', i, '*') == 1 for g in range(datos.m))

    # De todos los segmentos hay que salir y entrar menos de aquel que se toma como entrada al grafo y como salida del grafo
    MODEL.addConstrs(mugi[g, i] - ugit.sum(g, i, '*') == zgij.sum(g, '*', i) for g, i in rhogi.keys())
    MODEL.addConstrs(mugi[g, i] - vgit.sum(g, i, '*') == zgij.sum(g, i, '*') for g, i in rhogi.keys())

    MODEL.addConstrs(ugit.sum(g, '*', t) - vgit.sum(g, '*', t) == 0 for t in T_index_prima for g in T_index_prima)

    MODEL.addConstrs(pgi[g, i] >= mugi[g, i] + alphagi[g, i] - 1 for g, i in rhogi.keys())
    MODEL.addConstrs(pgi[g, i] <= mugi[g, i] for g, i in rhogi.keys())
    MODEL.addConstrs(pgi[g, i] <= alphagi[g, i] for g, i in rhogi.keys())

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

    # for g in range(datos.m):
    #     MODEL.addConstr(sgi[g, grafos[g].aristas[0]] == 0)

    for g in T_index_prima:
        for i in grafos[g-1].aristas[0:]:
            MODEL.addConstr(sgi[g, i] >= 0)
            MODEL.addConstr(sgi[g, i] <= grafos[g-1].num_aristas - 1)


    # Restricciones de distancias y producto
    MODEL.addConstrs((difgLit[g, i, t, dim] >=   xLt[t, dim] - Rgi[g, i, dim]) for g, i, t, dim in difgLit.keys())
    MODEL.addConstrs((difgLit[g, i, t, dim] >= - xLt[t, dim] + Rgi[g, i, dim]) for g, i, t, dim in difgLit.keys())

    MODEL.addConstrs(difgLit[g, i, t, 0]*difgLit[g, i, t, 0] + difgLit[g, i, t, 1] * difgLit[g, i, t, 1] <= dgLit[g, i, t] * dgLit[g, i, t] for g, i, t in ugit.keys())

    # SmallM = 0
    # BigM = 10000
    SmallM = 0
    BigM = 10000
    # SmallM = 10000
    # BigM = 0
    # for g in T_index_prima:
    #     for v in grafos[g-1].V:
    #         for w in ciclo.V:
    #             BigM = max([np.linalg.norm(w - v), BigM])
    #             SmallM = min([np.linalg.norm(w - v), SmallM])


    #BigM += 5
    #BigM = max([np.linalg.norm(origin-grafos[g].V) for g in range(datos.m)])
    MODEL.addConstrs((pgLit[g, i, t] >= SmallM * ugit[g, i, t]) for g, i, t in ugit.keys())
    MODEL.addConstrs((pgLit[g, i, t] >= dgLit[g, i, t] - BigM * (1 - ugit[g, i, t])) for g, i, t in ugit.keys())

    MODEL.addConstrs((difgij[g, i, j, dim] >=   Lgi[g, i, dim] - Rgi[g, j, dim]) for g, i, j, dim in difgij.keys())
    MODEL.addConstrs((difgij[g, i, j, dim] >= - Lgi[g, i, dim] + Rgi[g, j, dim]) for g, i, j, dim in difgij.keys())

    MODEL.addConstrs(difgij[g, i, j, 0]*difgij[g, i, j, 0] + difgij[g, i, j, 1] * difgij[g, i, j, 1] <= dgij[g, i, j] * dgij[g, i, j] for g, i, j in dgij.keys())


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

    MODEL.addConstrs((difgRit[g, i, t, dim] >=   Lgi[g, i, dim] - xRt[t, dim]) for g, i, t, dim in difgRit.keys())
    MODEL.addConstrs((difgRit[g, i, t, dim] >= - Lgi[g, i, dim] + xRt[t, dim]) for g, i, t, dim in difgRit.keys())

    MODEL.addConstrs(difgRit[g, i, t, 0]*difgRit[g, i, t, 0] + difgRit[g, i, t, 1] * difgRit[g, i, t, 1] <= dgRit[g, i, t] * dgRit[g, i, t] for g, i, t in vgit.keys())


    # SmallM = 0
    #BigM = 10000
    MODEL.addConstrs((pgRit[g, i, t] >= SmallM * vgit[g, i, t]) for g, i, t in vgit.keys())
    MODEL.addConstrs((pgRit[g, i, t] >= dgRit[g, i, t] - BigM * (1 - vgit[g, i, t])) for g, i, t in vgit.keys())

    # MODEL.addConstrs((difRLt[t, dim] >= xRt[t, dim] - xLt[t + 1, dim])
    #                  for t in range(datos.m+1) for dim in range(2))
    # MODEL.addConstrs((difRLt[t, dim] >= - xRt[t, dim] + xLt[t + 1, dim])
    #                  for t in range(datos.m+1) for dim in range(2))
    # MODEL.addConstrs(difRLt[t, 0]*difRLt[t, 0] + difRLt[t, 1]
    #                  * difRLt[t, 1] <= dRLt[t] * dRLt[t] for t in range(datos.m+1))
    #
    # MODEL.addConstrs((difLRt[t, dim] >= xLt[t, dim] - xRt[t, dim])
    #                  for t, dim in difLRt.keys())
    # MODEL.addConstrs((difRLt[t, dim] >= - xLt[t, dim] + xRt[t, dim])
    #                  for t, dim in difLRt.keys())
    # MODEL.addConstrs(difLRt[t, 0]*difLRt[t, 0] + difLRt[t, 1]
    #                  * difLRt[t, 1] <= dLRt[t] * dLRt[t] for t in dLRt.keys())

    BigM= 10000
    MODEL.addConstrs((gp.quicksum(pgLit[g, i, t] for i in grafos[g-1].aristas) + pgij.sum(g, '*', '*') +  gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas) + gp.quicksum(pgRit[g, i, t] for i in grafos[g-1].aristas))/vD <= dLRt[t]/vC + BigM*(1- gp.quicksum(ugit[g, i, t] for i in grafos[g-1].aristas)) for t in T_index_prima for g in T_index_prima)

    # longitudes = []
    # for g in T_index_prima:
    #     longitudes.append(sum([grafos[g-1].A[i // 100 - 1, i % 100]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas]))

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
            MODEL.addConstr(pgi[g, i] >= grafos[g-1].A[first, second])
        MODEL.addConstr(maxgi[g, i] <= 1 - entrygi[g, i])
        MODEL.addConstr(mingi[g, i] <= entrygi[g, i])
        MODEL.addConstr(Rgi[g, i, 0] == rhogi[g, i] * grafos[g-1].V[first, 0] + (1 - rhogi[g, i]) * grafos[g-1].V[second, 0])
        MODEL.addConstr(Rgi[g, i, 1] == rhogi[g, i] * grafos[g-1].V[first, 1] + (1 - rhogi[g, i]) * grafos[g-1].V[second, 1])
        MODEL.addConstr(Lgi[g, i, 0] == landagi[g, i] * grafos[g-1].V[first, 0] + (1 - landagi[g, i]) * grafos[g-1].V[second, 0])
        MODEL.addConstr(Lgi[g, i, 1] == landagi[g, i] * grafos[g-1].V[first, 1] + (1 - landagi[g, i]) * grafos[g-1].V[second, 1])

    if not(datos.alpha):
        for g in T_index_prima:
            MODEL.addConstr(gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas) >= grafos[g-1].alpha*grafos[g-1].longitud)

    for i, t in entryit.keys():
        MODEL.addConstr(maxit[i, t] - minit[i, t] == muxRt[i, t] - muxLt[i, t])
        MODEL.addConstr(maxit[i, t] <= entryit[i, t])
        MODEL.addConstr(minit[i, t] <= 1 - entryit[i, t])

    for t in T_index:
        for dim in range(2):
            MODEL.addConstr(xLt[t, dim] == gp.quicksum(muLit[j, t]*ciclo.V[j][dim] + pLit[j, t]*(ciclo.V[j+1][dim] - ciclo.V[j][dim]) for j in range(ciclo.num_segmentos)))
            MODEL.addConstr(xRt[t, dim] == gp.quicksum(muRit[j, t]*ciclo.V[j][dim] + pRit[j, t]*(ciclo.V[j+1][dim] - ciclo.V[j][dim]) for j in range(ciclo.num_segmentos)))

            #MODEL.addConstr(xRt[t, j] == gp.quicksum(muxRt[v, t] * ciclo.V[v][j] for v in range(ciclo.num_puntos)))
        for j in range(ciclo.num_segmentos):
            MODEL.addConstr(landaLt[t] - j >= muxLt[j, t] - ciclo.num_segmentos * (1 - muLit[j, t]))
            MODEL.addConstr(landaLt[t] - j <= muxLt[j, t] + ciclo.num_segmentos * (1 - muLit[j, t]))
            MODEL.addConstr(landaRt[t] - j >= muxRt[j, t] - ciclo.num_segmentos * (1 - muRit[j, t]))
            MODEL.addConstr(landaRt[t] - j <= muxRt[j, t] + ciclo.num_segmentos * (1 - muRit[j, t]))

        for i in range(ciclo.num_segmentos):
            MODEL.addConstr(pLit[i, t] >= muLit[i, t] + muxLt[i, t] - 1)
            MODEL.addConstr(pLit[i, t] <= muxLt[i, t])
            MODEL.addConstr(pLit[i, t] <= muLit[i, t])
            MODEL.addConstr(pRit[i, t] >= muRit[i, t] + muxRt[i, t] - 1)
            MODEL.addConstr(pRit[i, t] <= muxRt[i, t])
            MODEL.addConstr(pRit[i, t] <= muRit[i, t])


        for i in range(ciclo.num_segmentos):
            for j in range(ciclo.num_segmentos):
                MODEL.addConstr(zijLRt[i, j, t] <= muLit[i, t])
                MODEL.addConstr(zijLRt[i, j, t] <= muRit[j, t])
                MODEL.addConstr(zijLRt[i, j, t] >= muLit[i, t] + muRit[j, t] - 1)

                if i == j:
                    MODEL.addConstr(dijLRt[i, j, t] == (maxit[i, t] + minit[j, t])*ciclo.longitudes[i])
                if i < j:
                    MODEL.addConstr(dijLRt[i, j, t] == (1 - muxLt[i, t])*ciclo.longitudes[i] + sum([ciclo.longitudes[k] for k in range(i+1, j)]) + muxRt[j, t]*ciclo.longitudes[j])
                if i > j:
                    MODEL.addConstr(dijLRt[i, j, t] == muxLt[i, t]*ciclo.longitudes[i] + sum([ciclo.longitudes[k] for k in range(j+1, i)]) + (1-muxRt[j, t])*ciclo.longitudes[j])

                MODEL.addConstr(pijLRt[i, j, t] <= ciclo.longitud* zijLRt[i, j, t])
                MODEL.addConstr(pijLRt[j, i, t] <= dijLRt[j, i, t])
                MODEL.addConstr(pijLRt[i, j, t] >= dijLRt[i, j, t] - ciclo.longitud*(1 - zijLRt[i, j, t]))

        MODEL.addConstr(dLRt[t] == gp.quicksum(pijLRt[j, k, t] for j in range(ciclo.num_segmentos) for k in range(ciclo.num_segmentos)))

    for i, t in entryit1.keys():
        MODEL.addConstr(maxit1[i, t] - minit1[i, t] == muxRt[i, t] - muxLt[i, t+1])
        MODEL.addConstr(maxit1[i, t] <= entryit1[i, t])
        MODEL.addConstr(minit1[i, t] <= 1 - entryit1[i, t])

    for t in T_index_primaprima:
        for i in range(ciclo.num_segmentos):
            for j in range(ciclo.num_segmentos):
                MODEL.addConstr(zijRLt[i, j, t] <= muRit[i, t])
                MODEL.addConstr(zijRLt[i, j, t] <= muLit[j, t+1])
                MODEL.addConstr(zijRLt[i, j, t] >= muRit[i, t] + muLit[j, t+1] - 1)

                if i == j:
                    MODEL.addConstr(dijRLt[i, j, t] == (maxit1[i, t] + minit1[j, t])*ciclo.longitudes[i])
                if i < j:
                    MODEL.addConstr(dijRLt[i, j, t] == (1-muxRt[i, t])*ciclo.longitudes[i] + sum([ciclo.longitudes[k] for k in range(i+1, j)]) + muxLt[j, t+1]*ciclo.longitudes[j])
                if i > j:
                    MODEL.addConstr(dijRLt[i, j, t] == muxRt[i, t]*ciclo.longitudes[i] + sum([ciclo.longitudes[k] for k in range(j+1, i)]) + (1-muxLt[j, t+1])*ciclo.longitudes[j])

                MODEL.addConstr(pijRLt[i, j, t] <= ciclo.longitud* zijRLt[i, j, t])
                MODEL.addConstr(pijRLt[j, i, t] <= dijRLt[j, i, t])
                MODEL.addConstr(pijRLt[i, j, t] >= dijRLt[i, j, t] - ciclo.longitud*(1 - zijRLt[i, j, t]))

        MODEL.addConstr(dRLt[t] == gp.quicksum(pijRLt[j, k, t] for j in range(ciclo.num_segmentos) for k in range(ciclo.num_segmentos)))

    MODEL.addConstrs(muLit.sum('*', t) == 1 for t in T_index)
    MODEL.addConstrs(muRit.sum('*', t) == 1 for t in T_index)

        # MODEL.addConstr(gp.quicksum((pRit[j, t] - pLit[j, t])*ciclo.longitudes[j] for j in range(ciclo.num_segmentos)) == dLRt[t])

        # MODEL.addConstr(muxRt[0, t] <= muRit[0, t])
        # MODEL.addConstrs(muxRt[k, t] <= muRit[k-1, t] + muRit[k, t] for k in range(1, ciclo.num_segmentos - 1))
        # MODEL.addConstr(muxRt[ciclo.num_segmentos-1, t] <= muRit[ciclo.num_segmentos-2, t])

    # MODEL.addConstrs(landaRt[t] >= landaLt[t] for t in T_index_prima)

    # MODEL.addConstr(dLRt.sum('*') <= ciclo.longitud)
    # MODEL.addConstrs(vjt[k, t] >= muLit[j, t+1] + muRit[i, t] - 1 for t in range(datos.m+1) for k in range(ciclo.num_segmentos) for i in range(ciclo.num_segmentos) for j in range(ciclo.num_segmentos) if i < j and i < k and k < j)
    # MODEL.addConstrs(landaLt[t+1] >= landaRt[t] for t in T_index_primaprima)
    #
    # Distancias
    # MODEL.addConstrs(gp.quicksum((muRit[j, t] -pRit[j, t])*ciclo.longitudes[j] for j in range(ciclo.num_segmentos))
    # + gp.quicksum(vjt[j, t]*ciclo.longitudes[j] for j in range(ciclo.num_segmentos))
    # + gp.quicksum(pLit[j, t + 1]*ciclo.longitudes[j] for j in range(ciclo.num_segmentos)) == dRLt[t] for t in range(datos.m+1))

    # Distancias
    # MODEL.addConstrs(gp.quicksum((muLit[j, t] - pLit[j, t])*ciclo.longitudes[j] for j in range(ciclo.num_segmentos))
    # + gp.quicksum(ujt[j, t]*ciclo.longitudes[j] for j in range(ciclo.num_segmentos))
    # + gp.quicksum(pRit[j, t]*ciclo.longitudes[j] for j in range(ciclo.num_segmentos)) == dLRt[t] for t in T_index_prima)





    # Distancias
    # MODEL.addConstrs(gp.quicksum(landaLt[t] - j*(muLit[j, t] - pLit[j, t])*ciclo.longitudes[j] for j in range(ciclo.num_segmentos))
    # + gp.quicksum(ujt[j, t]*ciclo.longitudes[j] for j in range(ciclo.num_segmentos))
    # + gp.quicksum((landaRt[t] - j*pRit[j, t])*ciclo.longitudes[j] for j in range(ciclo.num_segmentos)) == dLRt[t] for t in T_index_prima)

    #MODEL.addConstrs(gp.quicksum((landaLt[t] - j)*mu + ujt[j, t] + landaRt[t] - j)*ciclo.longitudes[j] for j in range(ciclo.num_segmentos)) == dLRt[t] for t in T_index_prima)

    # MODEL.addConstrs(ciclo.longitud*(landaRt[t] - landaLt[t])/ciclo.num_segmentos == dLRt[t] for t in T_index)

    # Origen y destino
    MODEL.addConstrs(xLt[0, dim] == ciclo.V[0][dim] for dim in range(2))
    MODEL.addConstrs(xRt[0, dim] == ciclo.V[0][dim] for dim in range(2))

    MODEL.addConstrs(xRt[datos.m+1, dim] == ciclo.V[0][dim] for dim in range(2))
    MODEL.addConstrs(xLt[datos.m+1, dim] == ciclo.V[0][dim] for dim in range(2))

    MODEL.update()

    objective = gp.quicksum(pgLit[g, i, t] + pgRit[g, i, t] for g, i, t in pgRit.keys()) + gp.quicksum(pgij[g, i, j] for g, i, j in pgij.keys()) + gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for g in T_index_prima for i in grafos[g-1].aristas) + gp.quicksum(3*dLRt[t] for t in dLRt.keys()) + gp.quicksum(3*dRLt[t] for t in dRLt.keys())
    # + gp.quicksum(1*dLRt[t] for t in dLRt.keys())

    # objective = gp.quicksum(dRLt[t] for t in dRLt.keys()) + gp.quicksum(dLRt[t] for t in dLRt.keys())

    MODEL.setObjective(objective, GRB.MINIMIZE)
    MODEL.Params.Threads = 6
    # MODEL.Params.NonConvex = 2
    MODEL.Params.timeLimit = datos.tmax

    MODEL.update()
    MODEL.write('lp_TDST.lp')
    #MODEL.write('modelo.lp')
    MODEL.optimize()

    if MODEL.Status == 3:
        MODEL.computeIIS()
        MODEL.write('casa.ilp')
        result =  [np.nan, np.nan, np.nan, np.nan]
        if datos.grid:
            result.append('Grid')
        else:
            result.append('Delauney')

        result.append('Stages')

        return result

    if MODEL.SolCount == 0:
        result =  [np.nan, np.nan, np.nan, np.nan]
        if datos.grid:
            result.append('Grid')
        else:
            result.append('Delauney')

        result.append('Stages')

        return result

    result = []

    result.append(MODEL.getAttr('MIPGap'))
    result.append(MODEL.Runtime)
    result.append(MODEL.getAttr('NodeCount'))
    result.append(MODEL.ObjVal)

    if datos.grid:
        result.append('Grid')
    else:
        result.append('Delauney')

    result.append('Stages')

    # fig, ax = plt.subplots()
    # plt.axis([0, 100, 0, 100])
    # ax.set_aspect('equal')

    # vals_u = MODEL.getAttr('x', ugit)
    # selected_u = gp.tuplelist((g, i, t) for g, i, t in vals_u.keys() if vals_u[g, i, t] > 0.5)
    # print('ugit = ' + str(selected_u))

    # vals_z = MODEL.getAttr('x', zgij)
    # selected_z = gp.tuplelist((g, i, j) for g, i, j in vals_z.keys() if vals_z[g, i, j] > 0.5)
    # print('zgij = ' + str(selected_z))

    # vals_p = MODEL.getAttr('x', pgij)
    # selected_p = gp.tuplelist(vals_p[(g, i, j)] for g, i, j in vals_p.keys() if vals_z[g, i, j] > 0.5)
    # print('pgij = ' + str(selected_p))


    # vals_v = MODEL.getAttr('x', vgit)
    # selected_v = gp.tuplelist((g, i, t) for g, i, t in vals_v.keys() if vals_v[g, i, t] > 0.5)
    # print('vgit = ' + str(selected_v))
    #
    #
    # vals_muLit = MODEL.getAttr('x', muLit)
    # selected_muLit = gp.tuplelist((j, t) for j, t in vals_muLit.keys() if vals_muLit[j, t] > 0.5)
    # print('muLit = ' + str(selected_muLit))
    #
    # vals_muRit = MODEL.getAttr('x', muRit)
    # selected_muRit = gp.tuplelist((j, t) for j, t in vals_muRit.keys() if vals_muRit[j, t] > 0.5)
    # print('muRit = ' + str(selected_muRit))
    #
    # print('landaLt' + str(MODEL.getAttr('x', landaLt)))
    # print('landaRt' + str(MODEL.getAttr('x', landaRt)))


    # vals_dRLt = MODEL.getAttr('x', dRLt)
    # print('dRLt = ' + str(vals_dRLt))

    # vals_dLRt = MODEL.getAttr('x', dLRt)
    # print('dLRt = ' + str(vals_dLRt))
    # muxLtprint = np.zeros((ciclo.num_puntos, datos.m+1))
    # muxRtprint = np.zeros((ciclo.num_puntos, datos.m+1))

    # print('entrygi = ' + str(MODEL.getAttr('x', entrygi)))

    # with np.printoptions(precision=3, suppress=True):
    #     print('muxLt = ' + str(np.array(MODEL.getAttr('x', muxLt).values()).reshape((ciclo.num_segmentos, datos.m+2))).replace('\n', '\n\t'))
    #     print('muLit = ' + str(np.array(MODEL.getAttr('x', muLit).values()).reshape((ciclo.num_segmentos, datos.m+2))).replace('\n', '\n\t'))
    #     print('pLit = ' + str(np.array(MODEL.getAttr('x', pLit).values()).reshape((ciclo.num_segmentos, datos.m+2))).replace('\n', '\n\t'))
    #     print('muxRt = ' + str(np.array(MODEL.getAttr('x', muxRt).values()).reshape((ciclo.num_segmentos, datos.m+2))).replace('\n', '\n\t'))
    #     print('muRit = ' + str(np.array(MODEL.getAttr('x', muRit).values()).reshape((ciclo.num_segmentos, datos.m+2))).replace('\n', '\n\t'))
    #     print('pRit = ' + str(np.array(MODEL.getAttr('x', pRit).values()).reshape((ciclo.num_segmentos, datos.m+2))).replace('\n', '\n\t'))
    # for t, v in pLit.keys():
    #     muxLtprint[t, v] = np.round(pLit[t, v].X, 2)
    #     muxRtprint[t, v] = np.round(pRit[t, v].X, 2)
    #
    #
    # print(str(muxLtprint.T))
    # print(str(muxRtprint.T))
    # print('muxLt = ' + str(MODEL.getAttr('x', muxLt)))
    # print('muxRt = ' + str(MODEL.getAttr('x', muxRt)))

    # muxLtprint = np.zeros((ciclo.num_puntos, datos.m+2))
    # muxRtprint = np.zeros((ciclo.num_puntos, datos.m+2))
    #
    #
    # for t, v in pLit.keys():
    #     muxLtprint[t, v] = muxLt[t, v].X
    #     muxRtprint[t, v] = muxRt[t, v].X
    #
    #
    # print(str(muxLtprint.T))
    # print(str(muxRtprint.T))

    # muxLtprint = np.zeros((ciclo.num_puntos, datos.m+1))
    # muxRtprint = np.zeros((ciclo.num_puntos, datos.m+1))
    #
    #
    # for t, v in pLit.keys():
    #     muxLtprint[t, v] = muLit[t, v].X
    #     muxRtprint[t, v] = muRit[t, v].X
    #
    #
    # print(str(muxLtprint.T))
    # print(str(muxRtprint.T))
    # print('muLit = ' + str(MODEL.getAttr('x', muLit)))
    # print('muRit = ' + str(MODEL.getAttr('x', muRit)))


    # muLitprint = np.zeros((ciclo.num_segmentos, datos.m+2))
    # muRitprint = np.zeros((ciclo.num_segmentos, datos.m+2))
    #
    # for v, t in muLit.keys():
    #     muLitprint[v, t] = np.round(muLit[v, t].X, 2)
    #     muRitprint[v, t] = np.round(muRit[v, t].X, 2)
    #
    #
    # print(str(muLitprint.T))
    # print(str(muRitprint.T))
    #
    # vals_dRLt = MODEL.getAttr('x', ujt)
    # print(vals_dRLt)
    # print(pgij)
    # ind = 0
    # path_C = []
    # paths_D = []
    #
    # #path_C.append(origin)
    # path_C.append([xLt[0, 0].X, xLt[0, 1].X])
    # for t in T_index_prima:
    #     #    if ind < datos.m:
    #     path_C.append([xLt[t, 0].X, xLt[t, 1].X])
    #     if ind < datos.m:
    #         path_D = []
    #         path_D.append([xLt[t, 0].X, xLt[t, 1].X])
    #         index_g = 0
    #         index_i = 0
    #         for g, i, ti in selected_u:
    #             if ti == t:
    #                 index_g = g
    #                 index_i = i
    #
    #         count = 0
    #         path_D.append([Rgi[index_g, index_i, 0].X, Rgi[index_g, index_i, 1].X])
    #         path_D.append([Lgi[index_g, index_i, 0].X, Lgi[index_g, index_i, 1].X])
    #         while count < grafos[index_g-1].num_aristas-1:
    #             for g, i, j in selected_z:
    #                 if index_g == g and index_i == i:
    #                     count += 1
    #                     index_i = j
    #                     path_D.append([Rgi[index_g, index_i, 0].X, Rgi[index_g, index_i, 1].X])
    #                     path_D.append([Lgi[index_g, index_i, 0].X, Lgi[index_g, index_i, 1].X])
    #
    #         ind += 1
    #         path_D.append([xRt[t, 0].X, xRt[t, 1].X])
    #     paths_D.append(path_D)
    #     path_C.append([xRt[t, 0].X, xRt[t, 1].X])
    #
    # path_C.append([xLt[datos.m+1, 0].X, xLt[datos.m+1, 1].X])
    #
    # for t in range(len(path_C) - 1):
    #     if t % 2 == 1 and t > 0:
    #         print('xL[' + str(int((t + 1)/2)) + '] = ' + str(path_C[t]))
    #     else:
    #         print('xR[' + str(int(t / 2)) + '] = ' + str(path_C[t]))
    #
    # fig, ax = plt.subplots()
    # plt.axis([0, 100, 0, 100])
    # ax.set_aspect('equal')
    #
    # for g, i in rhogi.keys():
    #     plt.plot(Rgi[g, i, 0].X, Rgi[g, i, 1].X, 'ko', markersize=0.5, color='cyan')
    #     plt.plot(Lgi[g, i, 0].X, Lgi[g, i, 1].X, 'ko', markersize=0.5, color='cyan')
    # #
    # # path_C = []
    # for t in T_index_prima:
    #     # path_C.append([xLt[t, 0].X, xLt[t, 1].X])
    #     # path_C.append([xRt[t, 0].X, xRt[t, 1].X])
    #     plt.plot(xLt[t, 0].X, xLt[t, 1].X, 'ko', alpha = 0.3, markersize=10, color='green')
    #     ax.annotate("L" + str(t), xy = (xLt[t, 0].X, xLt[t, 1].X))
    #     plt.plot(xRt[t, 0].X, xRt[t, 1].X, 'ko', markersize=5, color='blue')
    #     ax.annotate("R" + str(t), xy = (xRt[t, 0].X, xRt[t, 1].X))
    # # ax.add_artist(Polygon(path_C, fill=False, animated=False,
    # #               linestyle='-', alpha=1, color='blue'))
    #
    # for path in paths_D:
    #     ax.add_artist(Polygon(path, fill=False, closed=False,
    #                   animated=False, alpha=1, color='red'))
    #
    # # ax.add_artist(Polygon(path_D, fill=False, animated=False,
    # #               linestyle='dotted', alpha=1, color='red'))
    #
    # ax.add_artist(ciclo.artist)
    #
    # for g in range(datos.m):
    #     grafo = grafos[g]
    #     nx.draw(grafo.G, grafo.pos, node_size=10,
    #             node_color='black', alpha=0.3, edge_color='gray')
    #     # nx.draw_networkx_labels(grafo.G, grafo.pos)
    #
    # plt.savefig('contraintuitivo2.png')
    #
    # plt.show()
    print(result)
    return result
