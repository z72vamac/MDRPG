"""Tenemos un conjunto E de entornos y un conjunto de poligonales P de las que queremos recorrer un porcentaje alfa p . Buscamos
   un tour de mínima distancia que alterne poligonal-entorno y que visite todas las poligonales"""


# Incluimos primero los paqmuRetes
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

def NDST(datos, grafo_data):

    grafos = datos.mostrar_datos()

    T_index = range(datos.m + 2)
    T_index_prima = range(1, datos.m+1)
    T_index_primaprima = range(datos.m+1)

    vD = grafo_data.vD
    vC = grafo_data.vC

    orig = grafo_data.orig
    dest = grafo_data.dest

    grafo = grafo_data.data[0]

    # Creamos el modelo8
    MODEL = gp.Model("NDST")

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
        dgLit_index, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='difgLit')

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
        dgRit_index, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='difgRit')


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

    # ParamuLetrizacion del grafo
    muxLt = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'muxLt')
    mineLRt = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'mineLRt')
    maxeLRt = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'maxeLRt')

    mineRLt = MODEL.addVars(grafo.aristas, T_index_primaprima, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'mineRLt')
    maxeRLt = MODEL.addVars(grafo.aristas, T_index_primaprima, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'maxeRLt')

    #aet = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'aet')
    entryeLRt = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.BINARY, name = 'entryeLRt')
    entryeRLt = MODEL.addVars(grafo.aristas, T_index_primaprima, vtype = GRB.BINARY, name = 'entryeRLt')

    # Variables binaria que es 1 si en la etapa t el punto xL se encuentra en el segmento e
    muLet = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.BINARY, name = 'muLet')

    pLet = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'pLet')

    # Variable binaria que toma el valor 1 cuando el camion entra al vertice i
    biLRt = MODEL.addVars(grafo.num_puntos, T_index, vtype = GRB.BINARY, name = 'biLRt')
    biRLt = MODEL.addVars(grafo.num_puntos, T_index_primaprima, vtype = GRB.BINARY, name = 'biRLt')

    # Variable que modela el producto de muxLt con biLRt
    pxLbiLRt = MODEL.addVars(grafo.aristas, grafo.num_puntos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pxLbiLRt')
    pxLbiRLt = MODEL.addVars(grafo.aristas, grafo.num_puntos, T_index_primaprima, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pxLbiRLt')


    # Variable binaria que toma el valor 1 cuando el camion atraviesa la arista e
    qeLRt = MODEL.addVars(grafo.num_puntos, grafo.num_puntos, T_index, vtype = GRB.INTEGER, name = 'qeLRt')
    qeRLt = MODEL.addVars(grafo.num_puntos, grafo.num_puntos, T_index_primaprima, vtype = GRB.INTEGER, name = 'qeRLt')

    muxRt = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'muxRt')

    # Variables binaria que es 1 si en la etapa t el punto xR se encuentra en el segmento e
    muRet = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.BINARY, name = 'muRet')

    pRet = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'pRet')

    # Variable binaria que toma el valor 1 cuando el camion entra al vertice i
    ciLRt = MODEL.addVars(grafo.num_puntos, T_index, vtype = GRB.BINARY, name = 'ciLRt')
    ciRLt = MODEL.addVars(grafo.num_puntos, T_index_primaprima, vtype = GRB.BINARY, name = 'ciRLt')

    # Variable que modela el producto de muxLt con muLjt
    pxRciLRt = MODEL.addVars(grafo.aristas, grafo.num_puntos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pxRciLRt')
    pxRciRLt = MODEL.addVars(grafo.aristas, grafo.num_puntos, T_index_primaprima, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pxRciRLt')


    zeeLRt = MODEL.addVars(grafo.aristas, grafo.aristas, T_index, vtype = GRB.BINARY, name = 'zeeLRt')
    deeLRt = MODEL.addVars(grafo.aristas, grafo.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'deeLRt')
    peeLRt = MODEL.addVars(grafo.aristas, grafo.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'peeLRt')

    zeeRLt = MODEL.addVars(grafo.aristas, grafo.aristas, T_index_primaprima, vtype = GRB.BINARY, name = 'zeeRLt')
    deeRLt = MODEL.addVars(grafo.aristas, grafo.aristas, T_index_primaprima, vtype = GRB.CONTINUOUS, name = 'deeRLt')
    peeRLt = MODEL.addVars(grafo.aristas, grafo.aristas, T_index_primaprima, vtype = GRB.CONTINUOUS, name = 'peeRLt')

    # Variable continua no negativa dLRt que indica la distancia que recorre el camión en la etapa t mientras el dron se mueve
    dLRt_index = T_index

    dLRt = MODEL.addVars(dLRt_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dLRt')
    dRLt = MODEL.addVars(T_index_primaprima, vtype=GRB.CONTINUOUS, lb=0.0, name='dRLt')


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

    # MODEL.addConstrs(ugit.sum('*', i, '*') == 1 for i in range(datos.m))
    # MODEL.addConstrs(vgit.sum('*', i, '*') == 1 for g in range(datos.m))

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

    # for g in range(datos.m):
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
    # BigM = 0
    # for g in T_index_prima:
    #     for v in grafos[g-1].V:
    #         for w in grafo.V:
    #             BigM = max([np.linalg.norm(w - v), BigM])
    #             SmallM = min([np.linalg.norm(w - v), SmallM])
    SmallM = 0
    BigM = 10000

    #BigM += 5
    #BigM = max([np.linalg.norm(origin-grafos[g].V) for g in range(datos.m)])
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
    #                  for t in range(datos.m+1) for dim in range(2))
    # MODEL.addConstrs((auxRLt[t, dim] >= - xRt[t, dim] + xLt[t + 1, dim])
    #                  for t in range(datos.m+1) for dim in range(2))
    # MODEL.addConstrs(auxRLt[t, 0]*auxRLt[t, 0] + auxRLt[t, 1]
    #                  * auxRLt[t, 1] <= dRLt[t] * dRLt[t] for t in range(datos.m+1))
    #
    # MODEL.addConstrs((auxLRt[t, dim] >= xLt[t, dim] - xRt[t, dim])
    #                  for t, dim in auxLRt.keys())
    # MODEL.addConstrs((auxRLt[t, dim] >= - xLt[t, dim] + xRt[t, dim])
    #                  for t, dim in auxLRt.keys())
    # MODEL.addConstrs(auxLRt[t, 0]*auxLRt[t, 0] + auxLRt[t, 1]
    #                  * auxLRt[t, 1] <= dLRt[t] * dLRt[t] for t in dLRt.keys())

    BigM= 100000
    MODEL.addConstrs((gp.quicksum(pgLit[g, i, t] for i in grafos[g-1].aristas) + pgij.sum(g, '*', '*') +  gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas) + gp.quicksum(pgRit[g, i, t] for i in grafos[g-1].aristas))/vD <= dLRt[t]/vC + BigM*(1- gp.quicksum(ugit[g, i, t] for i in grafos[g-1].aristas)) for t in T_index_prima for g in T_index_prima)

    longitudes = []
    for g in T_index_prima:
        longitudes.append(sum([grafos[g-1].A[i // 100 - 1, i % 100]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas]))


    for g, i in rhogi.keys():
        first = i // 100 - 1
        second = i % 100
        MODEL.addConstr(rhogi[g, i] - landagi[g, i] == maxgi[g, i] - mingi[g, i])
        MODEL.addConstr(maxgi[g, i] <= 1 - entrygi[g, i])
        MODEL.addConstr(mingi[g, i] <= entrygi[g, i])
        MODEL.addConstr(maxgi[g, i] + mingi[g, i] == alphagi[g, i])
        if datos.alpha:
            MODEL.addConstr(pgi[g, i] >= grafos[g-1].A[first, second])
        MODEL.addConstr(Rgi[g, i, 0] == rhogi[g, i] * grafos[g-1].V[first, 0] + (1 - rhogi[g, i]) * grafos[g-1].V[second, 0])
        MODEL.addConstr(Rgi[g, i, 1] == rhogi[g, i] * grafos[g-1].V[first, 1] + (1 - rhogi[g, i]) * grafos[g-1].V[second, 1])
        MODEL.addConstr(Lgi[g, i, 0] == landagi[g, i] * grafos[g-1].V[first, 0] + (1 - landagi[g, i]) * grafos[g-1].V[second, 0])
        MODEL.addConstr(Lgi[g, i, 1] == landagi[g, i] * grafos[g-1].V[first, 1] + (1 - landagi[g, i]) * grafos[g-1].V[second, 1])

    if not(datos.alpha):
        for g in T_index_prima:
            MODEL.addConstr(gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas) >= grafos[g-1].alpha*grafos[g-1].longitud)

    # Restricciones para ir de L a R
    for e, t in muRet.keys():

        first = e // 100 - 1
        second = e % 100

        MODEL.addConstr(pLet[e, t] >= muLet[e, t] + muxLt[e, t] - 1)
        MODEL.addConstr(pLet[e, t] <= muLet[e, t])
        MODEL.addConstr(pLet[e, t] <= muxLt[e, t])

        MODEL.addConstr(pRet[e, t] >= muRet[e, t] + muxRt[e, t] - 1)
        MODEL.addConstr(pRet[e, t] <= muRet[e, t])
        MODEL.addConstr(pRet[e, t] <= muxRt[e, t])

        MODEL.addConstr(maxeLRt[e, t] - mineLRt[e, t] == muxRt[e, t] - muxLt[e, t])
        MODEL.addConstr(maxeLRt[e, t] <= entryeLRt[e, t])
        MODEL.addConstr(mineLRt[e, t] <= 1 - entryeLRt[e, t])

        # MODEL.addConstr(biLRt[first, t] + biLRt[second, t] >= gp.quicksum(zeeLRt[e, e1, t] for e1 in grafo.aristas if e1 != e))
        # MODEL.addConstr(ciLRt[first, t] + ciLRt[second, t] >= gp.quicksum(zeeLRt[e1, e, t] for e1 in grafo.aristas if e1 != e))
        MODEL.addConstr(biLRt[first, t] + biLRt[second, t] >= muLet[e, t])
        MODEL.addConstr(ciLRt[first, t] + ciLRt[second, t] >= muRet[e, t])

    for e, i, t in pxLbiLRt.keys():
        # for i in range(grafo.num_puntos):
        MODEL.addConstr(pxLbiLRt[e, i, t] >= muxLt[e, t] + biLRt[i, t] - 1)
        MODEL.addConstr(pxLbiLRt[e, i, t] <= muxLt[e, t])
        MODEL.addConstr(pxLbiLRt[e, i, t] <= biLRt[i, t])

        MODEL.addConstr(pxRciLRt[e, i, t] >= muxRt[e, t] + ciLRt[i, t] - 1)
        MODEL.addConstr(pxRciLRt[e, i, t] <= muxRt[e, t])
        MODEL.addConstr(pxRciLRt[e, i, t] <= ciLRt[i, t])

    for i, t in biLRt.keys():
        MODEL.addConstr(gp.quicksum(muLet[e, t] for e in grafo.aristas if e // 100 - 1 == i or e % 100 == i) >= biLRt[i, t])
        MODEL.addConstr(gp.quicksum(muRet[e, t] for e in grafo.aristas if e // 100 - 1 == i or e % 100 == i) >= ciLRt[i, t])

    for e1, e2, t in deeLRt.keys():
        if e1 == e2:
            first_e1 = e1 // 100 - 1
            second_e1 = e1 % 100

            MODEL.addConstr(deeLRt[e1, e2, t] == (maxeLRt[e1, t] + mineLRt[e1, t])*grafo.longaristas[first_e1, second_e1])

        else:
            first_e1 = e1 // 100 - 1
            second_e1 = e1 % 100

            first_e2 = e2 // 100 - 1
            second_e2 = e2 % 100

            MODEL.addConstr(deeLRt[e1, e2, t] == (pxLbiLRt[e1, first_e1, t] + biLRt[second_e1, t] - pxLbiLRt[e1, second_e1, t])*grafo.longaristas[first_e1, second_e1] + gp.quicksum(qeLRt[grafo.aristas[ind]//100-1, grafo.aristas[ind]%100, t]*grafo.longitudes[ind] for ind in range(grafo.num_aristas)) + gp.quicksum(qeLRt[grafo.aristas[ind]%100, grafo.aristas[ind]//100-1, t]*grafo.longitudes[ind] for ind in range(grafo.num_aristas)) + (pxRciLRt[e2, first_e2, t] + ciLRt[second_e2, t] - pxRciLRt[e2, second_e2, t])*grafo.longaristas[first_e2, second_e2])

        # MODEL.addConstr(peeLRt[e1, e2, t] <= 10000* zeeLRt[e1, e2, t])
        # MODEL.addConstr(peeLRt[e1, e2, t] <= deeLRt[e1, e2, t])
        MODEL.addConstr(peeLRt[e1, e2, t] >= SmallM * zeeLRt[e1, e2, t])
        MODEL.addConstr(peeLRt[e1, e2, t] >= deeLRt[e1, e2, t] - 10000*(1 - zeeLRt[e1, e2, t]))

        MODEL.addConstr(zeeLRt[e1, e2, t] >= muLet[e1, t] + muRet[e2, t] - 1)
        MODEL.addConstr(zeeLRt[e1, e2, t] <= muLet[e1, t])
        MODEL.addConstr(zeeLRt[e1, e2, t] <= muRet[e2, t])

    for t in T_index:
        for i in range(grafo.num_puntos):
            MODEL.addConstr(biLRt[i, t] + gp.quicksum(qeLRt[j, i, t] for j in range(grafo.num_puntos) if ((j+1) * 100 + i) in grafo.aristas or ((i+1) * 100 + j) in grafo.aristas) == gp.quicksum(qeLRt[i, j, t] for j in range(grafo.num_puntos) if ((j+1) * 100 + i) in grafo.aristas or ((i+1) * 100 + j) in grafo.aristas) + ciLRt[i, t])

    MODEL.addConstrs(dLRt[t] == gp.quicksum(peeLRt[e1, e2, t] for e1 in grafo.aristas for e2 in grafo.aristas) for t in T_index)

    # Restricciones para ir de R a L
    for e, t in entryeRLt.keys():
        # for t in T_index_primaprima:

        first = e // 100 - 1
        second = e % 100

        # if t < datos.m + 1:
        MODEL.addConstr(maxeRLt[e, t] - mineRLt[e, t] == muxRt[e, t] - muxLt[e, t+1])

        MODEL.addConstr(maxeRLt[e, t] <= entryeRLt[e, t])
        MODEL.addConstr(mineRLt[e, t] <= 1 - entryeRLt[e, t])

        # MODEL.addConstr(biRLt[first, t] + biRLt[second, t] >= gp.quicksum(zeeRLt[e, e1, t] for e1 in grafo.aristas if e1 != e))
        # MODEL.addConstr(ciRLt[first, t] + ciRLt[second, t] >= gp.quicksum(zeeRLt[e1, e, t] for e1 in grafo.aristas if e1 != e))

        MODEL.addConstr(biRLt[first, t] + biRLt[second, t] >= muRet[e, t])
        MODEL.addConstr(ciRLt[first, t] + ciRLt[second, t] >= muLet[e, t+1])


        # MODEL.addConstr(biRLt[first, t] + biRLt[second, t] <= zeeRLt[e, e, t])
        # MODEL.addConstr(ciRLt[first, t] + ciRLt[second, t] <= zeeRLt[e, e, t])

    for e, i, t in pxRciRLt.keys():
        # for i in range(grafo.num_puntos):
        MODEL.addConstr(pxRciRLt[e, i, t] >= muxRt[e, t] + biRLt[i, t] - 1)
        MODEL.addConstr(pxRciRLt[e, i, t] <= muxRt[e, t])
        MODEL.addConstr(pxRciRLt[e, i, t] <= biRLt[i, t])

        # if t < datos.m + 1:

        MODEL.addConstr(pxLbiRLt[e, i, t] >= muxLt[e, t+1] + ciRLt[i, t] - 1)
        MODEL.addConstr(pxLbiRLt[e, i, t] <= muxLt[e, t+1])
        MODEL.addConstr(pxLbiRLt[e, i, t] <= ciRLt[i, t])

    for i, t in biRLt.keys():
        MODEL.addConstr(gp.quicksum(muRet[e, t] for e in grafo.aristas if e // 100 - 1 == i or e % 100 == i) >= biRLt[i, t])
        MODEL.addConstr(gp.quicksum(muLet[e, t] for e in grafo.aristas if e // 100 - 1 == i or e % 100 == i) >= ciRLt[i, t])


    for t in T_index:
        for dim in range(2):
            MODEL.addConstr(xLt[t, dim] == gp.quicksum(muLet[e, t]*grafo.V[e // 100 - 1][dim] + pLet[e, t]*(grafo.V[e % 100][dim] - grafo.V[e // 100 - 1][dim]) for e in grafo.aristas))
            MODEL.addConstr(xRt[t, dim] == gp.quicksum(muRet[e, t]*grafo.V[e // 100 - 1][dim] + pRet[e, t]*(grafo.V[e % 100][dim] - grafo.V[e // 100 - 1][dim]) for e in grafo.aristas))

    MODEL.addConstrs(muLet.sum('*', t) == 1 for t in T_index)
    MODEL.addConstrs(muRet.sum('*', t) == 1 for t in T_index)

    for e1, e2, t in deeRLt.keys():
        if e1 == e2:
            first_e1 = e1 // 100 - 1
            second_e1 = e1 % 100

            MODEL.addConstr(deeRLt[e1, e2, t] == (maxeRLt[e1, t] + mineRLt[e1, t])*grafo.longaristas[first_e1, second_e1])

        else:
            first_e1 = e1 // 100 - 1
            second_e1 = e1 % 100

            first_e2 = e2 // 100 - 1
            second_e2 = e2 % 100

            MODEL.addConstr(deeRLt[e1, e2, t] == (pxRciRLt[e1, first_e1, t] + biRLt[second_e1, t] - pxRciRLt[e1, second_e1, t])*grafo.longaristas[first_e1, second_e1] + gp.quicksum(qeRLt[grafo.aristas[ind]//100-1, grafo.aristas[ind]%100, t]*grafo.longitudes[ind] for ind in range(grafo.num_aristas)) + gp.quicksum(qeRLt[grafo.aristas[ind]%100, grafo.aristas[ind]//100-1, t]*grafo.longitudes[ind] for ind in range(grafo.num_aristas)) + (pxLbiRLt[e2, first_e2, t] + ciRLt[second_e2, t] - pxLbiRLt[e2, second_e2, t])*grafo.longaristas[first_e2, second_e2])

        # MODEL.addConstr(peeRLt[e1, e2, t] <= 10000* zeeRLt[e1, e2, t])
        # MODEL.addConstr(peeRLt[e1, e2, t] <= deeRLt[e1, e2, t])
        MODEL.addConstr(peeRLt[e1, e2, t] >= SmallM * zeeRLt[e1, e2, t])
        MODEL.addConstr(peeRLt[e1, e2, t] >= deeRLt[e1, e2, t] - 10000*(1 - zeeRLt[e1, e2, t]))

        MODEL.addConstr(zeeRLt[e1, e2, t] >= muRet[e1, t] + muLet[e2, t+1] - 1)
        MODEL.addConstr(zeeRLt[e1, e2, t] <= muRet[e1, t])
        MODEL.addConstr(zeeRLt[e1, e2, t] <= muLet[e2, t+1])

    for t in T_index_primaprima:
        for i in range(grafo.num_puntos):
            MODEL.addConstr(biRLt[i, t] + gp.quicksum(qeRLt[j, i, t] for j in range(grafo.num_puntos) if ((j+1) * 100 + i) in grafo.aristas or ((i+1) * 100 + j) in grafo.aristas) == gp.quicksum(qeRLt[i, j, t] for j in range(grafo.num_puntos) if ((j+1) * 100 + i) in grafo.aristas or ((i+1) * 100 + j) in grafo.aristas) + ciRLt[i, t])

    MODEL.addConstrs(dRLt[t] == gp.quicksum(peeRLt[e1, e2, t] for e1 in grafo.aristas for e2 in grafo.aristas) for t in T_index_primaprima)


    # Origen y destino
    MODEL.addConstrs(xLt[0, dim] == orig[dim] for dim in range(2))
    MODEL.addConstrs(xRt[0, dim] == orig[dim] for dim in range(2))

    MODEL.addConstrs(xRt[datos.m+1, dim] == dest[dim] for dim in range(2))
    MODEL.addConstrs(xLt[datos.m+1, dim] == dest[dim] for dim in range(2))

    MODEL.update()

    objective = gp.quicksum(pgLit[g, i, t] + pgRit[g, i, t] for g, i, t in pgRit.keys()) + gp.quicksum(pgij[g, i, j] for g, i, j in pgij.keys()) + gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for g in T_index_prima for i in grafos[g-1].aristas) + gp.quicksum(3*dLRt[t] for t in dLRt.keys()) + gp.quicksum(3*dRLt[t] for t in dRLt.keys())

    # objective = gp.quicksum(dLRt[t] for t in dLRt.keys()) + gp.quicksum(dRLt[t] for t in dRLt.keys())


    MODEL.setObjective(objective, GRB.MINIMIZE)
    MODEL.Params.Threads = 6
    MODEL.Params.NumericFocus = 2
    # MODEL.Params.NonConvex = 2
    MODEL.Params.timeLimit = datos.tmax

    MODEL.update()

    MODEL.write('modelo.lp')
    MODEL.optimize()

    MODEL.write('solucion.json')
    MODEL.update()

    if MODEL.Status == 3:
        MODEL.computeIIS()
        MODEL.write('casa.ilp')
        result =  [np.nan, np.nan, np.nan, np.nan]
        if datos.grid:
            result.append('Grid')
        else:
            result.append('Delauney')

        if datos.alpha:
            result.append('Alpha-e')
        else:
            result.append('Alpha-g')
        result.append('Stages')
        result.append(grafo_data.mode)

        return result

    if MODEL.SolCount == 0:
        result =  [np.nan, np.nan, np.nan, np.nan]
        if datos.grid:
            result.append('Grid')
        else:
            result.append('Delauney')

        if datos.alpha:
            result.append('Alpha-e')
        else:
            result.append('Alpha-g')
        result.append('Stages')
        result.append(grafo_data.mode)

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

    if datos.alpha:
        result.append('Alpha-e')
    else:
        result.append('Alpha-g')

    result.append('Stages')
    result.append(grafo_data.mode)

    MODEL.write('solution_Stages.sol')

    vals_u = MODEL.getAttr('x', ugit)
    selected_u = gp.tuplelist((g, i, t) for g, i, t in vals_u.keys() if vals_u[g, i, t] > 0.5)
    # print('ugit = ' + str(selected_u))

    vals_z = MODEL.getAttr('x', zgij)
    selected_z = gp.tuplelist((g, i, j) for g, i, j in vals_z.keys() if vals_z[g, i, j] > 0.5)
    # print('zgij = ' + str(selected_z))


    vals_v = MODEL.getAttr('x', vgit)
    selected_v = gp.tuplelist((g, i, t) for g, i, t in vals_v.keys() if vals_v[g, i, t] > 0.5)
    # print('vgit = ' + str(selected_v))

    # print(dLRt)
    # print(dRLt)

    vals_muLet = MODEL.getAttr('x', muLet)
    selected_muLet = gp.tuplelist((e, t) for e, t in vals_muLet.keys() if vals_muLet[e, t] > 0.5)
    # print('muLet = ' + str(selected_muLet))

    vals_muRet = MODEL.getAttr('x', muRet)
    selected_muRet = gp.tuplelist((e, t) for e, t in vals_muRet.keys() if vals_muRet[e, t] > 0.5)
    # print('muRet = ' + str(selected_muRet))

    vals_zeeRLt = MODEL.getAttr('x', zeeRLt)
    selected_zeeRLt = gp.tuplelist((e1, e2, t) for e1, e2, t in vals_zeeRLt.keys() if vals_zeeRLt[e1, e2, t] > 0.5)
    # print('zeeRLt = ' + str(selected_zeeRLt))


    ind = 0
    path_C = []
    paths_D = []

    #path_C.append(origin)
    path_C.append([xLt[0, 0].X, xLt[0, 1].X])
    for t in T_index_prima:
        #    if ind < datos.m:
        path_C.append([xLt[t, 0].X, xLt[t, 1].X])
        if ind < datos.m:
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

    path_C.append([xLt[datos.m+1, 0].X, xLt[datos.m+1, 1].X])

    fig, ax = plt.subplots()
    plt.axis([0, 100, 0, 100])
    ax.set_aspect('equal')

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

    nx.draw(grafo.G, grafo.pos, node_size=100, node_color='black', alpha=0.3, edge_color='gray')
    nx.draw_networkx_labels(grafo.G, grafo.pos)

    for g in range(datos.m):
        grafo = grafos[g]
        nx.draw(grafo.G, grafo.pos, node_size=10,
                node_color='black', alpha=0.3, edge_color='gray')
        # nx.draw_networkx_labels(grafo.G, grafo.pos)

    plt.savefig('./images/NDST-grid-mode{0}.png'.format(grafo_data.mode))

    return result
    # plt.show()
