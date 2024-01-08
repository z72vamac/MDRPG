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

# np.random.seed(4)

# lista = [4, 4]
# nG = len(lista)
# datos = Data([], m=nG, r=1, grid = False, tmax=600, alpha = True,
#              init=True,
#              show=True,
#              seed=2)

# datos.generar_grid()



# datos.generar_grafos(lista)

def PDST(datos): #, vals_xL, vals_xR):

    orig = datos.orig
    dest = datos.dest

    grafos = datos.mostrar_datos()

    result = []

    nG = datos.m
    T_index = range(datos.m + 2)
    T_index_prima = range(1, datos.m+1)
    T_index_primaprima = range(datos.m+1)


    vD = datos.vD
    vC = datos.vC

    # Creamos el modelo8
    MODEL = gp.Model("PD-Stages")

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
    dRLt_index = T_index_primaprima

    dRLt = MODEL.addVars(dRLt_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dRLt')
    difRLt = MODEL.addVars(dRLt_index, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='difRLt')

    # Variable continua no negativa dLRt que indica la distancia que recorre el camión en la etapa t mientras el dron se mueve
    dLRt_index = T_index

    dLRt = MODEL.addVars(dLRt_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dLRt')
    difLRt = MODEL.addVars(dLRt_index, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='difLRt')

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

    # MODEL.addConstrs(ugit.sum('*', i, '*') == 1 for i in range(nG))
    # MODEL.addConstrs(vgit.sum('*', i, '*') == 1 for g in range(nG))

    # De todos los segmentos hay que salir y entrar menos de aquel que se toma como entrada al grafo y como salida del grafo
    MODEL.addConstrs(mugi[g, i] - ugit.sum(g, i, '*') == zgij.sum(g, '*', i) for g, i in rhogi.keys())
    MODEL.addConstrs(mugi[g, i] - vgit.sum(g, i, '*') == zgij.sum(g, i, '*') for g, i in rhogi.keys())

    MODEL.addConstrs(ugit.sum(g, '*', t) - vgit.sum(g, '*', t) == 0 for t in T_index_prima for g in T_index_prima)

    # MODEL.addConstrs(ugit[g, i, t] == vgi)
    # MODEL.addConstrs((ugit.sum(g, i, '*') + zgij.sum(g, '*', i))
    #                  - (vgit.sum(g, i, '*') + zgij.sum(g, i, '*')) == 0 for g in T_index_prima for i in grafos[g-1].aristas)

    MODEL.addConstrs(pgi[g, i] >= mugi[g, i] + alphagi[g, i] - 1 for g, i in rhogi.keys())
    MODEL.addConstrs(pgi[g, i] <= mugi[g, i] for g, i in rhogi.keys())
    MODEL.addConstrs(pgi[g, i] <= alphagi[g, i] for g, i in rhogi.keys())

    # MODEL.addConstr(ugit[0, 101, 0] == 0)
    # MODEL.addConstr(ugit[0, 101, 1] == 0)


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
    MODEL.addConstrs((difgLit[g, i, t, dim] >=   xLt[t, dim] - Rgi[g, i, dim]) for g, i, t, dim in difgLit.keys())
    MODEL.addConstrs((difgLit[g, i, t, dim] >= - xLt[t, dim] + Rgi[g, i, dim]) for g, i, t, dim in difgLit.keys())

    MODEL.addConstrs((difgLit[g, i, t, 0]*difgLit[g, i, t, 0] + difgLit[g, i, t, 1] * difgLit[g, i, t, 1] <= dgLit[g, i, t] * dgLit[g, i, t] for g, i, t in ugit.keys()), name = 'difgLit')

    SmallM = 0
    BigM = 10000
    #
    # BigM = 0
    # for g in T_index_prima:
    #     for v in grafos[g-1].V:
    #         BigM = max([np.linalg.norm(orig - v), BigM])
    #
    # BigM += 5
    #BigM = max([np.linalg.norm(orig-grafos[g].V) for g in range(nG)])
    MODEL.addConstrs((pgLit[g, i, t] >= SmallM * ugit[g, i, t]) for g, i, t in ugit.keys())
    MODEL.addConstrs((pgLit[g, i, t] >= dgLit[g, i, t] - BigM * (1 - ugit[g, i, t])) for g, i, t in ugit.keys())

    MODEL.addConstrs((difgij[g, i, j, dim] >=   Lgi[g, i, dim] - Rgi[g, j, dim]) for g, i, j, dim in difgij.keys())
    MODEL.addConstrs((difgij[g, i, j, dim] >= - Lgi[g, i, dim] + Rgi[g, j, dim]) for g, i, j, dim in difgij.keys())

    MODEL.addConstrs((difgij[g, i, j, 0]*difgij[g, i, j, 0] + difgij[g, i, j, 1] * difgij[g, i, j, 1] <= dgij[g, i, j] * dgij[g, i, j] for g, i, j in dgij.keys()), name = 'difgij')


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

    MODEL.addConstrs((difgRit[g, i, t, 0]*difgRit[g, i, t, 0] + difgRit[g, i, t, 1] * difgRit[g, i, t, 1] <= dgRit[g, i, t] * dgRit[g, i, t] for g, i, t in vgit.keys()), name = 'difgRit')


    #SmallM = 0
    #BigM = 10000
    MODEL.addConstrs((pgRit[g, i, t] >= SmallM * vgit[g, i, t]) for g, i, t in vgit.keys())
    MODEL.addConstrs((pgRit[g, i, t] >= dgRit[g, i, t] - BigM * (1 - vgit[g, i, t])) for g, i, t in vgit.keys())

    MODEL.addConstrs((difRLt[t, dim] >=   xRt[t, dim] - xLt[t + 1, dim] for t in dRLt.keys() for dim in range(2)), name = 'error')
    MODEL.addConstrs((difRLt[t, dim] >= - xRt[t, dim] + xLt[t + 1, dim] for t in dRLt.keys() for dim in range(2)), name = 'error2')
    MODEL.addConstrs((difRLt[t, 0]*difRLt[t, 0] + difRLt[t, 1] * difRLt[t, 1] <= dRLt[t] * dRLt[t] for t in dRLt.keys()), name = 'difRLt')

    MODEL.addConstrs((difLRt[t, dim] >=   xLt[t, dim] - xRt[t, dim]) for t, dim in difLRt.keys())
    MODEL.addConstrs((difLRt[t, dim] >= - xLt[t, dim] + xRt[t, dim]) for t, dim in difLRt.keys())
    MODEL.addConstrs((difLRt[t, 0]*difLRt[t, 0] + difLRt[t, 1] * difLRt[t, 1] <= dLRt[t] * dLRt[t] for t in dLRt.keys()), name = 'difLRt')


    # longitudes = []
    # for g in T_index_prima:
    #     longitudes.append(sum([grafos[g-1].A[i // 100 - 1, i % 100]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas]))

    BigM = 1000000

    MODEL.addConstrs((gp.quicksum(pgLit[g, i, t] for i in grafos[g-1].aristas) + pgij.sum(g, '*', '*') +  gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas) + gp.quicksum(pgRit[g, i, t] for i in grafos[g-1].aristas))/vD <= dLRt[t]/vC + BigM*(1- gp.quicksum(ugit[g, i, t] for i in grafos[g-1].aristas)) for t in T_index_prima for g in T_index_prima)

    #MODEL.addConstrs(zgij[g, i, j] <= ugit.sum(g, '*', '*') for g, i, j in zgij.keys())
    #MODEL.addConstrs(mugi[g, i] <= ugit.sum(g, i, '*') for g, i in mugi.keys())

    # MODEL.addConstrs((pgLit.sum('*', '*', t) +
    #                   pgij.sum(g, '*', '*') +
    #                   ugit.sum(g, '*', '*')*longitudes[g-1] +
    #                   pgRit.sum('*', '*', t))/vD <= dLRt[t]/vC for t in T_index_prima for g in T_index_prima)
    # MODEL.addConstrs((dLRt[t]/vD <= 50) for t in T_index_prima)
    # MODEL.addConstrs((pgLit[g, i, t]
    #                   + pgij.sum(g, '*', '*') + grafos[g-1].A[i // 100 - 1, i % 100]*grafos[g-1].longaristas[i // 100 - 1, i % 100]
    #                   + pgRit[g, i, t])/vD <= dLRt[t]/vC for g, i, t in pgLit.keys())

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

    # [0, 2, 1, 3, 4]
    # MODEL.addConstr(ugit[2, 102, 1] >= 0.5)
    # MODEL.addConstr(ugit[1, 101, 2] >= 0.5)
    # MODEL.addConstr(ugit[3, 102, 3] >= 0.5)
    #
    # MODEL.addConstr(vgit[2, 303, 1] >= 0.5)
    # MODEL.addConstr(vgit[1, 203, 2] >= 0.5)
    # MODEL.addConstr(vgit[3, 101, 3] >= 0.5)
    #
    # MODEL.addConstr(dLRt[1] <= 3.1490912899469254e+01 + 1e-5)
    # MODEL.addConstr(dLRt[2] <= 8.1903472383647316e+00 + 1e-5)
    # MODEL.addConstr(dLRt[3] <= 3.2352819808730416e+01 + 1e-5)
    #
    #
    # MODEL.addConstr(xLt[1, 0] == 5.0000000005633247e+01)
    # MODEL.addConstr(xLt[1, 1] == 4.9999999994606391e+01)
    # MODEL.addConstr(xLt[2, 0] == 45.0173764234826)
    # MODEL.addConstr(xLt[2, 1] == 7.3845420113314660e+01)
    # MODEL.addConstr(xLt[3, 0] == 4.5080138866746275e+01)
    # MODEL.addConstr(xLt[3, 1] == 8.0665242773352233e+01)
    #
    # MODEL.addConstr(xRt[1, 0] == 4.5017376430147451e+01)
    # MODEL.addConstr(xRt[1, 1] == 7.3845420106065134e+01)
    # MODEL.addConstr(xRt[2, 0] == 45.0801388601561)
    # MODEL.addConstr(xRt[2, 1] == 80.6652427666625)
    # MODEL.addConstr(xRt[3, 0] == 5.0000000002577934e+01)
    # MODEL.addConstr(xRt[3, 1] == 5.0000000007343687e+01)

    # Origen y destino
    MODEL.addConstrs(xLt[0, dim] == orig[dim] for dim in range(2))
    MODEL.addConstrs(xRt[0, dim] == orig[dim] for dim in range(2))

    MODEL.addConstrs(xLt[datos.m+1, dim] == dest[dim] for dim in range(2))
    MODEL.addConstrs(xRt[datos.m+1, dim] == dest[dim] for dim in range(2))

    # print(vals_xL)
    # for g in T_index_prima:
    #     MODEL.addConstrs(xLt[g, dim] == vals_xL[g][dim] for dim in range(2))
    #     MODEL.addConstrs(xRt[g, dim] == vals_xR[g][dim] for dim in range(2))

    MODEL.update()

    objective = gp.quicksum(pgLit[g, i, t] + pgRit[g, i, t] for g, i, t in pgRit.keys()) + gp.quicksum(pgij[g, i, j] for g, i, j in pgij.keys()) + gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for g in T_index_prima for i in grafos[g-1].aristas) + gp.quicksum(3*dLRt[t] for t in dLRt.keys()) + gp.quicksum(3*dRLt[t] for t in dRLt.keys())

    # objective = gp.quicksum(1*dLRt[g] for g in dLRt.keys()) + gp.quicksum(1*dRLt[g1] for g1 in dRLt.keys())

    # objective = gp.quicksum(dRLt[t] + dLRt[t] for t in T_index)

    MODEL.setObjective(objective, GRB.MINIMIZE)
    MODEL.Params.Threads = 8
    # MODEL.Params.NonConvex = 2
    MODEL.Params.timeLimit = datos.tmax

    MODEL.update()

    MODEL.write('AMDRPG-Stages.lp')
    MODEL.write('AMDRPG-Stages.mps')
    MODEL.optimize()


    # MODEL.update()

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

    MODEL.write('solution_Stages.sol')

    # vals_u = MODEL.getAttr('x', ugit)
    # selected_u = gp.tuplelist((g, i, t)
    #                           for g, i, t in vals_u.keys() if vals_u[g, i, t] > 0.5)
    # # print(selected_u)
    #
    # vals_z = MODEL.getAttr('x', zgij)
    # selected_z = gp.tuplelist((g, i, j)
    #                           for g, i, j in vals_z.keys() if vals_z[g, i, j] > 0.5)
    # # print(selected_z)
    #
    # vals_v = MODEL.getAttr('x', vgit)
    # selected_v = gp.tuplelist((g, i, t)
    #                           for g, i, t in vals_v.keys() if vals_v[g, i, t] > 0.5)
    # # print(selected_v)
    #
    # path = []
    # path.append(0)
    #
    # for t in T_index_prima:
    #     tripleta = selected_u.select('*', '*', t)
    #     path.append(tripleta[0][0])
    #
    # path.append(nG+1)
    # print(path)
    #
    # ind = 0
    # path_C = []
    # paths_D = []
    #
    # #path_C.append(orig)
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
    #         limite = sum([1 for g, i, j in selected_z if g == index_g])
    #         while count < limite:
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


    # ind = 0
    # path_C = []
    # paths_D = []
    #
    # #path_C.append(orig)
    # path_C.append([xLt[0, 0].X, xLt[0, 1].X])
    # for t in T_index_prima:
    #     #    if ind < nG:
    #     path_C.append([xLt[t, 0].X, xLt[t, 1].X])
    #     if ind < nG:
    #         path_D = []
    #         path_D.append([xLt[t, 0].X, xLt[t, 1].X])
    #         index_g = 0
    #         index_i = 0
    #         for g, i, ti in selected_u:
    #             if ti == t:
    #                 index_g = g
    #                 index_i = i
    #         count = 0
    #         while count < grafos[index_g-1].num_aristas-1:
    #             for g, i, j in selected_z:
    #                 if index_g == g and index_i == i:
    #                     path_D.append([Rgi[g, i, 0].X, Rgi[g, i, 1].X])
    #                     path_D.append([Lgi[g, i, 0].X, Lgi[g, i, 1].X])
    #                     count += 1
    #                     index_i = j
    #         path_D.append([Rgi[index_g, index_i, 0].X, Rgi[index_g, index_i, 1].X])
    #         path_D.append([Lgi[index_g, index_i, 0].X, Lgi[index_g, index_i, 1].X])
    #         ind += 1
    #         path_D.append([xRt[t, 0].X, xRt[t, 1].X])
    #     paths_D.append(path_D)
    #     path_C.append([xRt[t, 0].X, xRt[t, 1].X])

    # path_C.append([xLt[datos.m+1, 0].X, xLt[datos.m+1, 1].X])
    #
    # fig, ax = plt.subplots()
    # plt.axis([0, 100, 0, 100])
    #
    # for g, i in rhogi.keys():
    #     plt.plot(Rgi[g, i, 0].X, Rgi[g, i, 1].X, 'kD', markersize=1, color='cyan')
    #     plt.plot(Lgi[g, i, 0].X, Lgi[g, i, 1].X, 'kD', markersize=1, color='cyan')
    # #
    # # path_C = []
    # for t in T_index:
    #     # path_C.append([xLt[t, 0].X, xLt[t, 1].X])
    #     # path_C.append([xRt[t, 0].X, xRt[t, 1].X])
    #     plt.plot(xLt[t, 0].X, xLt[t, 1].X, 'ko', alpha = 0.3, markersize=10, color='green')
    #     #ax.annotate("L" + str(t), xy = (xLt[t, 0].X, xLt[t, 1].X))
    #     plt.plot(xRt[t, 0].X, xRt[t, 1].X, 'ko', markersize=5, color='blue')
    #     #ax.annotate("R" + str(t), xy = (xRt[t, 0].X, xRt[t, 1].X))
    #
    # ax.add_artist(Polygon(path_C, fill=False, animated=False,
    #               linestyle='-', alpha=1, color='blue'))
    #
    # for pathd in paths_D:
    #     ax.add_artist(Polygon(pathd, fill=False, closed=False, lw = 0.1,
    #                   animated=False, alpha=1, color='red'))
    # #
    # # ax.add_artist(Polygon(path_D, fill=False, animated=False,
    # #               linestyle='dotted', alpha=1, color='red'))
    #
    # for g in T_index_prima:
    #     grafo = grafos[g-1]
    #     centroide = np.mean(grafo.V, axis = 0)
    #     nx.draw(grafo.G, grafo.pos, node_size=30,
    #             node_color='black', alpha=1, edge_color='gray')
    #     ax.annotate(g, xy = (centroide[0], centroide[1]))
    #     nx.draw_networkx_labels(grafo.G, grafo.pos, font_color = 'red', font_size=5)
    #
    # plt.savefig('PDST-' + str(result[4]) +  '.png')
    #
    # plt.show()
    print(result)
    return result
    # plt.show()

# PDST(datos)
