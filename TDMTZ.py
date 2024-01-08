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
import auxiliar_functions as af

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

# np.random.seed(1)
#
# lista = [4, 4, 4, 4]
# nG = len(lista)
# datos = Data([], m=nG, r=1, grid = False, tmax=120, alpha = True,
#              init=True,
#              show=True,
#              seed=2)
#
# datos.generar_grid()
# #
# # # path = [0, 2, 3, 4, 1, 5]
# # # z = af.path2matrix(path)
# #
# datos.generar_grafos(lista)

def TDMTZ(datos, ciclo):

    origin= [50,50]
    dest = origin
    grafos = datos.mostrar_datos()
    # print(grafos[0].V)

    result = []

    nG = datos.m
    T_index = range(datos.m + 2)
    T_index_prima = range(1, datos.m+1)
    T_index_primaprima = range(datos.m+1)

    # ciclos = Data([], m = 1, r = 6, grid = False, tmax=120, alpha = True,
    #                 init = True,
    #                 show = True,
    #                 seed = 2)
    # ciclos.generar_ciclo()
    #
    # ciclo = ciclos.mostrar_datos()[0]

    vD = 2

    vC = 1
    # Creamos el modelo8
    MODEL = gp.Model("PD-Mtz")

    # Variables que modelan las distancias
    # Variable binaria ugi = 1 si en la etapa t entramos por el segmento sgi
    ugi_index = []

    for g in T_index_prima:
        for i in grafos[g-1].aristas:
            ugi_index.append((g, i))


    ugi = MODEL.addVars(ugi_index, vtype=GRB.BINARY, name='ugi')

    # Variable continua no negativa dgLi que indica la distancia desde el punto de lanzamiento hasta el segmento
    # sgi.
    dgLi_index = ugi_index

    dgLi = MODEL.addVars(dgLi_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dgLi')
    difgLi = MODEL.addVars(dgLi_index, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='difgLi')

    # Variable continua no negativa pgLi = ugi * dgLi
    pgLi_index = ugi_index

    pgLi = MODEL.addVars(pgLi_index, vtype=GRB.CONTINUOUS, lb=0.0, name='pgLi')


    # Variable binaria vgi = 1 si en la etapa t salimos por el segmento sgi
    vgi_index = ugi_index

    vgi = MODEL.addVars(vgi_index, vtype=GRB.BINARY, name='vgi')

    # Variable continua no negativa dgRi que indica la distancia desde el punto de salida del segmento sgi hasta el
    # punto de recogida del camion
    dgRi_index = ugi_index

    dgRi = MODEL.addVars(dgRi_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dgRi')
    difgRi = MODEL.addVars(dgRi_index, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='difgRi')


    # Variable continua no negativa pgRi = vgi * dgRi
    pgRi_index = ugi_index

    pgRi = MODEL.addVars(pgRi_index, vtype=GRB.CONTINUOUS, lb=0.0, name='pgRi')


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
    difgij = MODEL.addVars(
        dgij_index, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='dgij')

    # Variable continua no negativa pgij = zgij * dgij
    pgij_index = zgij_index

    pgij = MODEL.addVars(pgij_index, vtype=GRB.CONTINUOUS, lb=0.0, name='pgij')

    # Parametrizacion del ciclo
    muxLg = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'muxLg')
    muxRg = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'muxRg')

    # Variables binaria que es 1 si en la etapa t el punto xL se encuentra en el segmento j
    muLig = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.BINARY, name = 'muLig')

    # Variables binaria que es 1 si en la etapa t el punto xR se encuentra en el segmento j
    muRig = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.BINARY, name = 'muRig')

    minig = MODEL.addVars(ciclo.num_segmentos, T_index, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'minig')
    maxig = MODEL.addVars(ciclo.num_segmentos, T_index, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'maxig')

    entryig = MODEL.addVars(ciclo.num_segmentos, T_index, T_index, vtype = GRB.BINARY, name = 'entryig')

    # Distancia del punto de lanzamiento al punto de recogida
    dLR = MODEL.addVars(T_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dLR')

    zijLRg = MODEL.addVars(ciclo.num_segmentos, ciclo.num_segmentos, T_index, vtype = GRB.BINARY, lb = 0.0, name = 'zijLRg')
    dijLRg = MODEL.addVars(ciclo.num_segmentos, ciclo.num_segmentos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dijLRg')
    pijLRg = MODEL.addVars(ciclo.num_segmentos, ciclo.num_segmentos, T_index, vtype = GRB.CONTINUOUS, name = 'pijLRg')

    zijRLgg = MODEL.addVars(ciclo.num_segmentos, ciclo.num_segmentos, T_index, T_index, vtype = GRB.BINARY, name = 'zijRLgg')
    dijRLgg = MODEL.addVars(ciclo.num_segmentos, ciclo.num_segmentos, T_index, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dijRLgg')
    pijRLgg = MODEL.addVars(ciclo.num_segmentos, ciclo.num_segmentos, T_index, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'pijRLgg')

    # Variable que modela el producto de muxL con muLig
    pLig = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pLig')

    # Variable que modela el producto de muxR con muRig
    pRig = MODEL.addVars(ciclo.num_segmentos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pRig')


    # Variable binaria z que vale uno si se va del grafo g al grafo g'
    z_index = []

    for v in T_index:
        for w in T_index:
            if v != w:
                z_index.append((v, w))

    z = MODEL.addVars(z_index, vtype=GRB.BINARY, name='z')
    s = MODEL.addVars(T_index, vtype=GRB.CONTINUOUS, lb=0, name='s')

    dRL = MODEL.addVars(z_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dRL')
    pRL = MODEL.addVars(z_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'pRL')


    # Variables que modelan los puntos de entrada o recogida
    # xL: punto de salida del dron del camion en la etapa t
    xL_index = []

    for g in T_index:
        for dim in range(2):
            xL_index.append((g, dim))

    xL = MODEL.addVars(xL_index, vtype=GRB.CONTINUOUS, name='xL')

    # xR: punto de recogida del dron del camion en la etapa t
    xR_index = []

    for t in T_index:
        for dim in range(2):
            xR_index.append((t, dim))

    xR = MODEL.addVars(xR_index, vtype=GRB.CONTINUOUS, name='xR')

    # Rgi: punto de recogida del dron para el segmento sgi
    Rgi_index = []
    rhogi_index = []

    for g in T_index_prima:
        for i in grafos[g-1].aristas:
            rhogi_index.append((g, i))
            for dim in range(2):
                Rgi_index.append((g, i, dim))

    Rgi = MODEL.addVars(Rgi_index, vtype=GRB.CONTINUOUS, name='Rgi')
    rhogi = MODEL.addVars(rhogi_index, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='rhogi')

    # Lgi: punto de lanzamiento del dron del segmento sgi
    Lgi_index = Rgi_index
    landagi_index = rhogi_index

    Lgi = MODEL.addVars(Lgi_index, vtype=GRB.CONTINUOUS, name='Lgi')
    landagi = MODEL.addVars(landagi_index, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='landagi')

    # Variables difiliares para modelar el valor absoluto
    mingi = MODEL.addVars(rhogi_index, vtype=GRB.CONTINUOUS, lb=0.0, ub = 1.0, name='mingi')
    maxgi = MODEL.addVars(rhogi_index, vtype=GRB.CONTINUOUS, lb=0.0, ub = 1.0, name='maxgi')
    entrygi = MODEL.addVars(rhogi_index, vtype=GRB.BINARY, name='entrygi')
    mugi = MODEL.addVars(rhogi_index, vtype = GRB.BINARY, name = 'mugi')
    pgi = MODEL.addVars(rhogi_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pgi')
    alphagi = MODEL.addVars(rhogi_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'alphagi')

    MODEL.update()


    # Para cada grafo g, existe un segmento i y una etapa t donde hay que recoger al dron
    MODEL.addConstrs((ugi.sum(g, '*') == 1 for g in T_index_prima), name = 'entrag')
    MODEL.addConstrs((vgi.sum(g, '*') == 1 for g in T_index_prima), name = 'saleg')

    # MODEL.addConstrs(ugi.sum('*', i, '*') == 1 for i in range(nG))
    # MODEL.addConstrs(vgi.sum('*', i, '*') == 1 for g in range(nG))

    # De todos los segmentos hay que salir y entrar menos de aquel que se toma como entrada al grafo y como salida del grafo
    MODEL.addConstrs((mugi[g, i] - ugi[g, i] == zgij.sum(g, '*', i) for g, i, j in zgij.keys()), name = 'flujou')
    MODEL.addConstrs((mugi[g, i] - vgi[g, i] == zgij.sum(g, i, '*') for g, i, j in zgij.keys()), name = 'flujov')

    MODEL.addConstrs((pgi[g, i] >= mugi[g, i] + alphagi[g, i] - 1 for g, i in rhogi.keys()), name = 'pgi1')
    MODEL.addConstrs((pgi[g, i] <= mugi[g, i] for g, i in rhogi.keys()), name = 'pgi2')
    MODEL.addConstrs((pgi[g, i] <= alphagi[g, i] for g, i in rhogi.keys()), name = 'pgi3')

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
    MODEL.addConstrs((difgLi[g, i, dim] >=   xL[g, dim] - Rgi[g, i, dim]) for g, i, dim in difgLi.keys())
    MODEL.addConstrs((difgLi[g, i, dim] >= - xL[g, dim] + Rgi[g, i, dim]) for g, i, dim in difgLi.keys())

    MODEL.addConstrs((difgLi[g, i, 0]*difgLi[g, i, 0] + difgLi[g, i, 1] * difgLi[g, i, 1] <= dgLi[g, i] * dgLi[g, i] for g, i in ugi.keys()), name = 'u-conic')

    # SmallM = 10000
    # BigM = 0
    # for g in T_index_prima:
    #     for v in grafos[g-1].V:
    #         for w in ciclo.V:
    #             BigM = max([np.linalg.norm(w - v), BigM])
    #             SmallM = min([np.linalg.norm(w - v), SmallM])

    SmallM = 0
    BigM = 10000
    # MODEL.addConstr(ugi[1, 101] == 1)
    # MODEL.addConstr(vgi[1, 203] == 1)

    MODEL.addConstrs((pgLi[g, i] >= SmallM * ugi[g, i]) for g, i in ugi.keys())
    MODEL.addConstrs((pgLi[g, i] >= dgLi[g, i] - BigM * (1 - ugi[g, i])) for g, i in ugi.keys())

    MODEL.addConstrs((difgij[g, i, j, dim] >=   Lgi[g, i, dim] - Rgi[g, j, dim]) for g, i, j, dim in difgij.keys())
    MODEL.addConstrs((difgij[g, i, j, dim] >= - Lgi[g, i, dim] + Rgi[g, j, dim]) for g, i, j, dim in difgij.keys())

    MODEL.addConstrs((difgij[g, i, j, 0]*difgij[g, i, j, 0] + difgij[g, i, j, 1] * difgij[g, i, j, 1] <= dgij[g, i, j] * dgij[g, i, j] for g, i, j in dgij.keys()), name = 'zgij-conic')


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

    MODEL.addConstrs((difgRi[g, i, dim] >=   Lgi[g, i, dim] - xR[g, dim]) for g, i, dim in difgRi.keys())
    MODEL.addConstrs((difgRi[g, i, dim] >= - Lgi[g, i, dim] + xR[g, dim]) for g, i, dim in difgRi.keys())

    MODEL.addConstrs((difgRi[g, i, 0]*difgRi[g, i, 0] + difgRi[g, i, 1] * difgRi[g, i, 1] <= dgRi[g, i] * dgRi[g, i] for g, i in vgi.keys()), name = 'v-conic')


    # SmallM = 0
    #BigM = 10000
    MODEL.addConstrs((pgRi[g, i] >= SmallM * vgi[g, i]) for g, i in vgi.keys())
    MODEL.addConstrs((pgRi[g, i] >= dgRi[g, i] - BigM * (1 - vgi[g, i])) for g, i in vgi.keys())

    MODEL.addConstrs((pRL[g1, g2] >= SmallM * z[g1, g2] for g1, g2 in z_index))
    MODEL.addConstrs((pRL[g1, g2] >= dRL[g1, g2] - BigM * (1 - z[g1, g2]) for g1, g2 in z_index))

    # Restricciones para formar un tour
    MODEL.addConstr(gp.quicksum(z[v, 0] for v in T_index_prima) == 0)
    MODEL.addConstr(gp.quicksum(z[nG+1, w] for w in T_index_prima) == 0)
    MODEL.addConstrs(gp.quicksum(z[v , w] for w in T_index if w != v) == 1 for v in T_index)
    MODEL.addConstrs(gp.quicksum(z[w , v] for w in T_index if w != v) == 1 for v in T_index)

    # MODEL.addConstr(gp.quicksum(z[v, 0] for v in T_index if v != 0) == 0)
    # MODEL.addConstr(gp.quicksum(z[nG+1, w] for w in T_index_prima) == 0)
    # MODEL.addConstrs(gp.quicksum(z[v , w] for w in T_index_prima if w != v) == 1 for v in T_index_prima)
    # MODEL.addConstrs(gp.quicksum(z[w , v] for w in T_index_prima if w != v) == 1 for v in T_index_prima)
    # MODEL.addConstr(gp.quicksum(z[0, w] for w in T_index if w != 0) == 1)
    # MODEL.addConstr(gp.quicksum(z[v, nG+1] for v in T_index if v != nG+1) == 0)

    # Conectividad
    for v in T_index_prima:
        for w in T_index_prima:
            if v != w:
                MODEL.addConstr(len(T_index) - 1 >= (s[v] - s[w]) + len(T_index) * z[v, w])

    # for v in range(1, nG+1):
    #     MODEL.addConstr(s[v] - s[0] + (nG+1 - 2)*z[0, v] <= len(T_index) - 1)
    #
    # for v in range(1, nG+1):
    #     MODEL.addConstr(s[0] - s[v] + (nG+1 - 1)*z[v, 0] <= 0)

    # for v in range(1, nG+1):
    #     MODEL.addConstr(-z[0,v] - s[v] + (nG+1-3)*z[v,0] <= -2, name="LiftedLB(%s)"%v)
    #     MODEL.addConstr(-z[v,0] + s[v] + (nG+1-3)*z[0,v] <= nG+1-2, name="LiftedUB(%s)"%v)

    for v in T_index_prima:
        MODEL.addConstr(s[v] >= 1)
        MODEL.addConstr(s[v] <= len(T_index) - 1)

    MODEL.addConstr(s[0] == 0)
    MODEL.addConstr(s[nG + 1] == nG+1)


    MODEL.addConstrs((pgLi.sum(g, '*') + pgij.sum(g, '*', '*') +  gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas) + pgRi.sum(g, '*'))/vD <= dLR[g]/vC for g in T_index_prima)

    # MODEL.addConstrs(dLR[g] <= 150 for g in dLR.keys())
    # MODEL.addConstrs((pgLi.sum('*', '*', t) +
    #                   pgij.sum(g, '*', '*') +
    #                   ugi.sum(g, '*', '*')*longitudes[g-1] +
    #                   pgRi.sum('*', '*', t))/vD <= dLR[t]/vC for t in T_index_prima for g in T_index_prima)
    # MODEL.addConstrs((dLR[t]/vD <= 50) for t in T_index_prima)
    # MODEL.addConstrs((pgLi[g, i, t]
    #                   + pgij.sum(g, '*', '*') + grafos[g-1].A[i // 100 - 1, i % 100]*grafos[g-1].longaristas[i // 100 - 1, i % 100]
    #                   + pgRi[g, i, t])/vD <= dLR[t]/vC for g, i, t in pgLi.keys())

    # MODEL.addConstr(z[0, 2] + z[1, 3] + z[2, 1] + z[3, 4] == 4)

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

    for i, g1, g2 in entryig.keys():
        MODEL.addConstr(maxig[i, g1, g2] - minig[i, g1, g2] == muxRg[i, g1] - muxLg[i, g2])
        MODEL.addConstr(maxig[i, g1, g2] <= entryig[i, g1, g2])
        MODEL.addConstr(minig[i, g1, g2] <= 1 - entryig[i, g1, g2])

    for g in T_index:
        for dim in range(2):
            MODEL.addConstr(xL[g, dim] == gp.quicksum(muLig[j, g]*ciclo.V[j][dim] + pLig[j, g]*(ciclo.V[j+1][dim] - ciclo.V[j][dim]) for j in range(ciclo.num_segmentos)))
            MODEL.addConstr(xR[g, dim] == gp.quicksum(muRig[j, g]*ciclo.V[j][dim] + pRig[j, g]*(ciclo.V[j+1][dim] - ciclo.V[j][dim]) for j in range(ciclo.num_segmentos)))


        for i in range(ciclo.num_segmentos):
            MODEL.addConstr(pLig[i, g] >= muLig[i, g] + muxLg[i, g] - 1)
            MODEL.addConstr(pLig[i, g] <= muxLg[i, g])
            MODEL.addConstr(pLig[i, g] <= muLig[i, g])
            MODEL.addConstr(pRig[i, g] >= muRig[i, g] + muxRg[i, g] - 1)
            MODEL.addConstr(pRig[i, g] <= muxRg[i, g])
            MODEL.addConstr(pRig[i, g] <= muRig[i, g])

        for i in range(ciclo.num_segmentos):
            for j in range(ciclo.num_segmentos):
                MODEL.addConstr(zijLRg[i, j, g] <= muLig[i, g])
                MODEL.addConstr(zijLRg[i, j, g] <= muRig[j, g])
                MODEL.addConstr(zijLRg[i, j, g] >= muLig[i, g] + muRig[j, g] - 1)

                if i == j:
                    MODEL.addConstr(dijLRg[i, j, g] == (maxig[i, g, g] + minig[i, g, g])*ciclo.longitudes[i])
                if i < j:
                    MODEL.addConstr(dijLRg[i, j, g] == (1 - muxLg[i, g])*ciclo.longitudes[i] + sum([ciclo.longitudes[k] for k in range(i+1, j)]) + muxRg[j, g]*ciclo.longitudes[j])
                if i > j:
                    MODEL.addConstr(dijLRg[i, j, g] == muxLg[i, g]*ciclo.longitudes[i] + sum([ciclo.longitudes[k] for k in range(j+1, i)]) + (1-muxRg[j, g])*ciclo.longitudes[j])

                MODEL.addConstr(pijLRg[i, j, g] <= ciclo.longitud* zijLRg[i, j, g])
                MODEL.addConstr(pijLRg[i, j, g] <= dijLRg[i, j, g])
                MODEL.addConstr(pijLRg[i, j, g] >= dijLRg[i, j, g] - ciclo.longitud*(1 - zijLRg[i, j, g]))

        MODEL.addConstr(dLR[g] == gp.quicksum(pijLRg[j, k, g] for j in range(ciclo.num_segmentos) for k in range(ciclo.num_segmentos)))

    for g1, g2 in z.keys():
        for i in range(ciclo.num_segmentos):
            for j in range(ciclo.num_segmentos):
                MODEL.addConstr(zijRLgg[i, j, g1, g2] <= muRig[i, g1])
                MODEL.addConstr(zijRLgg[i, j, g1, g2] <= muLig[j, g2])
                MODEL.addConstr(zijRLgg[i, j, g1, g2] >= muRig[i, g1] + muLig[j, g2] - 1)

                if i == j:
                    MODEL.addConstr(dijRLgg[i, j, g1, g2] == (maxig[j, g1, g2] + minig[j, g1, g2])*ciclo.longitudes[i])
                if i < j:
                    MODEL.addConstr(dijRLgg[i, j, g1, g2] == (1 - muxRg[i, g1])*ciclo.longitudes[i] + sum([ciclo.longitudes[k] for k in range(i+1, j)]) + muxLg[j, g2]*ciclo.longitudes[j])
                if i > j:
                    MODEL.addConstr(dijRLgg[i, j, g1, g2] == muxRg[i, g1]*ciclo.longitudes[i] + sum([ciclo.longitudes[k] for k in range(j+1, i)]) + (1-muxLg[j, g2])*ciclo.longitudes[j])

                MODEL.addConstr(pijRLgg[i, j, g1, g2] <= ciclo.longitud* zijRLgg[i, j, g1, g2])
                MODEL.addConstr(pijRLgg[i, j, g1, g2] <= dijRLgg[i, j, g1, g2])
                MODEL.addConstr(pijRLgg[i, j, g1, g2] >= dijRLgg[i, j, g1, g2] - ciclo.longitud*(1 - zijRLgg[i, j, g1, g2]))

        MODEL.addConstr(dRL[g1, g2] == gp.quicksum(pijRLgg[j, k, g1, g2] for j in range(ciclo.num_segmentos) for k in range(ciclo.num_segmentos)))

    MODEL.addConstrs(muLig.sum('*', g) == 1 for g in T_index)
    MODEL.addConstrs(muRig.sum('*', g) == 1 for g in T_index)


    # Origen y destino
    MODEL.addConstrs(xL[0, dim] == ciclo.V[0][dim] for dim in range(2))
    MODEL.addConstrs(xR[0, dim] == ciclo.V[0][dim] for dim in range(2))

    MODEL.addConstrs(xL[datos.m+1, dim] == ciclo.V[0][dim] for dim in range(2))
    MODEL.addConstrs(xR[datos.m+1, dim] == ciclo.V[0][dim] for dim in range(2))

    MODEL.update()

    objective = gp.quicksum(pgLi[g, i] + pgRi[g, i] for g, i in pgRi.keys()) + gp.quicksum(pgij[g, i, j] for g, i, j in pgij.keys()) + gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for g in T_index_prima for i in grafos[g-1].aristas) + gp.quicksum(3*dLR[g] for g in dLR.keys()) + gp.quicksum(3*pRL[g1, g2] for g1, g2 in dRL.keys())

    # objective = gp.quicksum(1*dLR[g] for g in dLR.keys()) + gp.quicksum(1*pRL[g1, g2] for g1, g2 in dRL.keys()) # + gp.quicksum(1e5*holg[g] for g in holg.keys())

    # objective = gp.quicksum(dRL[t] + dLR[t] for t in T_index)

    MODEL.setObjective(objective, GRB.MINIMIZE)
    MODEL.Params.Threads = 6
    # MODEL.Params.NonConvex = 2
    MODEL.Params.timeLimit = datos.tmax
    # MODEL.Params.FeasibilityTol = 1e-2

    MODEL.update()
    # MODEL.computeIIS()
    # MODEL.write('casa.ilp')

    MODEL.write('AMDRPG-MTZ.lp')
    MODEL.write('AMDRPG-MTZ.mps')
    MODEL.optimize()

    # MODEL.write('solucion_TDMTZ.sol')


    # MODEL.update()

    if MODEL.Status == 3:
        MODEL.computeIIS()
        MODEL.write('casa.ilp')
        result =  [np.nan, np.nan, np.nan, np.nan]
        if datos.grid:
            result.append('Grid')
        else:
            result.append('Delauney')

        result.append('MTZ')

        return result

    if MODEL.SolCount == 0:
        result =  [np.nan, np.nan, np.nan, np.nan]
        if datos.grid:
            result.append('Grid')
        else:
            result.append('Delauney')

        result.append('MTZ')

        return result


    result.append(MODEL.getAttr('MIPGap'))
    result.append(MODEL.Runtime)
    result.append(MODEL.getAttr('NodeCount'))
    result.append(MODEL.ObjVal)

    if datos.grid:
        result.append('Grid')
    else:
        result.append('Delauney')

    result.append('MTZ')


    # vals_u = MODEL.getAttr('x', ugi)
    # selected_u = gp.tuplelist((g, i)
    #                           for g, i in vals_u.keys() if vals_u[g, i] > 0.5)
    # # print(selected_u)
    #
    # vals_zgij = MODEL.getAttr('x', zgij)
    # selected_zgij = gp.tuplelist((g, i, j)
    #                           for g, i, j in vals_zgij.keys() if vals_zgij[g, i, j] > 0.5)
    # # print(selected_zgij)
    #
    # vals_v = MODEL.getAttr('x', vgi)
    # selected_v = gp.tuplelist((g, i)
    #                           for g, i in vals_v.keys() if vals_v[g, i] > 0.5)
    # # print(selected_v)
    #
    # valsz = MODEL.getAttr('x', z)
    #
    # selected_z = gp.tuplelist(e for e in valsz if valsz[e] > 0)
    # # print(selected_z)
    # path = af.subtour(selected_z)
    # print(path)
    #
    # ind = 0
    # path_C = []
    # paths_D = []
    #
    # for p in path:
    #     path_C.append([xL[p, 0].X, xL[p, 1].X])
    #     path_C.append([xR[p, 0].X, xR[p, 1].X])
    #
    #
    # for p in path[1:]:
    #     #    if ind < nG:
    #     if ind < nG:
    #         path_D = []
    #         path_D.append([xL[p, 0].X, xL[p, 1].X])
    #         index_g = 0
    #         index_i = 0
    #         for g, i in selected_u:
    #             if g == p:
    #                 index_g = g
    #                 index_i = i
    #
    #         count = 0
    #         path_D.append([Rgi[index_g, index_i, 0].X, Rgi[index_g, index_i, 1].X])
    #         path_D.append([Lgi[index_g, index_i, 0].X, Lgi[index_g, index_i, 1].X])
    #         limite = sum([1 for g, i, j in selected_zgij if g == index_g])
    #         while count < limite:
    #             for g, i, j in selected_zgij:
    #                 if index_g == g and index_i == i:
    #                     count += 1
    #                     index_i = j
    #                     path_D.append([Rgi[index_g, index_i, 0].X, Rgi[index_g, index_i, 1].X])
    #                     path_D.append([Lgi[index_g, index_i, 0].X, Lgi[index_g, index_i, 1].X])
    #
    #         ind += 1
    #         path_D.append([xR[p, 0].X, xR[p, 1].X])
    #     paths_D.append(path_D)
    #
    # # path_C.append([xL[nG+1, 0].X, xL[nG+1, 1].X])
    # fig, ax = plt.subplots()
    # plt.axis([0, 100, 0, 100])
    # ax.set_aspect('equal')
    #
    # for g, i in rhogi.keys():
    #     first = i // 100 - 1
    #     second = i % 100
    #     if mugi[g, i].X > 0.5:
    #         plt.plot(Rgi[g, i, 0].X, Rgi[g, i, 1].X, 'kD', markersize=1, color='red')
    #         # ax.annotate("$R_" + str(g) + "^{" + str((first, second)) + "}$", xy = (Rgi[g, i, 0].X+1, Rgi[g, i, 1].X+1))
    #         plt.plot(Lgi[g, i, 0].X, Lgi[g, i, 1].X, 'kD', markersize=1, color='red')
    #         # ax.annotate("$L_" + str(g) + "^{" + str((first, second)) + "}$", xy = (Lgi[g, i, 0].X-3, Lgi[g, i, 1].X-3))
    #
    # #
    # for p, i in zip(path, range(len(path))):
    #     # path_C.append([xL[t, 0].X, xL[t, 1].X])
    #     # path_C.append([xR[t, 0].X, xR[t, 1].X])
    #     plt.plot(xL[p, 0].X, xL[p, 1].X, 'ko', alpha = 0.3, markersize=10, color='green')
    #     #ax.annotate("L" + str(i), xy = (xL[p, 0].X+0.5, xL[p, 1].X+0.5))
    #     plt.plot(xR[p, 0].X, xR[p, 1].X, 'ko', markersize=5, color='blue')
    #     #ax.annotate("R" + str(i), xy = (xR[p, 0].X-1.5, xR[p, 1].X-1.5))
    #
    #
    # ax.add_artist(Polygon(path_C, fill=False, closed = False, animated=False,
    #               linestyle='-', alpha=1, color='blue'))
    #
    # for pathd in paths_D:
    #     ax.add_artist(Polygon(pathd, fill=False, closed=False, lw = 0.1,
    #                   animated=False, alpha=1, color='red'))
    #
    # ax.add_artist(ciclo.artist)
    #
    # # ax.add_artist(Polygon(path_D, fill=False, animated=False,
    # #               linestyle='dotted', alpha=1, color='red'))
    #
    # for g in range(nG):
    #     grafo = grafos[g]
    #     nx.draw(grafo.G, grafo.pos, node_size=40,
    #             node_color='black', alpha=1, edge_color='black')
    #
    # plt.savefig('PDMTZ-' + str(result[4]) +  '.png')
    #
    # plt.show()

    print(result)
    return result

# PDMTZ(datos)
