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
from MTZ import *

# np.random.seed(2)
# orig = [50, 50]
# datos.dest = orig
#
# nG = 20
# datos = Data([], m=nG, r=3, modo=4, tmax=120, alpha = True,
#              init=True,
#              show=True,
#              seed=2)
# datos.generar_grafos()
# grafos = datos.mostrar_datos()
#
# T_index = range(datos.m + 2)
# T_index_prima = range(1, datos.m+1)
# T_index_primaprima = range(datos.m+1)
#
# vD = 2
#
# vC = 1

def NDMTZ_heuristic(datos, grafo_data):

    grafos = datos.mostrar_datos()

    T_index = range(datos.m + 2)
    T_index_prima = range(1, datos.m+1)
    T_index_primaprima = range(datos.m+1)

    orig = grafo_data.orig
    dest = grafo_data.dest

    first_time = time.time()
    results = []
    centroides = {}

    for g in T_index_prima:
        centroides[g] = np.mean(grafos[g-1].V, axis = 0)

    centros = []
    centros.append(orig)
    for g in centroides.values():
        centros.append(g)
    centros.append(datos.dest)

    elipses = []

    radio = 2

    for c in centros:
        P = np.identity(2)
        q = -2*np.array(c)
        r = c[0]**2 + c[1]**2 - radio**2
        elipse = Elipse(P, q, r)
        elipses.append(elipse)

    elipse = Data(elipses, m = 6, grid = True, alpha = datos.alpha, tmax = 1200, seed = 2)

    path_1, path_P, obj = MTZ(elipse)

    if datos.show:
        fig, ax = plt.subplots()
        plt.axis([0, 100, 0, 100])
        ax.set_aspect('equal')

        for e in elipse.data:
            ax.add_artist(e.artist)

        for g in T_index_prima:
            grafo = grafos[g-1]
            centroide = np.mean(grafo.V, axis = 0)
            nx.draw(grafo.G, grafo.pos, node_size=40,
                    node_color='black', alpha=1, edge_color='black')

        for c in centros:
            if c[0] == orig[0] and c[1] == orig[1]:
                ax.plot(c[0], c[1], 'o', color = 'orange', markersize = 3)
            else:
                ax.plot(c[0], c[1], 'o', color = 'red', markersize = 3)

        for p, i in zip(path_P, range(len(path_P)-1)):
            if i == 0:
                ax.annotate(str(i), xy = (p[0]-3.5, p[1]-1))
            elif i == 2:
                ax.annotate(str(i), xy = (p[0], p[1]-8.5))
            elif i == 3:
                ax.annotate(str(i), xy = (p[0]-8.5, p[1]))
            else:
                ax.annotate(str(i), xy = (p[0]-0.5, p[1]-8.5))


        ax.add_artist(Polygon(path_P, fill=False, animated=False, linestyle='-', alpha=1, color='blue'))

        plt.savefig('./imagenes/step0.png', dpi = 200)
        # plt.show()
    #
    z = af.path2matrix(path_1)
    path = path_1

    # fixed_data = [z]
    iter = 0

    xL_dict = {}
    xR_dict = {}

    xL_dict[0] = grafo_data.orig
    xR_dict[0] = grafo_data.orig


    for i in range(2, len(path)):
        g0 = path[i-2]
        g1 = path[i-1]
        g2 = path[i]
        print('Resolviendo para el grafo ' + str(g1))
        datos_partial = Data(data = datos.data, alpha = datos.alpha, m = len(datos.data), init = False)
        path_1 = [0, g1, path[-1]]
        print(path_1)
        # xL, xR, obj = PDMTZ_aux(datos_partial)
        xL, xR, obj = NDMTZ_aux(datos_partial, grafo_data, g1, xR_dict[g0])
        xL_dict[g1] = [xL[(g1, 0)], xL[(g1, 1)]]
        xR_dict[g1] = [xR[(g1, 0)], xR[(g1, 1)]]

    # xL_dict[g2] = [xL[(2, 0)], xL[(2, 1)]]
    # xR_dict[g2] = [xR[(2, 0)], xR[(2, 1)]]

    xL_dict[path[-1]] = grafo_data.dest
    xR_dict[path[-1]] = grafo_data.dest

    # [0, 1, 5, 4, 3, 2, 6]

    # print(xL_dict)
    #
    # for g in T_index:
    #     xL_dict[g] = [xL[(g, 0)], xL[(g, 1)]]
    #     xR_dict[g] = [xR[(g, 0)], xR[(g, 1)]]
    #
    fixed_data = [xL_dict, xR_dict]
    xL, xR, ugi, vgi, obj, path_2 = NDMTZ_aux2(datos, grafo_data, fixed_data, 0)

    obj_best = obj
    print(path_1)

    path_app = path_1.copy()
    path_app.reverse()

    print(path_app)
    iter = 1
    while not(all([i == j for i, j in zip(path_1, path_2)]) or all([i == j for i, j in zip(path_app, path_2)])) and obj_best < obj:
        path = path_2

        # fixed_data = [z]
        xL_dict = {}
        xR_dict = {}

        xL_dict[0] = datos.orig
        xR_dict[0] = datos.orig


        for i in range(2, len(path)):
            g0 = path[i-2]
            g1 = path[i-1]
            g2 = path[i]
            print('Resolviendo para el grafo ' + str(g1))
            datos_partial = Data(data = datos.data, alpha = datos.alpha, m = len(datos.data), init = False)
            path_1 = [0, g1, path[-1]]
            print(path_1)
            # xL, xR, obj = PDMTZ_aux(datos_partial)
            xL, xR, obj = NDMTZ_aux(datos_partial, g1, xR_dict[g0])
            xL_dict[g1] = [xL[(g1, 0)], xL[(g1, 1)]]
            xR_dict[g1] = [xR[(g1, 0)], xR[(g1, 1)]]

        # xL_dict[g2] = [xL[(2, 0)], xL[(2, 1)]]
        # xR_dict[g2] = [xR[(2, 0)], xR[(2, 1)]]

        xL_dict[path[-1]] = datos.dest
        xR_dict[path[-1]] = datos.dest

        # [0, 1, 5, 4, 3, 2, 6]

        # print(xL_dict)
        #
        # for g in T_index:
        #     xL_dict[g] = [xL[(g, 0)], xL[(g, 1)]]
        #     xR_dict[g] = [xR[(g, 0)], xR[(g, 1)]]
        #
        fixed_data = [xL_dict, xR_dict]
        xL, xR, ugi, vgi, obj, path_2 = NDMTZ_aux2(datos, grafo_data, fixed_data, iter)

        path_app = path.copy()
        path_app.reverse()
        # print(path_1)
        # print(path_2)
        # print(path_app)
        obj_best = max(obj_best, obj)
        iter += 1

    second_time = time.time()
    runtime = second_time - first_time

    if datos.init:
        # z = af.path2matrix(path)
        #
        # for g in T_index:
        #     xL_dict[g] = [xL[(g, 0)], xL[(g, 1)]]
        #     xR_dict[g] = [xR[(g, 0)], xR[(g, 1)]]

        results.append(obj)
        results.append(runtime)
        results.append(iter-1)
        return z, xL_dict, xR_dict, results
    else:

        results.append(obj)
        results.append(runtime)
        results.append(iter-1)

        # if datos.grid:
        #     results.append('Grid')
        # else:
        #     results.append('Delauney')

        return results


def NDMTZ_aux(datos, grafo_data, index_g, xR_dict):

    orig = grafo_data.orig
    dest = grafo_data.dest
    grafos = datos.mostrar_datos()
    # print(grafos[0].V)

    result = []

    nG = datos.m
    T_index = [0, index_g, nG+1]
    T_index_prima = [index_g]

    vD = datos.vD
    vC = datos.vC

    grafo = grafo_data.data[0]


    # Creamos el modelo NDMTZ
    MODEL = gp.Model("NDMTZ_aux")

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
    dLR = MODEL.addVars(T_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dLR')

    # ParamuLegrizacion del grafo
    muxLg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'muxLg')
    mineLRg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'mineLRg')
    maxeLRg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'maxeLRg')

    mineRLgg = MODEL.addVars(grafo.aristas, z_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'mineRLgg')
    maxeRLgg = MODEL.addVars(grafo.aristas, z_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'maxeRLgg')

    #aet = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'aet')
    entryeLRg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.BINARY, name = 'entryeLRg')
    entryeRLgg = MODEL.addVars(grafo.aristas, z_index, vtype = GRB.BINARY, name = 'entryeRLgg')

    # Variables binaria que es 1 si en la etapa t el punto xL se encuentra en el segmento e
    muLeg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.BINARY, name = 'muLeg')

    pLeg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'pLeg')

    # Variable binaria que toma el valor 1 cuando el camion entra al vertice i
    biLRg = MODEL.addVars(grafo.num_puntos, T_index, vtype = GRB.BINARY, name = 'biLRg')
    biRLg = MODEL.addVars(grafo.num_puntos, T_index, vtype = GRB.BINARY, name = 'biRLg')

    # Variable que modela el producto de muxLg con biLRg
    pxLbiLRg = MODEL.addVars(grafo.aristas, grafo.num_puntos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pxLbiLRg')
    pxLbiRLg = MODEL.addVars(grafo.aristas, grafo.num_puntos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pxLbiRLg')


    # Variable binaria que toma el valor 1 cuando el camion atraviesa la arista e
    qeLRg = MODEL.addVars(grafo.num_puntos, grafo.num_puntos, T_index, vtype = GRB.INTEGER, name = 'qeLRg')
    qeRLgg = MODEL.addVars(grafo.num_puntos, grafo.num_puntos, z_index, vtype = GRB.INTEGER, name = 'qeRLgg')

    muxRg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'muxRg')

    # Variables binaria que es 1 si en la etapa t el punto xR se encuentra en el segmento e
    muReg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.BINARY, name = 'muReg')

    pReg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'pReg')

    # Variable binaria que toma el valor 1 cuando el camion entra al vertice i
    ciLRg = MODEL.addVars(grafo.num_puntos, T_index, vtype = GRB.BINARY, name = 'ciLRg')
    ciRLg = MODEL.addVars(grafo.num_puntos, T_index, vtype = GRB.BINARY, name = 'ciRLg')

    # Variable que modela el producto de muxLg con muLjt
    pxRciLRg = MODEL.addVars(grafo.aristas, grafo.num_puntos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pxRciLRg')
    pxRciRLg = MODEL.addVars(grafo.aristas, grafo.num_puntos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pxRciRLg')


    zeeLRg = MODEL.addVars(grafo.aristas, grafo.aristas, T_index, vtype = GRB.BINARY, name = 'zeeLRg')
    deeLRg = MODEL.addVars(grafo.aristas, grafo.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'deeLRg')
    peeLRg = MODEL.addVars(grafo.aristas, grafo.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'peeLRg')

    zeeRLgg = MODEL.addVars(grafo.aristas, grafo.aristas, z_index, vtype = GRB.BINARY, name = 'zeeRLgg')
    deeRLgg = MODEL.addVars(grafo.aristas, grafo.aristas, z_index, vtype = GRB.CONTINUOUS, name = 'deeRLgg')
    peeRLgg = MODEL.addVars(grafo.aristas, grafo.aristas, z_index, vtype = GRB.CONTINUOUS, name = 'peeRLgg')





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

    # MODEL.addConstrs(ugi.sum('*', i, '*') == 1 for i in range(datos.m))
    # MODEL.addConstrs(vgi.sum('*', i, '*') == 1 for g in range(datos.m))

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

    # for g in range(datos.m):
    #     MODEL.addConstr(sgi[g, grafos[g].aristas[0]] == 0)

    for g in T_index_prima:
        for i in grafos[g-1].aristas[0:]:
            MODEL.addConstr(sgi[g, i] >= 0)
            MODEL.addConstr(sgi[g, i] <= grafos[g-1].num_aristas - 1)


    # Restricciones de distancias y producto
    MODEL.addConstrs((difgLi[g, i, dim] >=   xL[g, dim] - Rgi[g, i, dim]) for g, i, dim in difgLi.keys())
    MODEL.addConstrs((difgLi[g, i, dim] >= - xL[g, dim] + Rgi[g, i, dim]) for g, i, dim in difgLi.keys())

    MODEL.addConstrs((difgLi[g, i, 0]*difgLi[g, i, 0] + difgLi[g, i, 1] * difgLi[g, i, 1] <= dgLi[g, i] * dgLi[g, i] for g, i in ugi.keys()), name = 'u-conic')

    for g in T_index_prima:
        SmallM = 0
        BigM = 10000

        # BigM = max(max([np.linalg.norm(v - w) for v in grafo.V for w in grafos[g-1].V]), BigM)
        # SmallM = min(min([np.linalg.norm(v - w) for v in grafo.V for w in grafos[g-1].V]), SmallM)

        MODEL.addConstrs((pgLi[g, i] >= SmallM * ugi[g, i]) for i in grafos[g-1].aristas)
        MODEL.addConstrs((pgLi[g, i] >= dgLi[g, i] - BigM * (1 - ugi[g, i])) for i in grafos[g-1].aristas)

        MODEL.addConstrs((pgRi[g, i] >= SmallM * vgi[g, i]) for i in grafos[g-1].aristas)
        MODEL.addConstrs((pgRi[g, i] >= dgRi[g, i] - BigM * (1 - vgi[g, i])) for i in grafos[g-1].aristas)


    MODEL.addConstrs((difgij[g, i, j, dim] >=   Lgi[g, i, dim] - Rgi[g, j, dim]) for g, i, j, dim in difgij.keys())
    MODEL.addConstrs((difgij[g, i, j, dim] >= - Lgi[g, i, dim] + Rgi[g, j, dim]) for g, i, j, dim in difgij.keys())

    MODEL.addConstrs((difgij[g, i, j, 0]*difgij[g, i, j, 0] + difgij[g, i, j, 1] * difgij[g, i, j, 1] <= dgij[g, i, j] * dgij[g, i, j] for g, i, j in dgij.keys()), name = 'zgij-conic')


    for g, i, j in zgij.keys():
        first_i = i // 100 - 1
        second_i = i % 100
        first_j = j // 100 - 1
        second_j = j % 100

        segm_i = Poligonal(np.array([[grafos[g-1].V[first_i, 0], grafos[g-1].V[first_i, 1]], [grafos[g-1].V[second_i, 0], grafos[g-1].V[second_i, 1]]]), grafos[g-1].A[first_i, second_i])
        segm_j = Poligonal(np.array([[grafos[g-1].V[first_j, 0], grafos[g-1].V[first_j, 1]], [grafos[g-1].V[second_j, 0], grafos[g-1].V[second_j, 1]]]), grafos[g-1].A[first_j, second_j])

        BigM_local = eM.estima_BigM_local(segm_i, segm_j)
        SmallM_local = eM.estima_SmallM_local(segm_i, segm_j)
        MODEL.addConstr((pgij[g, i, j] >= SmallM_local * zgij[g, i, j]))
        MODEL.addConstr((pgij[g, i, j] >= dgij[g, i, j] - BigM_local * (1 - zgij[g, i, j])))

    MODEL.addConstrs((difgRi[g, i, dim] >=   Lgi[g, i, dim] - xR[g, dim]) for g, i, dim in difgRi.keys())
    MODEL.addConstrs((difgRi[g, i, dim] >= - Lgi[g, i, dim] + xR[g, dim]) for g, i, dim in difgRi.keys())

    MODEL.addConstrs((difgRi[g, i, 0]*difgRi[g, i, 0] + difgRi[g, i, 1] * difgRi[g, i, 1] <= dgRi[g, i] * dgRi[g, i] for g, i in vgi.keys()), name = 'v-conic')

    MODEL.addConstrs((pRL[g1, g2] >= SmallM * z[g1, g2] for g1, g2 in z_index))
    MODEL.addConstrs((pRL[g1, g2] >= dRL[g1, g2] - BigM * (1 - z[g1, g2]) for g1, g2 in z_index))

    # Restricciones para formar un tour
    MODEL.addConstr(z.sum('*', 0) == 0)
    MODEL.addConstr(z.sum(datos.m+1, '*') == 0)
    MODEL.addConstr(z.sum(0, '*') - z[0, datos.m+1] == 1)
    MODEL.addConstr(z.sum('*', datos.m+1) - z[0, datos.m+1] == 1)
    MODEL.addConstrs(z.sum(v, '*') == 1 for v in T_index_prima)
    MODEL.addConstrs(z.sum('*', v) == 1 for v in T_index_prima)


    # MODEL.addConstrs(gp.quicksum(z[v , w] for w in T_index_prima) == 1 for v in T_index)
    # MODEL.addConstrs(gp.quicksum(z[w , v] for w in T_index_prima if w != v) == 1 for v in T_index)

    # Conectividad
    for v in T_index_prima:
        for w in T_index_prima:
            if v != w:
                MODEL.addConstr(len(T_index) - 1 >= (s[v] - s[w]) + len(T_index) * z[v, w])

    # for v in range(1, datos.m+1):
    #     MODEL.addConstr(s[v] - s[0] + (datos.m+1 - 2)*z[0, v] <= len(T_index) - 1)
    #
    # for v in range(1, datos.m+1):
    #     MODEL.addConstr(s[0] - s[v] + (datos.m+1 - 1)*z[v, 0] <= 0)

    # for v in range(1, datos.m+1):
    #     MODEL.addConstr(-z[0,v] - s[v] + (datos.m+1-3)*z[v,0] <= -2, name="LiftedLB(%s)"%v)
    #     MODEL.addConstr(-z[v,0] + s[v] + (datos.m+1-3)*z[0,v] <= datos.m+1-2, name="LiftedUB(%s)"%v)

    for v in T_index_prima:
        MODEL.addConstr(s[v] >= 1)
        MODEL.addConstr(s[v] <= len(T_index) - 1)

    MODEL.addConstr(s[0] == 0)
    MODEL.addConstr(s[datos.m + 1] == datos.m+1)


    MODEL.addConstrs((pgLi.sum(g, '*') + pgij.sum(g, '*', '*') +  gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas) + pgRi.sum(g, '*'))/vD <= dLR[g]/vC for g in T_index_prima)

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

    for e, g in pLeg.keys():
        MODEL.addConstr(pLeg[e, g] >= muLeg[e, g] + muxLg[e, g] - 1)
        MODEL.addConstr(pLeg[e, g] <= muLeg[e, g])
        MODEL.addConstr(pLeg[e, g] <= muxLg[e, g])

        MODEL.addConstr(pReg[e, g] >= muReg[e, g] + muxRg[e, g] - 1)
        MODEL.addConstr(pReg[e, g] <= muReg[e, g])
        MODEL.addConstr(pReg[e, g] <= muxRg[e, g])

    for g, dim in xL.keys():
        MODEL.addConstr(xL[g, dim] == gp.quicksum(muLeg[e, g]*grafo.V[e // 100 - 1][dim] + pLeg[e, g]*(grafo.V[e % 100][dim] - grafo.V[e // 100 - 1][dim]) for e in grafo.aristas))
        MODEL.addConstr(xR[g, dim] == gp.quicksum(muReg[e, g]*grafo.V[e // 100 - 1][dim] + pReg[e, g]*(grafo.V[e % 100][dim] - grafo.V[e // 100 - 1][dim]) for e in grafo.aristas))

    MODEL.addConstrs(muLeg.sum('*', g) == 1 for g in T_index)
    MODEL.addConstrs(muReg.sum('*', g) == 1 for g in T_index)

    # Restricciones para ir de L a R
    for e, g in muReg.keys():

        first = e // 100 - 1
        second = e % 100

        MODEL.addConstr(maxeLRg[e, g] - mineLRg[e, g] == muxRg[e, g] - muxLg[e, g])
        MODEL.addConstr(maxeLRg[e, g] <= entryeLRg[e, g])
        MODEL.addConstr(mineLRg[e, g] <= 1 - entryeLRg[e, g])

        # MODEL.addConstr(biLRg[first, g] + biLRg[second, g] >= gp.quicksum(zeeLRg[e, e1, g] for e1 in grafo.aristas if e1 != e))
        # MODEL.addConstr(ciLRg[first, g] + ciLRg[second, g] >= gp.quicksum(zeeLRg[e1, e, g] for e1 in grafo.aristas if e1 != e))
        MODEL.addConstr(biLRg[first, g] + biLRg[second, g] >= muLeg[e, g])
        MODEL.addConstr(ciLRg[first, g] + ciLRg[second, g] >= muReg[e, g])



    for e, i, g in pxLbiLRg.keys():
        MODEL.addConstr(pxLbiLRg[e, i, g] >= muxLg[e, g] + biLRg[i, g] - 1)
        MODEL.addConstr(pxLbiLRg[e, i, g] <= muxLg[e, g])
        MODEL.addConstr(pxLbiLRg[e, i, g] <= biLRg[i, g])

        MODEL.addConstr(pxRciLRg[e, i, g] >= muxRg[e, g] + ciLRg[i, g] - 1)
        MODEL.addConstr(pxRciLRg[e, i, g] <= muxRg[e, g])
        MODEL.addConstr(pxRciLRg[e, i, g] <= ciLRg[i, g])


    for i, g in biLRg.keys():
        MODEL.addConstr(gp.quicksum(muLeg[e, g] for e in grafo.aristas if e // 100 - 1 == i or e % 100 == i) >= biLRg[i, g])
        MODEL.addConstr(gp.quicksum(muReg[e, g] for e in grafo.aristas if e // 100 - 1 == i or e % 100 == i) >= ciLRg[i, g])


    for e1, e2, g in deeLRg.keys():
        if e1 == e2:
            first_e1 = e1 // 100 - 1
            second_e1 = e1 % 100

            MODEL.addConstr(deeLRg[e1, e2, g] == (maxeLRg[e1, g] + mineLRg[e1, g])*grafo.longaristas[first_e1, second_e1])

        else:
            first_e1 = e1 // 100 - 1
            second_e1 = e1 % 100

            first_e2 = e2 // 100 - 1
            second_e2 = e2 % 100

            MODEL.addConstr(deeLRg[e1, e2, g] == (pxLbiLRg[e1, first_e1, g] + biLRg[second_e1, g] - pxLbiLRg[e1, second_e1, g])*grafo.longaristas[first_e1, second_e1] + gp.quicksum(qeLRg[grafo.aristas[ind]//100-1, grafo.aristas[ind]%100, g]*grafo.longitudes[ind] for ind in range(grafo.num_aristas)) + gp.quicksum(qeLRg[grafo.aristas[ind]%100, grafo.aristas[ind]//100-1, g]*grafo.longitudes[ind] for ind in range(grafo.num_aristas)) + (pxRciLRg[e2, first_e2, g] + ciLRg[second_e2, g] - pxRciLRg[e2, second_e2, g])*grafo.longaristas[first_e2, second_e2])

        # MODEL.addConstr(peeLRg[e1, e2, g] <= 10000* zeeLRg[e1, e2, g])
        # MODEL.addConstr(peeLRg[e1, e2, g] <= deeLRg[e1, e2, g])
        MODEL.addConstr(peeLRg[e1, e2, g] >= SmallM*zeeLRg[e1, e2, g])
        MODEL.addConstr(peeLRg[e1, e2, g] >= deeLRg[e1, e2, g] - 10000*(1 - zeeLRg[e1, e2, g]))

        MODEL.addConstr(zeeLRg[e1, e2, g] >= muLeg[e1, g] + muReg[e2, g] - 1)
        MODEL.addConstr(zeeLRg[e1, e2, g] <= muLeg[e1, g])
        MODEL.addConstr(zeeLRg[e1, e2, g] <= muReg[e2, g])

    for i, g in biLRg.keys():
        MODEL.addConstr(biLRg[i, g] + gp.quicksum(qeLRg[j, i, g] for j in range(grafo.num_puntos) if ((j+1) * 100 + i) in grafo.aristas or ((i+1) * 100 + j) in grafo.aristas) == gp.quicksum(qeLRg[i, j, g] for j in range(grafo.num_puntos) if ((j+1) * 100 + i) in grafo.aristas or ((i+1) * 100 + j) in grafo.aristas) + ciLRg[i, g])

    MODEL.addConstrs(dLR[g] == gp.quicksum(peeLRg[e1, e2, g] for e1 in grafo.aristas for e2 in grafo.aristas) for g in dLR.keys())





    # Restricciones para ir de R a L
    for e, g1, g2 in mineRLgg.keys():

        MODEL.addConstr(maxeRLgg[e, g1, g2] - mineRLgg[e, g1, g2] == muxRg[e, g1] - muxLg[e, g2])
        MODEL.addConstr(maxeRLgg[e, g1, g2] <= entryeRLgg[e, g1, g2])
        MODEL.addConstr(mineRLgg[e, g1, g2] <= 1 - entryeRLgg[e, g1, g2])

    for e, g in muLeg.keys():

        first = e // 100 - 1
        second = e % 100

        # MODEL.addConstr(biRLg[first, g] + biRLg[second, g] >= gp.quicksum(zeeRLgg[e, e1, g, g2] for e1 in grafo.aristas for g2 in T_index_prima if e1 != e and g != g2))
        MODEL.addConstr(biRLg[first, g] + biRLg[second, g] >= muReg[e, g])
        # MODEL.addConstr(biRLg[first, g2] + biRLg[second, g2] >= gp.quicksum(zeeRLgg[e, e1, g1, g2] for e1 in grafo.aristas if e1 != e))
        MODEL.addConstr(ciRLg[first, g] + ciRLg[second, g] >= muLeg[e, g])

        # MODEL.addConstr(ciRLg[first, g] + ciRLg[second, g] >= gp.quicksum(zeeRLgg[e1, e, g1, g] for e1 in grafo.aristas for g1 in T_index_prima if e1 != e and g != g1))
        # MODEL.addConstr(ciRLg[first, g2] + ciRLg[second, g2] >= gp.quicksum(zeeRLgg[e1, e, g1, g2] for e1 in grafo.aristas if e1 != e))


    for e, i, g in pxRciRLg.keys():
        MODEL.addConstr(pxRciRLg[e, i, g] >= muxRg[e, g] + biRLg[i, g] - 1)
        MODEL.addConstr(pxRciRLg[e, i, g] <= muxRg[e, g])
        MODEL.addConstr(pxRciRLg[e, i, g] <= biRLg[i, g])

        MODEL.addConstr(pxLbiRLg[e, i, g] >= muxLg[e, g] + ciRLg[i, g] - 1)
        MODEL.addConstr(pxLbiRLg[e, i, g] <= muxLg[e, g])
        MODEL.addConstr(pxLbiRLg[e, i, g] <= ciRLg[i, g])

    for i, g in biRLg.keys():
        MODEL.addConstr(gp.quicksum(muReg[e, g] for e in grafo.aristas if e // 100 - 1 == i or e % 100 == i) >= biRLg[i, g])
        MODEL.addConstr(gp.quicksum(muLeg[e, g] for e in grafo.aristas if e // 100 - 1 == i or e % 100 == i) >= ciRLg[i, g])

        # MODEL.addConstr(biRLg[first, g] + biRLg[second, g] <= zeeRLgg[e, e, g])
        # MODEL.addConstr(ciRLg[first, g] + ciRLg[second, g] <= zeeRLgg[e, e, g])



    for e1, e2, g1, g2 in deeRLgg.keys():

        MODEL.addConstr(zeeRLgg[e1, e2, g1, g2] >= muReg[e1, g1] + muLeg[e2, g2] - 1)
        MODEL.addConstr(zeeRLgg[e1, e2, g1, g2] <= muReg[e1, g1])
        MODEL.addConstr(zeeRLgg[e1, e2, g1, g2] <= muLeg[e2, g2])

        if e1 == e2:
            first_e1 = e1 // 100 - 1
            second_e1 = e1 % 100

            MODEL.addConstr(deeRLgg[e1, e2, g1, g2] == (maxeRLgg[e1, g1, g2] + mineRLgg[e1, g1, g2])*grafo.longaristas[first_e1, second_e1])

        else:
            first_e1 = e1 // 100 - 1
            second_e1 = e1 % 100

            first_e2 = e2 // 100 - 1
            second_e2 = e2 % 100

            MODEL.addConstr(deeRLgg[e1, e2, g1, g2] == (pxRciRLg[e1, first_e1, g1] + biRLg[second_e1, g1] - pxRciRLg[e1, second_e1, g1])*grafo.longaristas[first_e1, second_e1] + gp.quicksum(qeRLgg[grafo.aristas[ind]//100-1, grafo.aristas[ind]%100, g1, g2]*grafo.longitudes[ind] for ind in range(grafo.num_aristas)) + gp.quicksum(qeRLgg[grafo.aristas[ind]%100, grafo.aristas[ind]//100-1, g1, g2]*grafo.longitudes[ind] for ind in range(grafo.num_aristas)) + (pxLbiRLg[e2, first_e2, g2] + ciRLg[second_e2, g2] - pxLbiRLg[e2, second_e2, g2])*grafo.longaristas[first_e2, second_e2])

        MODEL.addConstr(peeRLgg[e1, e2, g1, g2] >= SmallM* zeeRLgg[e1, e2, g1, g2])
        # MODEL.addConstr(peeRLgg[e1, e2, g1, g2] <= deeRLgg[e1, e2, g1, g2])
        # MODEL.addConstr(peeRLgg[e1, e2, g1, g2] >= 0)
        MODEL.addConstr(peeRLgg[e1, e2, g1, g2] >= deeRLgg[e1, e2, g1, g2] - 10000*(1 - zeeRLgg[e1, e2, g1, g2]))



    for g1, g2 in z.keys():
        for i in range(grafo.num_puntos):
            MODEL.addConstr(biRLg[i, g1] + gp.quicksum(qeRLgg[j, i, g1, g2] for j in range(grafo.num_puntos) if ((j+1) * 100 + i) in grafo.aristas or ((i+1) * 100 + j) in grafo.aristas) == gp.quicksum(qeRLgg[i, j, g1, g2] for j in range(grafo.num_puntos) if ((j+1) * 100 + i) in grafo.aristas or ((i+1) * 100 + j) in grafo.aristas) + ciRLg[i, g2])

    # for g in T_index:
    #     for i in range(grafo.num_puntos):
    #         MODEL.addConstr(biRLg[i, g] + gp.quicksum(qeRLgg[j, i, g, g] for j in range(grafo.num_puntos) if ((j+1) * 100 + i) in grafo.aristas or ((i+1) * 100 + j) in grafo.aristas) == gp.quicksum(qeRLgg[i, j, g, g] for j in range(grafo.num_puntos) if ((j+1) * 100 + i) in grafo.aristas or ((i+1) * 100 + j) in grafo.aristas) + ciRLg[i, g])

    MODEL.addConstrs(dRL[g1, g2] == gp.quicksum(peeRLgg[e1, e2, g1, g2] for e1 in grafo.aristas for e2 in grafo.aristas) for g1, g2 in z.keys())



    # Origen y destino
    MODEL.addConstrs(xL[0, dim] == orig[dim] for dim in range(2))
    MODEL.addConstrs(xR[0, dim] == orig[dim] for dim in range(2))

    MODEL.addConstrs(xL[datos.m+1, dim] == dest[dim] for dim in range(2))
    MODEL.addConstrs(xR[datos.m+1, dim] == dest[dim] for dim in range(2))

    MODEL.update()

    objective = gp.quicksum(pgLi[g, i] + pgRi[g, i] for g, i in pgRi.keys()) + gp.quicksum(pgij[g, i, j] for g, i, j in pgij.keys()) + gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for g in T_index_prima for i in grafos[g-1].aristas) + gp.quicksum(3*dLR[g] for g in dLR.keys()) + gp.quicksum(3*pRL[g1, g2] for g1, g2 in dRL.keys())

    # objective = gp.quicksum(1*dLR[g] for g in dLR.keys()) + gp.quicksum(1*pRL[g1, g2] for g1, g2 in dRL.keys()) # + gp.quicksum(1e5*holg[g] for g in holg.keys())

    # objective = gp.quicksum(dRL[g] + dLR[g] for t in T_index)

    MODEL.setObjective(objective, GRB.MINIMIZE)
    MODEL.Params.Threads = 6
    # MODEL.Params.NumericFocus = 3
    # MODEL.Params.Heuristics = 0
    # MODEL.Params.NonConvex = 2
    MODEL.Params.OutputFlag = 1
    MODEL.Params.timeLimit = 15
    # MODEL.Params.FeasibilityTol = 1e-2

    MODEL.update()
    # MODEL.computeIIS()
    # MODEL.write('casa.ilp')

    # MODEL.write('AMDRPG-MTZ.lp')
    # MODEL.write('AMDRPG-MTZ.mps')
    MODEL.write('./models/NDMTZ-partial' + str(index_g)+ '.lp')

    MODEL.optimize()



    # MODEL.update()

    if MODEL.Status == 3:
        MODEL.computeIIS()
        MODEL.write('./models/IIS_NDMTZ.ilp')
        result =  [np.nan, np.nan, np.nan, np.nan]
        if datos.grid:
            result.append('Grid')
        else:
            result.append('Delauney')

        if datos.alpha:
            result.append('Alpha-e')
        else:
            result.append('Alpha-g')

        result.append('MTZ')
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
#
#        result.append('MTZ')
#        result.append(grafo_data.mode)

        return result


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

    result.append('MTZ')
    result.append(grafo_data.mode)

    MODEL.write('./soluciones/NDMTZ-Heur{0}.sol'.format(index_g))

    vals_u = MODEL.getAttr('x', ugi)
    selected_u = gp.tuplelist((g, i)
                              for g, i in vals_u.keys() if vals_u[g, i] > 0.5)
    # print(selected_u)

    vals_zgij = MODEL.getAttr('x', zgij)
    selected_zgij = gp.tuplelist((g, i, j)
                              for g, i, j in vals_zgij.keys() if vals_zgij[g, i, j] > 0.5)
    # print(selected_zgij)

    vals_v = MODEL.getAttr('x', vgi)
    selected_v = gp.tuplelist((g, i)
                              for g, i in vals_v.keys() if vals_v[g, i] > 0.5)
    # print(selected_v)

    # valsz = MODEL.getAttr('x', z)
    # selected_z = gp.tuplelist(e for e in valsz if valsz[e] > 0)
    # # print(selected_z)
    # path = af.subtour(selected_z)
    # path.append(datos.m+1)
    # print(path)

    path = [0, index_g, nG+1]

    # vals_dLR = MODEL.getAttr('x', dLR)
    # print(vals_dLR)

    # vals_pRL = MODEL.getAttr('x', pRL)
    # selected_pRL = gp.tuplelist(vals_pRL[e] for e in vals_pRL.keys() if vals_pRL[e] > 0.5)
    # print(selected_pRL)

    ind = 0
    path_C = []
    paths_D = []

    for p in path:
        path_C.append([xL[p, 0].X, xL[p, 1].X])
        path_C.append([xR[p, 0].X, xR[p, 1].X])


    path_D = []
    path_D.append([xL[index_g, 0].X, xL[index_g, 1].X])
    index_gr, index_i = selected_u[0]
    count = 0
    path_D.append([Rgi[index_gr, index_i, 0].X, Rgi[index_gr, index_i, 1].X])
    path_D.append([Lgi[index_gr, index_i, 0].X, Lgi[index_gr, index_i, 1].X])
    limite = sum([1 for g, i, j in selected_zgij if g == index_gr])
    while count < limite:
        for g, i, j in selected_zgij:
            if index_gr == g and index_i == i:
                count += 1
                index_i = j
                path_D.append([Rgi[index_gr, index_i, 0].X, Rgi[index_gr, index_i, 1].X])
                path_D.append([Lgi[index_gr, index_i, 0].X, Lgi[index_gr, index_i, 1].X])

    ind += 1
    path_D.append([xR[index_g, 0].X, xR[index_g, 1].X])
    paths_D.append(path_D)

    # path_C.append([xL[datos.m+1, 0].X, xL[datos.m+1, 1].X])
    fig, ax = plt.subplots()
    plt.axis([0, 100, 0, 100])

    for g, i in rhogi.keys():
        first = i // 100 - 1
        second = i % 100
        if mugi[g, i].X > 0.5:
            plt.plot(Rgi[g, i, 0].X, Rgi[g, i, 1].X, 'kD', markersize=1, color='red')
            # ax.annotate("$R_" + str(g) + "^{" + str((first, second)) + "}$", xy = (Rgi[g, i, 0].X+1, Rgi[g, i, 1].X+1))
            plt.plot(Lgi[g, i, 0].X, Lgi[g, i, 1].X, 'kD', markersize=1, color='red')
            # ax.annotate("$L_" + str(g) + "^{" + str((first, second)) + "}$", xy = (Lgi[g, i, 0].X-3, Lgi[g, i, 1].X-3))

    #
    for p, i in zip(path, range(len(path))):
        # path_C.append([xL[t, 0].X, xL[t, 1].X])
        # path_C.append([xR[t, 0].X, xR[t, 1].X])
        plt.plot(xL[p, 0].X, xL[p, 1].X, 'ko', alpha = 0.3, markersize=10, color='green')
        #ax.annotate("L" + str(i), xy = (xL[p, 0].X+0.5, xL[p, 1].X+0.5))
        plt.plot(xR[p, 0].X, xR[p, 1].X, 'ko', markersize=5, color='blue')
        #ax.annotate("R" + str(i), xy = (xR[p, 0].X-1.5, xR[p, 1].X-1.5))


    ax.add_artist(Polygon(path_C, fill=False, closed = False, animated=False,
                  linestyle='-', alpha=1, color='blue'))

    for pathd in paths_D:
        ax.add_artist(Polygon(pathd, fill=False, closed=False, lw = 0.1,
                      animated=False, alpha=1, color='red'))
    #
    # ax.add_artist(Polygon(path_D, fill=False, animated=False,
    #               linestyle='dotted', alpha=1, color='red'))

    for g in range(nG):
        grafo = grafos[g]
        nx.draw(grafo.G, grafo.pos, node_size=40,
                node_color='black', alpha=1, edge_color='black')

    plt.savefig('NDMTZ-heur{0}.png'.format(index_g))

    return MODEL.getAttr('x', xL), MODEL.getAttr('x', xR),  MODEL.ObjVal

def NDMTZ_aux2(datos, grafo_data, fixed_data, iter):

    mode = len(fixed_data)

    if mode == 1:
        print()
        print('--------------------------------------------')
        print('Exact Formulation: Fixing w. Iteration: {0}'.format(iter))
        print('--------------------------------------------')
        print()

    if mode == 2:
        print()
        print('--------------------------------------------')
        print('Exact Formulation: Fixing points. Iteration: {i})'.format(i = iter))
        print('--------------------------------------------')
        print()

    orig = grafo_data.orig
    dest = grafo_data.dest
    grafos = datos.mostrar_datos()
    # print(grafos[0].V)

    result = []

    nG = datos.m

    T_index = range(nG + 2)
    T_index_prima = range(1, nG+1)
    T_index_primaprima = range(nG+1)

    vD = datos.vD
    vC = datos.vC

    grafo = grafo_data.data[0]


    # Creamos el modelo NDMTZ
    MODEL = gp.Model("NDMTZ_aux2")

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
    dLR = MODEL.addVars(T_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dLR')

    # ParamuLegrizacion del grafo
    muxLg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'muxLg')
    mineLRg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'mineLRg')
    maxeLRg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'maxeLRg')

    mineRLgg = MODEL.addVars(grafo.aristas, z_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'mineRLgg')
    maxeRLgg = MODEL.addVars(grafo.aristas, z_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'maxeRLgg')

    #aet = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'aet')
    entryeLRg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.BINARY, name = 'entryeLRg')
    entryeRLgg = MODEL.addVars(grafo.aristas, z_index, vtype = GRB.BINARY, name = 'entryeRLgg')

    # Variables binaria que es 1 si en la etapa t el punto xL se encuentra en el segmento e
    muLeg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.BINARY, name = 'muLeg')

    pLeg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'pLeg')

    # Variable binaria que toma el valor 1 cuando el camion entra al vertice i
    biLRg = MODEL.addVars(grafo.num_puntos, T_index, vtype = GRB.BINARY, name = 'biLRg')
    biRLg = MODEL.addVars(grafo.num_puntos, T_index, vtype = GRB.BINARY, name = 'biRLg')

    # Variable que modela el producto de muxLg con biLRg
    pxLbiLRg = MODEL.addVars(grafo.aristas, grafo.num_puntos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pxLbiLRg')
    pxLbiRLg = MODEL.addVars(grafo.aristas, grafo.num_puntos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pxLbiRLg')


    # Variable binaria que toma el valor 1 cuando el camion atraviesa la arista e
    qeLRg = MODEL.addVars(grafo.num_puntos, grafo.num_puntos, T_index, vtype = GRB.INTEGER, name = 'qeLRg')
    qeRLgg = MODEL.addVars(grafo.num_puntos, grafo.num_puntos, z_index, vtype = GRB.INTEGER, name = 'qeRLgg')

    muxRg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'muxRg')

    # Variables binaria que es 1 si en la etapa t el punto xR se encuentra en el segmento e
    muReg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.BINARY, name = 'muReg')

    pReg = MODEL.addVars(grafo.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'pReg')

    # Variable binaria que toma el valor 1 cuando el camion entra al vertice i
    ciLRg = MODEL.addVars(grafo.num_puntos, T_index, vtype = GRB.BINARY, name = 'ciLRg')
    ciRLg = MODEL.addVars(grafo.num_puntos, T_index, vtype = GRB.BINARY, name = 'ciRLg')

    # Variable que modela el producto de muxLg con muLjt
    pxRciLRg = MODEL.addVars(grafo.aristas, grafo.num_puntos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pxRciLRg')
    pxRciRLg = MODEL.addVars(grafo.aristas, grafo.num_puntos, T_index, vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'pxRciRLg')


    zeeLRg = MODEL.addVars(grafo.aristas, grafo.aristas, T_index, vtype = GRB.BINARY, name = 'zeeLRg')
    deeLRg = MODEL.addVars(grafo.aristas, grafo.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'deeLRg')
    peeLRg = MODEL.addVars(grafo.aristas, grafo.aristas, T_index, vtype = GRB.CONTINUOUS, name = 'peeLRg')

    zeeRLgg = MODEL.addVars(grafo.aristas, grafo.aristas, z_index, vtype = GRB.BINARY, name = 'zeeRLgg')
    deeRLgg = MODEL.addVars(grafo.aristas, grafo.aristas, z_index, vtype = GRB.CONTINUOUS, name = 'deeRLgg')
    peeRLgg = MODEL.addVars(grafo.aristas, grafo.aristas, z_index, vtype = GRB.CONTINUOUS, name = 'peeRLgg')

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

    for g in T_index_prima:
        sols = MODEL.read('./soluciones/NDMTZ-Heur{0}.sol'.format(g))

    # Para cada grafo g, existe un segmento i y una etapa t donde hay que recoger al dron
    MODEL.addConstrs((ugi.sum(g, '*') == 1 for g in T_index_prima), name = 'entrag')
    MODEL.addConstrs((vgi.sum(g, '*') == 1 for g in T_index_prima), name = 'saleg')

    # MODEL.addConstrs(ugi.sum('*', i, '*') == 1 for i in range(datos.m))
    # MODEL.addConstrs(vgi.sum('*', i, '*') == 1 for g in range(datos.m))

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

    # for g in range(datos.m):
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
    #         for w in grafo.V:
    #             BigM = max([np.linalg.norm(w - v), BigM])
    #             SmallM = min([np.linalg.norm(w - v), SmallM])
    SmallM = 0
    BigM = 10000

    for g in T_index_prima:


        # BigM = max(max([np.linalg.norm(v - w) for v in grafo.V for w in grafos[g-1].V]), BigM)
        # SmallM = min(min([np.linalg.norm(v - w) for v in grafo.V for w in grafos[g-1].V]), SmallM)

        MODEL.addConstrs((pgLi[g, i] >= SmallM * ugi[g, i]) for i in grafos[g-1].aristas)
        MODEL.addConstrs((pgLi[g, i] >= dgLi[g, i] - BigM * (1 - ugi[g, i])) for i in grafos[g-1].aristas)

        MODEL.addConstrs((pgRi[g, i] >= SmallM * vgi[g, i]) for i in grafos[g-1].aristas)
        MODEL.addConstrs((pgRi[g, i] >= dgRi[g, i] - BigM * (1 - vgi[g, i])) for i in grafos[g-1].aristas)

    MODEL.addConstrs((difgij[g, i, j, dim] >=   Lgi[g, i, dim] - Rgi[g, j, dim]) for g, i, j, dim in difgij.keys())
    MODEL.addConstrs((difgij[g, i, j, dim] >= - Lgi[g, i, dim] + Rgi[g, j, dim]) for g, i, j, dim in difgij.keys())

    MODEL.addConstrs((difgij[g, i, j, 0]*difgij[g, i, j, 0] + difgij[g, i, j, 1] * difgij[g, i, j, 1] <= dgij[g, i, j] * dgij[g, i, j] for g, i, j in dgij.keys()), name = 'zgij-conic')


    for g, i, j in zgij.keys():
        first_i = i // 100 - 1
        second_i = i % 100
        first_j = j // 100 - 1
        second_j = j % 100

        segm_i = Poligonal(np.array([[grafos[g-1].V[first_i, 0], grafos[g-1].V[first_i, 1]], [grafos[g-1].V[second_i, 0], grafos[g-1].V[second_i, 1]]]), grafos[g-1].A[first_i, second_i])
        segm_j = Poligonal(np.array([[grafos[g-1].V[first_j, 0], grafos[g-1].V[first_j, 1]], [grafos[g-1].V[second_j, 0], grafos[g-1].V[second_j, 1]]]), grafos[g-1].A[first_j, second_j])

        BigM_local = eM.estima_BigM_local(segm_i, segm_j)
        SmallM_local = eM.estima_SmallM_local(segm_i, segm_j)
        MODEL.addConstr((pgij[g, i, j] >= SmallM_local * zgij[g, i, j]))
        MODEL.addConstr((pgij[g, i, j] >= dgij[g, i, j] - BigM_local * (1 - zgij[g, i, j])))

    MODEL.addConstrs((difgRi[g, i, dim] >=   Lgi[g, i, dim] - xR[g, dim]) for g, i, dim in difgRi.keys())
    MODEL.addConstrs((difgRi[g, i, dim] >= - Lgi[g, i, dim] + xR[g, dim]) for g, i, dim in difgRi.keys())

    MODEL.addConstrs((difgRi[g, i, 0]*difgRi[g, i, 0] + difgRi[g, i, 1] * difgRi[g, i, 1] <= dgRi[g, i] * dgRi[g, i] for g, i in vgi.keys()), name = 'v-conic')


    # SmallM = 0
    #BigM = 10000


    if mode == 2:
        # BigM_z = np.zeros((nG+2, nG+2))
        #
        # for g1, g2 in z_index:
        #     if g1 == 0 and g2 < nG+1:
        #         BigM_z[g1, g2] = max([np.linalg.norm(datos.orig - v) for v in grafos[g2-1].V])
        #     elif g2 == 0 and g1 < nG+1:
        #         BigM_z[g1, g2] = max([np.linalg.norm(datos.orig - v) for v in grafos[g1-1].V])
        #     elif g1 == nG+1 and g2 > 0:
        #         BigM_z[g1, g2] = max([np.linalg.norm(datos.dest - v) for v in grafos[g2-1].V])
        #     elif g2 == nG+1 and g1 > 0:
        #         BigM_z[g1, g2] = max([np.linalg.norm(datos.dest - v) for v in grafos[g1-1].V])
        #     if g1 > 0 and g2 > 0 and g1 < nG+1 and g2 < nG+1:
        #         BigM_z[g1, g2] = max([np.linalg.norm(u - v) for u in grafos[g1-1].V for v in grafos[g2-1].V])

        SmallM = 0
        BigM = 10000

        MODEL.addConstrs((pRL[g1, g2] >= SmallM * z[g1, g2] for g1, g2 in z_index))
        MODEL.addConstrs((pRL[g1, g2] >= dRL[g1, g2] - BigM * (1 - z[g1, g2]) for g1, g2 in z_index))

        # Restricciones para formar un tour
        MODEL.addConstr(z.sum('*', 0) == 0)
        MODEL.addConstr(z.sum(datos.m+1, '*') == 0)
        MODEL.addConstr(z.sum(0, '*') - z[0, datos.m+1] == 1)
        MODEL.addConstr(z.sum('*', datos.m+1) - z[0, datos.m+1] == 1)
        MODEL.addConstrs(z.sum(v, '*') == 1 for v in T_index_prima)
        MODEL.addConstrs(z.sum('*', v) == 1 for v in T_index_prima)

        # Conectividad
        for v in T_index_prima:
            for w in T_index_prima:
                if v != w:
                    MODEL.addConstr(len(T_index) - 1 >= (s[v] - s[w]) + len(T_index) * z[v, w])

        for v in T_index_prima:
            MODEL.addConstr(s[v] >= 1)
            MODEL.addConstr(s[v] <= len(T_index) - 1)

            MODEL.addConstr(s[0] == 0)
            MODEL.addConstr(s[datos.m + 1] == datos.m+1)


    MODEL.addConstrs((pgLi.sum(g, '*') + pgij.sum(g, '*', '*') +  gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas) + pgRi.sum(g, '*'))/vD <= dLR[g]/vC for g in T_index_prima)

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

    for e, g in pLeg.keys():
        MODEL.addConstr(pLeg[e, g] >= muLeg[e, g] + muxLg[e, g] - 1)
        MODEL.addConstr(pLeg[e, g] <= muLeg[e, g])
        MODEL.addConstr(pLeg[e, g] <= muxLg[e, g])

        MODEL.addConstr(pReg[e, g] >= muReg[e, g] + muxRg[e, g] - 1)
        MODEL.addConstr(pReg[e, g] <= muReg[e, g])
        MODEL.addConstr(pReg[e, g] <= muxRg[e, g])

    for g, dim in xL.keys():
        MODEL.addConstr(xL[g, dim] == gp.quicksum(muLeg[e, g]*grafo.V[e // 100 - 1][dim] + pLeg[e, g]*(grafo.V[e % 100][dim] - grafo.V[e // 100 - 1][dim]) for e in grafo.aristas))
        MODEL.addConstr(xR[g, dim] == gp.quicksum(muReg[e, g]*grafo.V[e // 100 - 1][dim] + pReg[e, g]*(grafo.V[e % 100][dim] - grafo.V[e // 100 - 1][dim]) for e in grafo.aristas))

    MODEL.addConstrs(muLeg.sum('*', g) == 1 for g in T_index)
    MODEL.addConstrs(muReg.sum('*', g) == 1 for g in T_index)

    # Restricciones para ir de L a R
    for e, g in muReg.keys():

        first = e // 100 - 1
        second = e % 100

        MODEL.addConstr(maxeLRg[e, g] - mineLRg[e, g] == muxRg[e, g] - muxLg[e, g])
        MODEL.addConstr(maxeLRg[e, g] <= entryeLRg[e, g])
        MODEL.addConstr(mineLRg[e, g] <= 1 - entryeLRg[e, g])

        # MODEL.addConstr(biLRg[first, g] + biLRg[second, g] >= gp.quicksum(zeeLRg[e, e1, g] for e1 in grafo.aristas if e1 != e))
        # MODEL.addConstr(ciLRg[first, g] + ciLRg[second, g] >= gp.quicksum(zeeLRg[e1, e, g] for e1 in grafo.aristas if e1 != e))
        MODEL.addConstr(biLRg[first, g] + biLRg[second, g] >= muLeg[e, g])
        MODEL.addConstr(ciLRg[first, g] + ciLRg[second, g] >= muReg[e, g])



    for e, i, g in pxLbiLRg.keys():
        MODEL.addConstr(pxLbiLRg[e, i, g] >= muxLg[e, g] + biLRg[i, g] - 1)
        MODEL.addConstr(pxLbiLRg[e, i, g] <= muxLg[e, g])
        MODEL.addConstr(pxLbiLRg[e, i, g] <= biLRg[i, g])

        MODEL.addConstr(pxRciLRg[e, i, g] >= muxRg[e, g] + ciLRg[i, g] - 1)
        MODEL.addConstr(pxRciLRg[e, i, g] <= muxRg[e, g])
        MODEL.addConstr(pxRciLRg[e, i, g] <= ciLRg[i, g])


    for i, g in biLRg.keys():
        MODEL.addConstr(gp.quicksum(muLeg[e, g] for e in grafo.aristas if e // 100 - 1 == i or e % 100 == i) >= biLRg[i, g])
        MODEL.addConstr(gp.quicksum(muReg[e, g] for e in grafo.aristas if e // 100 - 1 == i or e % 100 == i) >= ciLRg[i, g])


    for e1, e2, g in deeLRg.keys():
        if e1 == e2:
            first_e1 = e1 // 100 - 1
            second_e1 = e1 % 100

            MODEL.addConstr(deeLRg[e1, e2, g] == (maxeLRg[e1, g] + mineLRg[e1, g])*grafo.longaristas[first_e1, second_e1])

        else:
            first_e1 = e1 // 100 - 1
            second_e1 = e1 % 100

            first_e2 = e2 // 100 - 1
            second_e2 = e2 % 100

            MODEL.addConstr(deeLRg[e1, e2, g] == (pxLbiLRg[e1, first_e1, g] + biLRg[second_e1, g] - pxLbiLRg[e1, second_e1, g])*grafo.longaristas[first_e1, second_e1] + gp.quicksum(qeLRg[grafo.aristas[ind]//100-1, grafo.aristas[ind]%100, g]*grafo.longitudes[ind] for ind in range(grafo.num_aristas)) + gp.quicksum(qeLRg[grafo.aristas[ind]%100, grafo.aristas[ind]//100-1, g]*grafo.longitudes[ind] for ind in range(grafo.num_aristas)) + (pxRciLRg[e2, first_e2, g] + ciLRg[second_e2, g] - pxRciLRg[e2, second_e2, g])*grafo.longaristas[first_e2, second_e2])

        # MODEL.addConstr(peeLRg[e1, e2, g] <= 10000* zeeLRg[e1, e2, g])
        # MODEL.addConstr(peeLRg[e1, e2, g] <= deeLRg[e1, e2, g])
        MODEL.addConstr(peeLRg[e1, e2, g] >= SmallM*zeeLRg[e1, e2, g])
        MODEL.addConstr(peeLRg[e1, e2, g] >= deeLRg[e1, e2, g] - 10000*(1 - zeeLRg[e1, e2, g]))

        MODEL.addConstr(zeeLRg[e1, e2, g] >= muLeg[e1, g] + muReg[e2, g] - 1)
        MODEL.addConstr(zeeLRg[e1, e2, g] <= muLeg[e1, g])
        MODEL.addConstr(zeeLRg[e1, e2, g] <= muReg[e2, g])

    for i, g in biLRg.keys():
        MODEL.addConstr(biLRg[i, g] + gp.quicksum(qeLRg[j, i, g] for j in range(grafo.num_puntos) if ((j+1) * 100 + i) in grafo.aristas or ((i+1) * 100 + j) in grafo.aristas) == gp.quicksum(qeLRg[i, j, g] for j in range(grafo.num_puntos) if ((j+1) * 100 + i) in grafo.aristas or ((i+1) * 100 + j) in grafo.aristas) + ciLRg[i, g])

    MODEL.addConstrs(dLR[g] == gp.quicksum(peeLRg[e1, e2, g] for e1 in grafo.aristas for e2 in grafo.aristas) for g in dLR.keys())





    # Restricciones para ir de R a L
    for e, g1, g2 in mineRLgg.keys():

        MODEL.addConstr(maxeRLgg[e, g1, g2] - mineRLgg[e, g1, g2] == muxRg[e, g1] - muxLg[e, g2])
        MODEL.addConstr(maxeRLgg[e, g1, g2] <= entryeRLgg[e, g1, g2])
        MODEL.addConstr(mineRLgg[e, g1, g2] <= 1 - entryeRLgg[e, g1, g2])

    for e, g in muLeg.keys():

        first = e // 100 - 1
        second = e % 100

        # MODEL.addConstr(biRLg[first, g] + biRLg[second, g] >= gp.quicksum(zeeRLgg[e, e1, g, g2] for e1 in grafo.aristas for g2 in T_index_prima if e1 != e and g != g2))
        MODEL.addConstr(biRLg[first, g] + biRLg[second, g] >= muReg[e, g])
        # MODEL.addConstr(biRLg[first, g2] + biRLg[second, g2] >= gp.quicksum(zeeRLgg[e, e1, g1, g2] for e1 in grafo.aristas if e1 != e))
        MODEL.addConstr(ciRLg[first, g] + ciRLg[second, g] >= muLeg[e, g])

        # MODEL.addConstr(ciRLg[first, g] + ciRLg[second, g] >= gp.quicksum(zeeRLgg[e1, e, g1, g] for e1 in grafo.aristas for g1 in T_index_prima if e1 != e and g != g1))
        # MODEL.addConstr(ciRLg[first, g2] + ciRLg[second, g2] >= gp.quicksum(zeeRLgg[e1, e, g1, g2] for e1 in grafo.aristas if e1 != e))


    for e, i, g in pxRciRLg.keys():
        MODEL.addConstr(pxRciRLg[e, i, g] >= muxRg[e, g] + biRLg[i, g] - 1)
        MODEL.addConstr(pxRciRLg[e, i, g] <= muxRg[e, g])
        MODEL.addConstr(pxRciRLg[e, i, g] <= biRLg[i, g])

        MODEL.addConstr(pxLbiRLg[e, i, g] >= muxLg[e, g] + ciRLg[i, g] - 1)
        MODEL.addConstr(pxLbiRLg[e, i, g] <= muxLg[e, g])
        MODEL.addConstr(pxLbiRLg[e, i, g] <= ciRLg[i, g])

    for i, g in biRLg.keys():
        MODEL.addConstr(gp.quicksum(muReg[e, g] for e in grafo.aristas if e // 100 - 1 == i or e % 100 == i) >= biRLg[i, g])
        MODEL.addConstr(gp.quicksum(muLeg[e, g] for e in grafo.aristas if e // 100 - 1 == i or e % 100 == i) >= ciRLg[i, g])

        # MODEL.addConstr(biRLg[first, g] + biRLg[second, g] <= zeeRLgg[e, e, g])
        # MODEL.addConstr(ciRLg[first, g] + ciRLg[second, g] <= zeeRLgg[e, e, g])



    for e1, e2, g1, g2 in deeRLgg.keys():

        MODEL.addConstr(zeeRLgg[e1, e2, g1, g2] >= muReg[e1, g1] + muLeg[e2, g2] - 1)
        MODEL.addConstr(zeeRLgg[e1, e2, g1, g2] <= muReg[e1, g1])
        MODEL.addConstr(zeeRLgg[e1, e2, g1, g2] <= muLeg[e2, g2])

        if e1 == e2:
            first_e1 = e1 // 100 - 1
            second_e1 = e1 % 100

            MODEL.addConstr(deeRLgg[e1, e2, g1, g2] == (maxeRLgg[e1, g1, g2] + mineRLgg[e1, g1, g2])*grafo.longaristas[first_e1, second_e1])

        else:
            first_e1 = e1 // 100 - 1
            second_e1 = e1 % 100

            first_e2 = e2 // 100 - 1
            second_e2 = e2 % 100

            MODEL.addConstr(deeRLgg[e1, e2, g1, g2] == (pxRciRLg[e1, first_e1, g1] + biRLg[second_e1, g1] - pxRciRLg[e1, second_e1, g1])*grafo.longaristas[first_e1, second_e1] + gp.quicksum(qeRLgg[grafo.aristas[ind]//100-1, grafo.aristas[ind]%100, g1, g2]*grafo.longitudes[ind] for ind in range(grafo.num_aristas)) + gp.quicksum(qeRLgg[grafo.aristas[ind]%100, grafo.aristas[ind]//100-1, g1, g2]*grafo.longitudes[ind] for ind in range(grafo.num_aristas)) + (pxLbiRLg[e2, first_e2, g2] + ciRLg[second_e2, g2] - pxLbiRLg[e2, second_e2, g2])*grafo.longaristas[first_e2, second_e2])

        MODEL.addConstr(peeRLgg[e1, e2, g1, g2] >= SmallM* zeeRLgg[e1, e2, g1, g2])
        # MODEL.addConstr(peeRLgg[e1, e2, g1, g2] <= deeRLgg[e1, e2, g1, g2])
        # MODEL.addConstr(peeRLgg[e1, e2, g1, g2] >= 0)
        MODEL.addConstr(peeRLgg[e1, e2, g1, g2] >= deeRLgg[e1, e2, g1, g2] - 10000*(1 - zeeRLgg[e1, e2, g1, g2]))



    for g1, g2 in z.keys():
        for i in range(grafo.num_puntos):
            MODEL.addConstr(biRLg[i, g1] + gp.quicksum(qeRLgg[j, i, g1, g2] for j in range(grafo.num_puntos) if ((j+1) * 100 + i) in grafo.aristas or ((i+1) * 100 + j) in grafo.aristas) == gp.quicksum(qeRLgg[i, j, g1, g2] for j in range(grafo.num_puntos) if ((j+1) * 100 + i) in grafo.aristas or ((i+1) * 100 + j) in grafo.aristas) + ciRLg[i, g2])

    # for g in T_index:
    #     for i in range(grafo.num_puntos):
    #         MODEL.addConstr(biRLg[i, g] + gp.quicksum(qeRLgg[j, i, g, g] for j in range(grafo.num_puntos) if ((j+1) * 100 + i) in grafo.aristas or ((i+1) * 100 + j) in grafo.aristas) == gp.quicksum(qeRLgg[i, j, g, g] for j in range(grafo.num_puntos) if ((j+1) * 100 + i) in grafo.aristas or ((i+1) * 100 + j) in grafo.aristas) + ciRLg[i, g])

    MODEL.addConstrs(dRL[g1, g2] == gp.quicksum(peeRLgg[e1, e2, g1, g2] for e1 in grafo.aristas for e2 in grafo.aristas) for g1, g2 in z.keys())



    # Origen y destino
    MODEL.addConstrs(xL[0, dim] == orig[dim] for dim in range(2))
    MODEL.addConstrs(xR[0, dim] == orig[dim] for dim in range(2))

    MODEL.addConstrs(xL[datos.m+1, dim] == dest[dim] for dim in range(2))
    MODEL.addConstrs(xR[datos.m+1, dim] == dest[dim] for dim in range(2))

    MODEL.update()

    if mode == 1:
        z = fixed_data[0]
        objective = gp.quicksum(pgLi[g, i] + pgRi[g, i] for g, i in pgRi.keys()) + gp.quicksum(pgij[g, i, j] for g, i, j in pgij.keys()) + gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for g in T_index_prima for i in grafos[g-1].aristas) + gp.quicksum(3*dLR[g] for g in dLR.keys()) + gp.quicksum(3*dRL[g1, g2]*z[g1, g2] for g1, g2 in dRL.keys())

        MODEL.setObjective(objective, GRB.MINIMIZE)

    if mode == 2:
        objective = gp.quicksum(pgLi[g, i] + pgRi[g, i] for g, i in pgRi.keys()) + gp.quicksum(pgij[g, i, j] for g, i, j in pgij.keys()) + gp.quicksum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for g in T_index_prima for i in grafos[g-1].aristas) + gp.quicksum(3*dLR[g] for g in dLR.keys()) + gp.quicksum(3*pRL[g1, g2] for g1, g2 in dRL.keys())

    # objective = gp.quicksum(1*dLR[g] for g in dLR.keys()) + gp.quicksum(1*pRL[g1, g2] for g1, g2 in dRL.keys()) # + gp.quicksum(1e5*holg[g] for g in holg.keys())

    # objective = gp.quicksum(dRL[g] + dLR[g] for t in T_index)

    MODEL.setObjective(objective, GRB.MINIMIZE)
    # MODEL.Params.NumericFocus = 3
    MODEL.Params.Threads = 6
    # MODEL.Params.Heuristic = 0
    # MODEL.Params.NonConvex = 2
    MODEL.Params.timeLimit = 3600
    if mode == 1:
        MODEL.Params.SolutionLimit = 5
    if mode == 2:
        MODEL.Params.SolutionLimit = 5

    MODEL.update()
    # MODEL.computeIIS()
    # MODEL.write('casa.ilp')

    MODEL.write('NDMTZ-heuristic.lp')
    # MODEL.write('NDMTZ.mps')
    MODEL.optimize()



    # MODEL.update()

    if MODEL.Status == 3:
        MODEL.computeIIS()
        MODEL.write('IIS_NDMTZ.ilp')
        result =  [np.nan, np.nan, np.nan, np.nan]
        if datos.grid:
            result.append('Grid')
        else:
            result.append('Delauney')

        if datos.alpha:
            result.append('Alpha-e')
        else:
            result.append('Alpha-g')

        result.append('MTZ')
#        result.append(grafo_data.mode)

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

#        result.append('MTZ')
#        result.append(grafo_data.mode)

        return result


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

    result.append('MTZ')
    result.append(grafo_data.mode)

    MODEL.write('./soluciones/NDMTZ-HeurI{0}.sol'.format(iter))


    vals_u = MODEL.getAttr('x', ugi)
    selected_u = gp.tuplelist((g, i)
                              for g, i in vals_u.keys() if vals_u[g, i] > 0.5)
    # print(selected_u)

    vals_zgij = MODEL.getAttr('x', zgij)
    selected_zgij = gp.tuplelist((g, i, j)
                              for g, i, j in vals_zgij.keys() if vals_zgij[g, i, j] > 0.5)
    # print(selected_zgij)

    vals_v = MODEL.getAttr('x', vgi)
    selected_v = gp.tuplelist((g, i)
                              for g, i in vals_v.keys() if vals_v[g, i] > 0.5)
    # print(selected_v)

    if mode == 2:
        valsz = MODEL.getAttr('x', z)
        selected_z = gp.tuplelist(e for e in valsz if valsz[e] > 0)
        path = af.subtour(selected_z)
        path.append(datos.m+1)
        print(path)


    if mode == 1:
        selected_z = gp.tuplelist(e for e in dRL.keys() if z[e] > 0)
        path = af.subtour(selected_z)


#    ind = 0
#    path_C = []
#    paths_D = []
#
#    for p in path:
#        if p > 0 and p < len(path)-1:
#            path_C.append([xL[p, 0].X, xL[p, 1].X])
#            path_C.append([xR[p, 0].X, xR[p, 1].X])
#        else:
#            path_C.append([xL[p, 0].X, xL[p, 1].X])
#
#
#    for p in path[1:]:
#        #    if ind < nG:
#        if ind < nG:
#            path_D = []
#            path_D.append([xL[p, 0].X, xL[p, 1].X])
#            index_g = 0
#            index_i = 0
#            for g, i in selected_u:
#                if g == p:
#                    index_g = g
#                    index_i = i
#
#            count = 0
#            path_D.append([Rgi[index_g, index_i, 0].X, Rgi[index_g, index_i, 1].X])
#            path_D.append([Lgi[index_g, index_i, 0].X, Lgi[index_g, index_i, 1].X])
#            limite = sum([1 for g, i, j in selected_zgij if g == index_g])
#            while count < limite:
#                for g, i, j in selected_zgij:
#                    if index_g == g and index_i == i:
#                        count += 1
#                        index_i = j
#                        path_D.append([Rgi[index_g, index_i, 0].X, Rgi[index_g, index_i, 1].X])
#                        path_D.append([Lgi[index_g, index_i, 0].X, Lgi[index_g, index_i, 1].X])
#
#            ind += 1
#            path_D.append([xR[p, 0].X, xR[p, 1].X])
#        paths_D.append(path_D)
#
#    # path_C.append([xLt[nG+1, 0].X, xLt[nG+1, 1].X])
#
#    fig, ax = plt.subplots()
#    plt.axis([0, 100, 0, 100])
#    ax.set_aspect('equal')
#
#    for g, i in rhogi.keys():
#        first = i // 100 - 1
#        second = i % 100
#        if mugi[g, i].X > 0.5:
#            plt.plot(Rgi[g, i, 0].X, Rgi[g, i, 1].X, 'kD', markersize=1, color='red')
#            # ax.annotate("$R_" + str(g) + "^{" + str((first, second)) + "}$", xy = (Rgi[g, i, 0].X+1, Rgi[g, i, 1].X+1))
#            plt.plot(Lgi[g, i, 0].X, Lgi[g, i, 1].X, 'kD', markersize=1, color='red')
#
#
#    print(path_C)
#    def point_in_list(point, pos, lista):
#        return any([np.linalg.norm(np.array(point) - np.array(lista[u])) < 1e-4 for u in range(len(lista)-1) if u != i])
#
#    for p, i in zip(path_C, range(len(path_C))):
#        if point_in_list(p, i, path_C):
#            path_C.remove(p)
#
#    print(path_C)
#
#    for p, i in zip(path, range(len(path))):
#        # path_C.append([xL[t, 0].X, xL[t, 1].X])
#        # path_C.append([xR[t, 0].X, xR[t, 1].X])
#        if i == 0 or i == len(path)-1:
#            plt.plot(xL[p, 0].X, xL[p, 1].X, 'ko', markersize=5, color='orange')
#            plt.plot(xR[p, 0].X, xR[p, 1].X, 'ko', markersize=5, color='orange')
#        else:
#            plt.plot(xL[p, 0].X, xL[p, 1].X, 'ko', markersize=5, color='green')
#            # ax.annotate(str(i), xy = (xL[p, 0].X+0.5, xL[p, 1].X+0.5))
#            plt.plot(xR[p, 0].X, xR[p, 1].X, 'ko', markersize=5, color='blue')
#            # if np.linalg.norm(np.array([xL[p, 0].X, xL[p, 1].X]) - np.array([xR[p, 0].X, xR[p, 1].X])) > 1:
#            #     ax.annotate(str(i), xy = (xR[p, 0].X+0.5, xR[p, 1].X+0.5))
#
#    for p, i in zip(path_C, range(len(path_C))):
#        ax.annotate(str(i), xy = (p[0]+2.5, p[1]-0.5))
#
#    ax.add_artist(Polygon(path_C, fill=False, animated=False, linestyle='-', alpha=1, color='blue'))
#
#    for paths in paths_D:
#        ax.add_artist(Polygon(paths, fill=False, closed=False,
#                      animated=False, alpha=1, color='red'))
#    #
#    # ax.add_artist(Polygon(path_D, fill=False, animated=False,
#    #               linestyle='dotted', alpha=1, color='red'))
#
#    # if elipses:
#    #     for e in elipses:
#    #         ax.add_artist(e.artist)
#
#    for g in T_index_prima:
#        grafo = grafos[g-1]
#        centroide = np.mean(grafo.V, axis = 0)
#        nx.draw(grafo.G, grafo.pos, node_size=40,
#                node_color='black', alpha=1, edge_color='black')
#        # ax.annotate(g, xy = (centroide[0], centroide[1]))
#        # nx.draw_networkx_labels(grafo.G, grafo.pos, font_color = 'red', font_size=5)
#
#    plt.savefig('./imagenes/step' + str(mode) + '.png',dpi = 200)
#    # plt.show()

    # plt.show()
    if mode == 1:
        return MODEL.getAttr('x', xL), MODEL.getAttr('x', xR), MODEL.ObjVal

    if mode == 2:
        return MODEL.getAttr('x', xL), MODEL.getAttr('x', xR), MODEL.getAttr('x', ugi), MODEL.getAttr('x', vgi), MODEL.ObjVal, path#
