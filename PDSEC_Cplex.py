"""Tenemos un conjunto E de entornos y un conjunto de poligonales P de las que queremos recorrer un porcentaje alfa p . Buscamos
   un tour de mínima distancia que alterne poligonal-entorno y que visite todas las poligonales"""


# Incluimos primero los paquetes
import cplex
import docplex.mp.model as cpx
from docplex.mp.model import Model
from docplex.cp.parameters import *
from docplex.util.environment import get_environment
import pdb
import numpy as np
from itertools import permutations, chain, combinations
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
import time

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

# np.random.seed(3)
#
# lista = [3, 3, 3, 3]
# nG = len(lista)
# datos = Data([], m=nG, r=1, grid = False, tmax=120, alpha = True,
#              init=True,
#              show=True,
#              seed=2)
#
# datos.generar_grid()
#
#
# datos.generar_grafos(lista)

def PDSEC_cplex(datos):

    origin= [50,50]
    dest = origin
    grafos = datos.mostrar_datos()

    # print(grafos[0].V)

    result = []

    nG = datos.m
    T_index = range(datos.m + 2)
    T_index_prima = range(1, datos.m+1)
    T_index_primaprima = range(datos.m+1)

    vD = 2

    vC = 1
    # Creamos el modelo8
    MODEL = Model(name="PD-Graph", log_output = True)
    MODEL.context.cplex_parameters.threads = 6

    # Variables que modelan las distancias
    # Variable binaria ugi = 1 si en la etapa t entramos por el segmento sgi
    ugi_index = []
    auxgLi_index = []

    for g in T_index_prima:
        for i in grafos[g-1].aristas:
            ugi_index.append((g, i))
            for dim in range(2):
                auxgLi_index.append((g, i, dim))


    ugi = MODEL.binary_var_dict(ugi_index, name='ugi')

    # Variable continua no negativa dgLi que indica la distancia desde el punto de lanzamiento hasta el segmento
    # sgi.
    dgLi_index = ugi_index

    dgLi = MODEL.continuous_var_dict(dgLi_index, lb=0.0, name='dgLi')
    auxgLi = MODEL.continuous_var_dict(auxgLi_index, lb=0.0, name='difgLi')

    # Variable continua no negativa pgLi = ugi * dgLi
    pgLi_index = ugi_index

    pgLi = MODEL.continuous_var_dict(pgLi_index, lb=0.0, name='pgLi')


    # Variable binaria vgi = 1 si en la etapa t salimos por el segmento sgi
    vgi_index = ugi_index

    vgi = MODEL.binary_var_dict(vgi_index, name='vgi')

    # Variable continua no negativa dgRi que indica la distancia desde el punto de salida del segmento sgi hasta el
    # punto de recogida del camion
    dgRi_index = ugi_index

    dgRi = MODEL.continuous_var_dict(dgRi_index, lb=0.0, name='dgRi')
    auxgRi = MODEL.continuous_var_dict(auxgLi_index, lb=0.0, name='difgRi')


    # Variable continua no negativa pgRi = vgi * dgRi
    pgRi_index = ugi_index

    pgRi = MODEL.continuous_var_dict(pgRi_index, lb=0.0, name='pgRi')


    # Variable binaria zgij = 1 si voy del segmento i al segmento j del grafo g.
    zgij_index = []
    sgi_index = []
    auxgij_index = []

    for g in T_index_prima:
        for i in grafos[g-1].aristas:
            sgi_index.append((g, i))
            for j in grafos[g-1].aristas:
                if i != j:
                    zgij_index.append((g, i, j))
                    for dim in range(2):
                        auxgij_index.append((g, i, j, dim))

    zgij = MODEL.binary_var_dict(zgij_index, name='zgij')
    sgi = MODEL.continuous_var_dict(sgi_index, lb=0, name='sgi')

    # Variable continua no negativa dgij que indica la distancia entre los segmentos i j en el grafo g.
    dgij_index = zgij_index

    dgij = MODEL.continuous_var_dict(dgij_index, lb=0.0, name='dgij')
    auxgij = MODEL.continuous_var_dict(auxgij_index, lb=0.0, name='dgij')

    # Variable continua no negativa pgij = zgij * dgij
    pgij_index = zgij_index

    pgij = MODEL.continuous_var_dict(pgij_index, lb=0.0, name='pgij')

    auxLR_index = []
    for i in T_index_prima:
        for dim in range(2):
            auxLR_index.append((i, dim))

    # Distancia del punto de lanzamiento al punto de recogida
    dLR = MODEL.continuous_var_dict(T_index_prima, lb = 0.0, name = 'dLR')
    auxLR = MODEL.continuous_var_dict(auxLR_index, lb = 0.0, name = 'difLR')

    # Variable binaria z que vale uno si se va del grafo g al grafo g'
    z_index = []
    auxRL_index = []

    for v in T_index:
        for w in T_index:
            if v != w:
                z_index.append((v, w))
                for dim in range(2):
                    auxRL_index.append((v, w, dim))


    z = MODEL.binary_var_dict(z_index, name='z')
    # s = MODEL.continuous_var_dict(T_index,    lb=0, name='s')

    dRL = MODEL.continuous_var_dict(z_index, lb = 0.0, name = 'dRL')
    auxRL = MODEL.continuous_var_dict(auxRL_index, lb = 0.0, name = 'difRL')
    pRL = MODEL.continuous_var_dict(z_index, lb = 0.0, name = 'pRL')

    # Variables que modelan los puntos de entrada o recogida
    # xL: punto de salida del dron del camion en la etapa t
    xL_index = []

    for g in T_index:
        for dim in range(2):
            xL_index.append((g, dim))

    xL = MODEL.continuous_var_dict(xL_index, name='xL')

    # xR: punto de recogida del dron del camion en la etapa t
    xR_index = []

    for t in T_index:
        for dim in range(2):
            xR_index.append((t, dim))

    xR = MODEL.continuous_var_dict(xR_index, name='xR')

    # Rgi: punto de recogida del dron para el segmento sgi
    Rgi_index = []
    rhogi_index = []

    for g in T_index_prima:
        for i in grafos[g-1].aristas:
            rhogi_index.append((g, i))
            for dim in range(2):
                Rgi_index.append((g, i, dim))

    Rgi = MODEL.continuous_var_dict(Rgi_index, name='Rgi')
    rhogi = MODEL.continuous_var_dict(rhogi_index, lb=0.0, ub=1.0, name='rhogi')

    # Lgi: punto de lanzamiento del dron del segmento sgi
    Lgi_index = Rgi_index
    landagi_index = rhogi_index

    Lgi = MODEL.continuous_var_dict(Lgi_index, name='Lgi')
    landagi = MODEL.continuous_var_dict(landagi_index, lb=0.0, ub=1.0, name='landagi')

    # Variables auxiliares para modelar el valor absoluto
    mingi = MODEL.continuous_var_dict(rhogi_index, lb=0.0, ub = 1.0, name='mingi')
    maxgi = MODEL.continuous_var_dict(rhogi_index, lb=0.0, ub = 1.0, name='maxgi')
    entrygi = MODEL.binary_var_dict(rhogi_index, name='entrygi')
    mugi = MODEL.continuous_var_dict(rhogi_index, name = 'mugi')
    pgi = MODEL.continuous_var_dict(rhogi_index, lb = 0.0, ub = 1.0, name = 'pgi')
    alphagi = MODEL.continuous_var_dict(rhogi_index, lb = 0.0, ub = 1.0, name = 'alphagi')



    # Para cada grafo g, existe un segmento i y una etapa t donde hay que recoger al dron
    for g in T_index_prima:
        MODEL.add_constraint((MODEL.sum(ugi[g, i] for i in grafos[g-1].aristas) == 1), ctname = 'entrag')
        MODEL.add_constraint((MODEL.sum(vgi[g, i] for i in grafos[g-1].aristas) == 1), ctname = 'saleg')
    # MODEL.add_constraint(ugi.sum('*', i, '*') == 1 for i in range(nG))
    # MODEL.add_constraint(vgi.sum('*', i, '*') == 1 for g in range(nG))

    # De todos los segmentos hay que salir y entrar menos de aquel que se toma como entrada al grafo y como salida del grafo
    for g, i, j in zgij.keys():
        MODEL.add_constraint((mugi[g, i] - ugi[g, i] == MODEL.sum(zgij[g, u, i] for u in grafos[g-1].aristas if u != i)), ctname = 'flujou')
        MODEL.add_constraint((mugi[g, i] - vgi[g, i] == MODEL.sum(zgij[g, i, u] for u in grafos[g-1].aristas if u != i)), ctname = 'flujov')

    for g, i in rhogi.keys():
        MODEL.add_constraint((pgi[g, i] >= mugi[g, i] + alphagi[g, i] - 1), ctname = 'pgi1')
        MODEL.add_constraint((pgi[g, i] <= mugi[g, i]), ctname = 'pgi2')
        MODEL.add_constraint((pgi[g, i] <= alphagi[g, i]), ctname = 'pgi3')

    # MODEL.add_constraint(ugi[0, 101, 0] == 0)
    # MODEL.add_constraint(ugi[0, 101, 1] == 0)


    # Eliminación de subtours
    for g in T_index_prima:
        for i in grafos[g-1].aristas[0:]:
            for j in grafos[g-1].aristas[0:]:
                if i != j:
                    MODEL.add_constraint(grafos[g-1].num_aristas - 1 >= (sgi[g, i] - sgi[g, j]) + grafos[g-1].num_aristas * zgij[g, i, j])

    # for g in range(nG):
    #     MODEL.add_constraint(sgi[g, grafos[g].aristas[0]] == 0)

    for g in T_index_prima:
        for i in grafos[g-1].aristas[0:]:
            MODEL.add_constraint(sgi[g, i] >= 0)
            MODEL.add_constraint(sgi[g, i] <= grafos[g-1].num_aristas - 1)


    # Restricciones de distancias y producto
    for g, i, dim in auxgLi.keys():
        MODEL.add_constraint(auxgLi[g, i, dim] >=   xL[g, dim] - Rgi[g, i, dim])
        MODEL.add_constraint(auxgLi[g, i, dim] >= - xL[g, dim] + Rgi[g, i, dim])

    for g, i in ugi.keys():
        MODEL.add_constraint((auxgLi[g, i, 0]*auxgLi[g, i, 0] + auxgLi[g, i, 1] * auxgLi[g, i, 1] <= dgLi[g, i] * dgLi[g, i]), ctname = 'u-conic')

    SmallM = 0
    BigM = 10000

    # BigM = 0
    # for g in T_index_prima:
    #     for v in grafos[g-1].V:
    #         BigM = max([np.linalg.norm(origin - v), BigM])
    #
    # BigM += 5
    #BigM = max([np.linalg.norm(origin-grafos[g].V) for g in range(nG)])

    # MODEL.add_constraint(ugi[1, 101] == 1)
    # MODEL.add_constraint(vgi[1, 203] == 1)

    for g, i in ugi.keys():
        MODEL.add_constraint(pgLi[g, i] >= SmallM * ugi[g, i])
        MODEL.add_constraint(pgLi[g, i] >= dgLi[g, i] - BigM * (1 - ugi[g, i]))

    for g, i, j, dim in auxgij.keys():
        MODEL.add_constraint(auxgij[g, i, j, dim] >=   Lgi[g, i, dim] - Rgi[g, j, dim])
        MODEL.add_constraint(auxgij[g, i, j, dim] >= - Lgi[g, i, dim] + Rgi[g, j, dim])

    for g, i, j in dgij.keys():
        MODEL.add_constraint((auxgij[g, i, j, 0]*auxgij[g, i, j, 0] + auxgij[g, i, j, 1] * auxgij[g, i, j, 1] <= dgij[g, i, j] * dgij[g, i, j]), ctname = 'zgij-conic')


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
        MODEL.add_constraint((pgij[g, i, j] >= SmallM_local * zgij[g, i, j]))
        MODEL.add_constraint((pgij[g, i, j] >= dgij[g, i, j] - BigM_local * (1 - zgij[g, i, j])))

    for g, i, dim in auxgRi.keys():
        MODEL.add_constraint(auxgRi[g, i, dim] >=   Lgi[g, i, dim] - xR[g, dim])
        MODEL.add_constraint(auxgRi[g, i, dim] >= - Lgi[g, i, dim] + xR[g, dim])

    for g, i in vgi.keys():
        MODEL.add_constraint((auxgRi[g, i, 0]*auxgRi[g, i, 0] + auxgRi[g, i, 1] * auxgRi[g, i, 1] <= dgRi[g, i] * dgRi[g, i]), ctname = 'v-conic')
        MODEL.add_constraint(pgRi[g, i] >= SmallM * vgi[g, i])
        MODEL.add_constraint(pgRi[g, i] >= dgRi[g, i] - BigM * (1 - vgi[g, i]))
    # SmallM = 0
    #BigM = 10000

    for g1, g2, dim in auxRL.keys():
        MODEL.add_constraint(auxRL[g1, g2, dim] >=   xR[g1, dim] - xL[g2, dim])
        MODEL.add_constraint(auxRL[g1, g2, dim] >= - xR[g1, dim] + xL[g2, dim])

    for g1, g2 in dRL.keys():
        MODEL.add_constraint((auxRL[g1, g2, 0]*auxRL[g1, g2, 0] + auxRL[g1, g2, 1] * auxRL[g1, g2, 1] <= dRL[g1, g2] * dRL[g1, g2]), ctname = 'RL-conic')
        MODEL.add_constraint(pRL[g1, g2] >= SmallM * z[g1, g2])
        MODEL.add_constraint(pRL[g1, g2] >= dRL[g1, g2] - BigM * (1 - z[g1, g2]))

    # Restricciones para formar un tour
    MODEL.add_constraint(MODEL.sum(z[v, 0] for v in T_index_prima) == 0)
    MODEL.add_constraint(MODEL.sum(z[nG+1, w] for w in T_index_prima) == 0)

    for v in T_index:
        MODEL.add_constraint(MODEL.sum(z[v , w] for w in T_index if w != v) == 1)
        MODEL.add_constraint(MODEL.sum(z[w , v] for w in T_index if w != v) == 1)

    def powerset(iterable):
        "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
        s = list(iterable)
        return chain.from_iterable(combinations(s, r) for r in range(2, len(s)+1))

    MODEL.add_lazy_constraints(MODEL.sum(z[v, w] for v, w in permutations(set, 2)) <= len(set) - 1 for set in list(powerset(T_index_prima)))
    # restricciones.Lazy = 3

    # MODEL.add_constraint(MODEL.sum(z[v, 0] for v in T_index if v != 0) == 0)
    # MODEL.add_constraint(MODEL.sum(z[nG+1, w] for w in T_index_prima) == 0)
    # MODEL.add_constraint(MODEL.sum(z[v , w] for w in T_index_prima if w != v) == 1 for v in T_index_prima)
    # MODEL.add_constraint(MODEL.sum(z[w , v] for w in T_index_prima if w != v) == 1 for v in T_index_prima)
    # MODEL.add_constraint(MODEL.sum(z[0, w] for w in T_index if w != 0) == 1)
    # MODEL.add_constraint(MODEL.sum(z[v, nG+1] for v in T_index if v != nG+1) == 0)

    # Conectividad
    # for v in T_index_prima:
    #     for w in T_index_prima:
    #         if v != w:
    #             MODEL.add_constraint(len(T_index) - 1 >= (s[v] - s[w]) + len(T_index) * z[v, w])
    # def powerset(iterable):
    #     "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    #     s = list(iterable)
    #     return chain.from_iterable(combinations(s, r) for r in range(2, len(s)+1))



    #restricciones = MODEL.add_constraint(MODEL.sum(z[v, w] for v, w in permutations(set, 2)) <= len(set) - 1 for set in list(powerset(T_index_prima)))
    #restricciones.Lazy = 3
    # for v in range(1, nG+1):
    #     MODEL.add_constraint(s[v] - s[0] + (nG+1 - 2)*z[0, v] <= len(T_index) - 1)
    #
    # for v in range(1, nG+1):
    #     MODEL.add_constraint(s[0] - s[v] + (nG+1 - 1)*z[v, 0] <= 0)

    # for v in range(1, nG+1):
    #     MODEL.add_constraint(-z[0,v] - s[v] + (nG+1-3)*z[v,0] <= -2, ctname="LiftedLB(%s)"%v)
    #     MODEL.add_constraint(-z[v,0] + s[v] + (nG+1-3)*z[0,v] <= nG+1-2, ctname="LiftedUB(%s)"%v)

    # for v in T_index_prima:
    #     MODEL.add_constraint(s[v] >= 1)
    #     MODEL.add_constraint(s[v] <= len(T_index) - 1)
    #
    # MODEL.add_constraint(s[0] == 0)
    # MODEL.add_constraint(s[nG + 1] == nG+1)

    for g, dim in auxLR.keys():
        MODEL.add_constraint(auxLR[g, dim] >=   xL[g, dim] - xR[g, dim])
        MODEL.add_constraint(auxLR[g, dim] >= - xL[g, dim] + xR[g, dim])

    for g in dLR.keys():
        MODEL.add_constraint((auxLR[g, 0]*auxLR[g, 0] + auxLR[g, 1] * auxLR[g, 1] <= dLR[g] * dLR[g]), ctname = 'LR-conic')

    for g in T_index_prima:
        MODEL.add_constraint((MODEL.sum(pgLi[g, i] for i in grafos[g-1].aristas)
                            + MODEL.sum(pgij[g, k, l] for k in grafos[g-1].aristas for l in grafos[g-1].aristas if k != l) + MODEL.sum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas)
                            + MODEL.sum(pgRi[g, i] for i in grafos[g-1].aristas))/vD <= dLR[g]/vC)

    # MODEL.add_constraint(dLR[g] <= 150 for g in dLR.keys())
    # MODEL.add_constraint((pgLi.sum('*', '*', t) +
    #                   pgij.sum(g, '*', '*') +
    #                   ugi.sum(g, '*', '*')*longitudes[g-1] +
    #                   pgRi.sum('*', '*', t))/vD <= dLR[t]/vC for t in T_index_prima for g in T_index_prima)
    # MODEL.add_constraint((dLR[t]/vD <= 50) for t in T_index_prima)
    # MODEL.add_constraint((pgLi[g, i, t]
    #                   + pgij.sum(g, '*', '*') + grafos[g-1].A[i // 100 - 1, i % 100]*grafos[g-1].longaristas[i // 100 - 1, i % 100]
    #                   + pgRi[g, i, t])/vD <= dLR[t]/vC for g, i, t in pgLi.keys())

    # MODEL.add_constraint(z[0, 2] + z[1, 3] + z[2, 1] + z[3, 4] == 4)

    for g, i in rhogi.keys():
        first = i // 100 - 1
        second = i % 100
        MODEL.add_constraint(rhogi[g, i] - landagi[g, i] == maxgi[g, i] - mingi[g, i])
        MODEL.add_constraint(maxgi[g, i] + mingi[g, i] == alphagi[g, i])
        if datos.alpha:
            MODEL.add_constraint(pgi[g, i] == grafos[g-1].A[first, second])
        MODEL.add_constraint(maxgi[g, i] <= 1 - entrygi[g, i])
        MODEL.add_constraint(mingi[g, i] <= entrygi[g, i])
        MODEL.add_constraint(Rgi[g, i, 0] == rhogi[g, i] * grafos[g-1].V[first, 0] + (1 - rhogi[g, i]) * grafos[g-1].V[second, 0])
        MODEL.add_constraint(Rgi[g, i, 1] == rhogi[g, i] * grafos[g-1].V[first, 1] + (1 - rhogi[g, i]) * grafos[g-1].V[second, 1])
        MODEL.add_constraint(Lgi[g, i, 0] == landagi[g, i] * grafos[g-1].V[first, 0] + (1 - landagi[g, i]) * grafos[g-1].V[second, 0])
        MODEL.add_constraint(Lgi[g, i, 1] == landagi[g, i] * grafos[g-1].V[first, 1] + (1 - landagi[g, i]) * grafos[g-1].V[second, 1])

    if not(datos.alpha):
        for g in T_index_prima:
            MODEL.add_constraint(MODEL.sum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas) >= grafos[g-1].alpha*grafos[g-1].longitud)


    # Origen y destino
    for dim in range(2):
        MODEL.add_constraint(xL[0, dim] == origin[dim])
        MODEL.add_constraint(xR[0, dim] == origin[dim])

        MODEL.add_constraint(xR[nG+1, dim] == dest[dim])
        MODEL.add_constraint(xL[nG+1, dim] == dest[dim])

    # MODEL.update()

    objective = MODEL.sum(pgLi[g, i] + pgRi[g, i] for g, i in pgRi.keys()) + MODEL.sum(pgij[g, i, j] for g, i, j in pgij.keys()) + MODEL.sum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for g in T_index_prima for i in grafos[g-1].aristas) + MODEL.sum(3*dLR[g] for g in dLR.keys()) + MODEL.sum(3*pRL[g1, g2] for g1, g2 in dRL.keys())

    #objective = MODEL.sum(1*dLR[g] for g in dLR.keys()) + MODEL.sum(1*pRL[g1, g2] for g1, g2 in dRL.keys())

    # objective = MODEL.sum(dRL[t] + dLR[t] for t in T_index)

    first_time = time.time()
    MODEL.minimize(objective)
    # MODEL.Params.Threads = 8
    # MODEL.Params.NonConvex = 2
    # MODEL.Params.timeLimit = datos.tmax

    # c = cplex.Cplex()
    MODEL.parameters.timelimit = datos.tmax
    # params.LogVerbosity = 'Verbose'
    # MODEL.update()

    # MODEL.write('AMDRPG-SEC.lp')
    # MODEL.write('AMDRPG-SEC.mps')
    MODEL.solve()

    second_time = time.time()
    z_val = []
    for c1, c2 in z.keys():
        if z[c1, c2].solution_value > 0.5:
            z_val.append((c1, c2))

    print(z_val)
    cycle = af.subtour_cplex(z_val)
    print(cycle)

    # iter = 1
    # print(MODEL.sum(z[c1, c2] for c1, c2 in permutations(cycle, 2)) <= len(cycle) - 1 )
    #
    # while second_time - first_time <= datos.tmax:
    #     if len(cycle) < nG + 2:
    #         print(cycle)
    #         print('Subtour found - Iteration: {u}'.format(u = iter))
    #         MODEL.add_lazy_constraint(MODEL.sum(z[c1, c2] for c1, c2 in permutations(cycle, 2)) <= len(cycle) - 1)
    #         MODEL.solve()
    #         second_time = time.time()
    #         iter += 1
    #     else:
    #         break

    # print(cycle)
    # print(MODEL.solve_details.time)
    # MODEL.update()

    # if MODEL.Status == 3:
    #     MODEL.computeIIS()
    #     MODEL.write('casa.ilp')
    #     result =  [np.nan, np.nan, np.nan, np.nan]
    #     if datos.grid:
    #         result.append('Grid')
    #     else:
    #         result.append('Delauney')
    #
    #     result.append('SEC')
    #
    #     return result
    #
    # if MODEL.SolCount == 0:
    #     result =  [np.nan, np.nan, np.nan, np.nan]
    #     if datos.grid:
    #         result.append('Grid')
    #     else:
    #         result.append('Delauney')
    #
    #     result.append('SEC')
    #
    #     return result

    result.append(MODEL.solve_details.mip_relative_gap)
    result.append(MODEL.solve_details.time)
    result.append(MODEL.solve_details.nb_nodes_processed)
    result.append(MODEL.objective_value)

    if datos.grid:
        result.append('Grid')
    else:
        result.append('Delauney')

    result.append('SEC')

    # print(result)
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
    # # print(path)
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
    # # path_C.append([xLt[nG+1, 0].X, xLt[nG+1, 1].X])
    # fig, ax = plt.subplots()
    # plt.axis([0, 100, 0, 100])
    #
    # for g, i in rhogi.keys():
    #     first = i // 100 - 1
    #     second = i % 100
    #     if mugi[g, i].X > 0.5:
    #         plt.plot(Rgi[g, i, 0].X, Rgi[g, i, 1].X, 'kD', markersize=1, color='red')
    #         ax.annotate("$R_" + str(g) + "^{" + str((first, second)) + "}$", xy = (Rgi[g, i, 0].X+1, Rgi[g, i, 1].X+1))
    #         plt.plot(Lgi[g, i, 0].X, Lgi[g, i, 1].X, 'kD', markersize=1, color='red')
    #         ax.annotate("$L_" + str(g) + "^{" + str((first, second)) + "}$", xy = (Lgi[g, i, 0].X-3, Lgi[g, i, 1].X-3))
    #
    # #
    # for p, i in zip(path, range(len(path))):
    #     # path_C.append([xL[t, 0].X, xL[t, 1].X])
    #     # path_C.append([xR[t, 0].X, xR[t, 1].X])
    #     plt.plot(xL[p, 0].X, xL[p, 1].X, 'ko', alpha = 0.3, markersize=10, color='green')
    #     ax.annotate("L" + str(i), xy = (xL[p, 0].X+0.5, xL[p, 1].X+0.5))
    #     plt.plot(xR[p, 0].X, xR[p, 1].X, 'ko', markersize=5, color='blue')
    #     ax.annotate("R" + str(i), xy = (xR[p, 0].X-1.5, xR[p, 1].X-1.5))
    #
    #
    # ax.add_artist(Polygon(path_C, fill=False, closed = False, animated=False,
    #               linestyle='-', alpha=1, color='blue'))
    #
    # for path in paths_D:
    #     ax.add_artist(Polygon(path, fill=False, closed=False, lw = 0.1,
    #                   animated=False, alpha=1, color='red'))
    # #
    # # ax.add_artist(Polygon(path_D, fill=False, animated=False,
    # #               linestyle='dotted', alpha=1, color='red'))
    #
    # for g in range(nG):
    #     grafo = grafos[g]
    #     nx.draw(grafo.G, grafo.pos, node_size=40,
    #             node_color='black', alpha=1, edge_color='black')
    #
    # plt.savefig('PDMTZ-' + str(result[4]) +  '.png')

    print(result)
    return result
    # plt.show()

# PDSEC_cplex(datos)
