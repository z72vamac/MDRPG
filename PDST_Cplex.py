"""Tenemos un conjunto E de entornos y un conjunto de poligonales P de las que queremos recorrer un porcentaje alfa p . Buscamos
   un tour de mínima distancia que alterne poligonal-entorno y que visite todas las poligonales"""


# Incluimos primero los paquetes
import pdb
import docplex.mp.model as cpx
from docplex.mp.model import Model
from docplex.util.environment import get_environment
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
#
# lista = [4, 4]
# nG = len(lista)
# datos = Data([], m=nG, r=1, grid = False, tmax=600, alpha = True,
#              init=True,
#              show=True,
#              seed=2)
#
# datos.generar_grid()
#
#
#
# datos.generar_grafos(lista)

def PDST_cplex(datos):


    origin = [50, 50]
    dest = origin

    grafos = datos.mostrar_datos()

    result = []

    T_index = range(datos.m + 2)
    T_index_prima = range(1, datos.m+1)
    T_index_primaprima = range(datos.m+1)


    vD = 2

    vC = 1

    # Creamos el modelo8
    MODEL = Model("PD-Stages-Cplex", log_output = True)
    MODEL.context.cplex_parameters.threads = 6

    # Variables que modelan las distancias
    # Variable binaria ugit = 1 si en la etapa t entramos por el segmento sgi
    ugit_index = []
    auxgLit_index = []

    for g in T_index_prima:
        for i in grafos[g-1].aristas:
            for t in T_index_prima:
                ugit_index.append((g, i, t))
                for dim in range(2):
                    auxgLit_index.append((g, i, t, dim))


    ugit = MODEL.binary_var_dict(ugit_index,  name='ugit')

    # Variable continua no negativa dgLit que indica la distancia desde el punto de lanzamiento hasta el segmento
    # sgi.
    dgLit_index = ugit_index

    dgLit = MODEL.continuous_var_dict(dgLit_index,  lb=0.0, name='dgLit')
    auxgLit = MODEL.continuous_var_dict(auxgLit_index,  lb=0.0, name='auxgLit')

    # Variable continua no negativa pgLit = ugit * dgLit
    pgLit_index = ugit_index

    pgLit = MODEL.continuous_var_dict(pgLit_index,  lb=0.0, name='pgLit')


    # Variable binaria vgit = 1 si en la etapa t salimos por el segmento sgi
    vgit_index = ugit_index

    vgit = MODEL.binary_var_dict(vgit_index,  name='vgit')

    # Variable continua no negativa dgRit que indica la distancia desde el punto de salida del segmento sgi hasta el
    # punto de recogida del camion
    dgRit_index = ugit_index

    dgRit = MODEL.continuous_var_dict(dgRit_index,  lb=0.0, name='dgRit')
    auxgRit = MODEL.continuous_var_dict(auxgLit_index,  lb=0.0, name='auxgRit')


    # Variable continua no negativa pgRit = vgit * dgRit
    pgRit_index = ugit_index

    pgRit = MODEL.continuous_var_dict(pgRit_index,  lb=0.0, name='pgRit')


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


    zgij = MODEL.binary_var_dict(zgij_index,  name='zgij')
    sgi = MODEL.continuous_var_dict(sgi_index,  lb=0, name='sgi')

    # Variable continua no negativa dgij que indica la distancia entre los segmentos i j en el grafo g.
    dgij_index = zgij_index

    dgij = MODEL.continuous_var_dict(dgij_index,  lb=0.0, name='dgij')
    auxgij = MODEL.continuous_var_dict(auxgij_index,  lb=0.0, name='auxgij')

    # Variable continua no negativa pgij = zgij * dgij
    pgij_index = zgij_index

    pgij = MODEL.continuous_var_dict(pgij_index,  lb=0.0, name='pgij')

    # Variable continua no negativa dRLt que indica la distancia entre el punto de recogida en la etapa t y el punto de
    # salida para la etapa t+1
    dRLt_index = []
    auxRLt_index = []

    for t in T_index_primaprima:
        dRLt_index.append(t)
        for dim in range(2):
            auxRLt_index.append((t, dim))

    dRLt = MODEL.continuous_var_dict(dRLt_index,  lb=0.0, name='dRLt')
    auxRLt = MODEL.continuous_var_dict(auxRLt_index,  lb=0.0, name='auxRLt')

    # Variable continua no negativa dLRt que indica la distancia que recorre el camión en la etapa t mientras el dron se mueve
    dLRt_index = []
    auxLRt_index = []

    for t in T_index:
        dLRt_index.append(t)
        for dim in range(2):
            auxLRt_index.append((t, dim))

    dLRt = MODEL.continuous_var_dict(dLRt_index,  lb=0.0, name='dLRt')
    auxLRt = MODEL.continuous_var_dict(auxLRt_index,  lb=0.0, name='auxLRt')

    # Variables que modelan los puntos de entrada o recogida
    # xLt: punto de salida del dron del camion en la etapa t
    xLt_index = []

    for t in T_index:
        for dim in range(2):
            xLt_index.append((t, dim))

    xLt = MODEL.continuous_var_dict(xLt_index,  name='xLt')

    # xRt: punto de recogida del dron del camion en la etapa t
    xRt_index = []

    for t in T_index:
        for dim in range(2):
            xRt_index.append((t, dim))

    xRt = MODEL.continuous_var_dict(xRt_index,  name='xRt')

    # Rgi: punto de recogida del dron para el segmento sgi
    Rgi_index = []
    rhogi_index = []

    for g in T_index_prima:
        for i in grafos[g-1].aristas:
            rhogi_index.append((g, i))
            for dim in range(2):
                Rgi_index.append((g, i, dim))

    Rgi = MODEL.continuous_var_dict(Rgi_index,  name='Rgi')
    rhogi = MODEL.continuous_var_dict(rhogi_index,
                          lb=0.0, ub=1.0, name='rhogi')

    # Lgi: punto de lanzamiento del dron del segmento sgi
    Lgi_index = Rgi_index
    landagi_index = rhogi_index

    Lgi = MODEL.continuous_var_dict(Lgi_index,  name='Lgi')
    landagi = MODEL.continuous_var_dict(
        landagi_index,  lb=0.0, ub=1.0, name='landagi')

    # Variables auxiliares para modelar el valor absoluto
    mingi = MODEL.continuous_var_dict(rhogi_index,  lb=0.0, ub = 1.0, name='mingi')
    maxgi = MODEL.continuous_var_dict(rhogi_index,  lb=0.0, ub = 1.0, name='maxgi')
    entrygi = MODEL.binary_var_dict(rhogi_index,  name='entrygi')
    mugi = MODEL.binary_var_dict(rhogi_index, name = 'mugi')
    pgi = MODEL.continuous_var_dict(rhogi_index, lb = 0.0, ub = 1.0, name = 'pgi')
    alphagi = MODEL.continuous_var_dict(rhogi_index, lb = 0.0, ub = 1.0, name = 'alphagi')



    # En cada etapa hay que visitar/salir un segmento de un grafo
    MODEL.add_constraints(MODEL.sum(ugit[g, i, t] for g, i in rhogi.keys()) == 1 for t in T_index_prima)
    MODEL.add_constraints(MODEL.sum(vgit[g, i, t] for g, i in rhogi.keys()) == 1 for t in T_index_prima)

    # # Para cada grafo g, existe un segmento i y una etapa t donde hay que recoger al dron
    MODEL.add_constraints(MODEL.sum(ugit[g, i, t] for i in grafos[g-1].aristas for t in T_index_prima) == 1 for g in T_index_prima)
    MODEL.add_constraints(MODEL.sum(vgit[g, i, t] for i in grafos[g-1].aristas for t in T_index_prima) == 1 for g in T_index_prima)

    # MODEL.add_constraints(ugit.sum('*', i, '*') == 1 for i in range(nG))
    # MODEL.add_constraints(vgit.sum('*', i, '*') == 1 for g in range(nG))

    # De todos los segmentos hay que salir y entrar menos de aquel que se toma como entrada al grafo y como salida del grafo
    MODEL.add_constraints(mugi[g, i] - MODEL.sum(ugit[g, i, t] for t in T_index_prima) == MODEL.sum(zgij[g, j, i] for j in grafos[g-1].aristas if j != i) for g in T_index_prima for i in grafos[g-1].aristas)

    MODEL.add_constraints(mugi[g, i] - MODEL.sum(vgit[g, i, t] for t in T_index_prima) == MODEL.sum(zgij[g, i, j] for j in grafos[g-1].aristas if j != i) for g in T_index_prima for i in grafos[g-1].aristas)

    MODEL.add_constraints(MODEL.sum(ugit[g, i, t] for i in grafos[g-1].aristas) - MODEL.sum(vgit[g, i, t] for i in grafos[g-1].aristas)
                     == 0 for t in T_index_prima for g in T_index_prima)

    # MODEL.add_constraints((ugit.sum(g, i, '*') + zgij.sum(g, '*', i))
    #                  - (vgit.sum(g, i, '*') + zgij.sum(g, i, '*')) == 0 for g in T_index_prima for i in grafos[g-1].aristas)

    MODEL.add_constraints(pgi[g, i] >= mugi[g, i] + alphagi[g, i] - 1 for g, i in rhogi.keys())
    MODEL.add_constraints(pgi[g, i] <= mugi[g, i] for g, i in rhogi.keys())
    MODEL.add_constraints(pgi[g, i] <= alphagi[g, i] for g, i in rhogi.keys())

    # MODEL.add_constraint(ugit[0, 101, 0] == 0)
    # MODEL.add_constraint(ugit[0, 101, 1] == 0)


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
    MODEL.add_constraints((auxgLit[g, i, t, dim] >=   xLt[t, dim] - Rgi[g, i, dim]) for g, i, t, dim in auxgLit.keys())
    MODEL.add_constraints((auxgLit[g, i, t, dim] >= - xLt[t, dim] + Rgi[g, i, dim]) for g, i, t, dim in auxgLit.keys())

    for g, i, t in ugit.keys():
        MODEL.add_constraint(auxgLit[g, i, t, 0]*auxgLit[g, i, t, 0] + auxgLit[g, i, t, 1] * auxgLit[g, i, t, 1] <= dgLit[g, i, t] * dgLit[g, i, t], ctname = 'auxgLit')

    SmallM = 0
    BigM = 10000
    #
    # BigM = 0
    # for g in T_index_prima:
    #     for v in grafos[g-1].V:
    #         BigM = max([np.linalg.norm(origin - v), BigM])
    #
    # BigM += 5
    #BigM = max([np.linalg.norm(origin-grafos[g].V) for g in range(nG)])
    MODEL.add_constraints((pgLit[g, i, t] >= SmallM * ugit[g, i, t]) for g, i, t in ugit.keys())
    MODEL.add_constraints((pgLit[g, i, t] >= dgLit[g, i, t] - BigM * (1 - ugit[g, i, t])) for g, i, t in ugit.keys())

    MODEL.add_constraints((auxgij[g, i, j, dim] >=   Lgi[g, i, dim] - Rgi[g, j, dim]) for g, i, j, dim in auxgij.keys())
    MODEL.add_constraints((auxgij[g, i, j, dim] >= - Lgi[g, i, dim] + Rgi[g, j, dim]) for g, i, j, dim in auxgij.keys())

    for g, i, j in dgij.keys():
        MODEL.add_constraint(auxgij[g, i, j, 0]*auxgij[g, i, j, 0] + auxgij[g, i, j, 1] * auxgij[g, i, j, 1] <= dgij[g, i, j] * dgij[g, i, j], ctname = 'auxgij')


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

    MODEL.add_constraints((auxgRit[g, i, t, dim] >=   Lgi[g, i, dim] - xRt[t, dim]) for g, i, t, dim in auxgRit.keys())
    MODEL.add_constraints((auxgRit[g, i, t, dim] >= - Lgi[g, i, dim] + xRt[t, dim]) for g, i, t, dim in auxgRit.keys())

    for g, i, t in vgit.keys():
        MODEL.add_constraint(auxgRit[g, i, t, 0]*auxgRit[g, i, t, 0] + auxgRit[g, i, t, 1] * auxgRit[g, i, t, 1] <= dgRit[g, i, t] * dgRit[g, i, t], ctname = "auxgRit")


    # SmallM = 0
    #BigM = 10000
    MODEL.add_constraints((pgRit[g, i, t] >= SmallM * vgit[g, i, t]) for g, i, t in vgit.keys())
    MODEL.add_constraints((pgRit[g, i, t] >= dgRit[g, i, t] - BigM * (1 - vgit[g, i, t])) for g, i, t in vgit.keys())

    MODEL.add_constraints((auxRLt[t, dim] >=   xRt[t, dim] - xLt[t + 1, dim]) for t in dRLt.keys() for dim in range(2))
    MODEL.add_constraints((auxRLt[t, dim] >= - xRt[t, dim] + xLt[t + 1, dim]) for t in range(datos.m+1) for dim in range(2))

    for t in range(datos.m+1):
        MODEL.add_constraint(auxRLt[t, 0]*auxRLt[t, 0] + auxRLt[t, 1] * auxRLt[t, 1] <= dRLt[t] * dRLt[t], ctname = 'auxRLt')

    for t, dim in auxLRt.keys():
        MODEL.add_constraint(auxLRt[t, dim] >=   xLt[t, dim] - xRt[t, dim], ctname = 'error')
        MODEL.add_constraint(auxLRt[t, dim] >= - xLt[t, dim] + xRt[t, dim], ctname = 'error2')

    for t in dLRt.keys():
        MODEL.add_constraint(auxLRt[t, 0]*auxLRt[t, 0] + auxLRt[t, 1] * auxLRt[t, 1] <= dLRt[t] * dLRt[t], ctname = 'auxLRt')


    longitudes = []
    for g in T_index_prima:
        longitudes.append(sum([grafos[g-1].A[i // 100 - 1, i % 100]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas]))

    # BigM = 10000
    BigM = 10000000
    MODEL.add_constraints((MODEL.sum(pgLit[g, i, t] for i in grafos[g-1].aristas)
                      + MODEL.sum(pgij[g, i, j] for i in grafos[g-1].aristas for j in grafos[g-1].aristas if i != j) + MODEL.sum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for i in grafos[g-1].aristas)
                      + MODEL.sum(pgRit[g, i, t] for i in grafos[g-1].aristas))/vD <= dLRt[t]/vC + BigM*(1- MODEL.sum(ugit[g, i, t] for i in grafos[g-1].aristas)) for t in T_index_prima for g in T_index_prima)


    # MODEL.add_constraints((pgLit.sum('*', '*', t) +
    #                   pgij.sum(g, '*', '*') +
    #                   ugit.sum(g, '*', '*')*longitudes[g-1] +
    #                   pgRit.sum('*', '*', t))/vD <= dLRt[t]/vC for t in T_index_prima for g in T_index_prima)
    # MODEL.add_constraints((dLRt[t]/vD <= 50) for t in T_index_prima)
    # MODEL.add_constraints((pgLit[g, i, t]
    #                   + pgij.sum(g, '*', '*') + grafos[g-1].A[i // 100 - 1, i % 100]*grafos[g-1].longaristas[i // 100 - 1, i % 100]
    #                   + pgRit[g, i, t])/vD <= dLRt[t]/vC for g, i, t in pgLit.keys())

    for g, i in rhogi.keys():
        first = i // 100 - 1
        second = i % 100
        MODEL.add_constraint(rhogi[g, i] - landagi[g, i] == maxgi[g, i] - mingi[g, i])
        MODEL.add_constraint(maxgi[g, i] + mingi[g, i] == alphagi[g, i])
        if datos.alpha:
            MODEL.add_constraint(pgi[g, i] >= grafos[g-1].A[first, second])
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
    MODEL.add_constraints(xLt[0, dim] == origin[dim] for dim in range(2))
    MODEL.add_constraints(xRt[0, dim] == origin[dim] for dim in range(2))

    MODEL.add_constraints(xRt[datos.m+1, dim] == dest[dim] for dim in range(2))
    MODEL.add_constraints(xLt[datos.m+1, dim] == dest[dim] for dim in range(2))



    objective = MODEL.sum(pgLit[g, i, t] + pgRit[g, i, t] for g, i, t in pgRit.keys()) + MODEL.sum(pgij[g, i, j] for g, i, j in pgij.keys()) + MODEL.sum(pgi[g, i]*grafos[g-1].longaristas[i // 100 - 1, i % 100] for g in T_index_prima for i in grafos[g-1].aristas) + MODEL.sum(3*dLRt[t] for t in dLRt.keys()) + MODEL.sum(3*dRLt[t] for t in dRLt.keys())

    # objective = MODEL.sum(1*dLRt[g] for g in dLRt.keys()) + MODEL.sum(1*dRLt[g1] for g1 in dRLt.keys())

    # objective = MODEL.sum(dRLt[t] + dLRt[t] for t in T_index)

    MODEL.minimize(objective)

    MODEL.export_as_lp('PDST_Cplex.lp')

    MODEL.parameters.timelimit = datos.tmax
    # MODEL.Params.NonConvex = 2
    # MODEL.Params.timeLimit = 150



    # MODEL.write('AMDRPG-Stages.lp')
    # MODEL.write('AMDRPG-Stages.mps')
    MODEL.solve()


    # if MODEL.Status == 3:
    #     MODEL.computeIIS()
    #     MODEL.write('casa.ilp')
    #     result =  [np.nan, np.nan, np.nan, np.nan]
    #     if datos.grid:
    #         result.append('Grid')
    #     else:
    #         result.append('Delauney')
    #
    #     result.append('Stages')
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
    #     result.append('Stages')
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

    result.append('Stages')

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
    #
    #
    # # ind = 0
    # # path_C = []
    # # paths_D = []
    # #
    # # #path_C.append(origin)
    # # path_C.append([xLt[0, 0].X, xLt[0, 1].X])
    # # for t in T_index_prima:
    # #     #    if ind < nG:
    # #     path_C.append([xLt[t, 0].X, xLt[t, 1].X])
    # #     if ind < nG:
    # #         path_D = []
    # #         path_D.append([xLt[t, 0].X, xLt[t, 1].X])
    # #         index_g = 0
    # #         index_i = 0
    # #         for g, i, ti in selected_u:
    # #             if ti == t:
    # #                 index_g = g
    # #                 index_i = i
    # #         count = 0
    # #         while count < grafos[index_g-1].num_aristas-1:
    # #             for g, i, j in selected_z:
    # #                 if index_g == g and index_i == i:
    # #                     path_D.append([Rgi[g, i, 0].X, Rgi[g, i, 1].X])
    # #                     path_D.append([Lgi[g, i, 0].X, Lgi[g, i, 1].X])
    # #                     count += 1
    # #                     index_i = j
    # #         path_D.append([Rgi[index_g, index_i, 0].X, Rgi[index_g, index_i, 1].X])
    # #         path_D.append([Lgi[index_g, index_i, 0].X, Lgi[index_g, index_i, 1].X])
    # #         ind += 1
    # #         path_D.append([xRt[t, 0].X, xRt[t, 1].X])
    # #     paths_D.append(path_D)
    # #     path_C.append([xRt[t, 0].X, xRt[t, 1].X])
    #
    # path_C.append([xLt[datos.m+1, 0].X, xLt[datos.m+1, 1].X])
    #
    # fig, ax = plt.subplots()
    # plt.axis([0, 100, 0, 100])
    #
    # for g, i in rhogi.keys():
    #     plt.plot(Rgi[g, i, 0].X, Rgi[g, i, 1].X, 'ko', markersize=3, color='cyan')
    #     plt.plot(Lgi[g, i, 0].X, Lgi[g, i, 1].X, 'ko', markersize=3, color='cyan')
    # #
    # # path_C = []
    # for t in T_index:
    #     # path_C.append([xLt[t, 0].X, xLt[t, 1].X])
    #     # path_C.append([xRt[t, 0].X, xRt[t, 1].X])
    #     plt.plot(xLt[t, 0].X, xLt[t, 1].X, 'ko', markersize=5, color='green')
    #     ax.annotate("L" + str(t), xy = (xLt[t, 0].X, xLt[t, 1].X))
    #     plt.plot(xRt[t, 0].X, xRt[t, 1].X, 'ko', markersize=5, color='blue')
    #     ax.annotate("R" + str(t), xy = (xRt[t, 0].X, xRt[t, 1].X))
    #
    # ax.add_artist(Polygon(path_C, fill=False, animated=False,
    #               linestyle='-', alpha=1, color='blue'))
    #
    # for path in paths_D:
    #     ax.add_artist(Polygon(path, fill=False, closed=False,
    #                   animated=False, alpha=1, color='red'))
    # #
    # # ax.add_artist(Polygon(path_D, fill=False, animated=False,
    # #               linestyle='dotted', alpha=1, color='red'))
    #
    # for g in range(datos.m):
    #     grafo = grafos[g]
    #     nx.draw(grafo.G, grafo.pos, node_size=20,
    #             node_color='black', alpha=0.3, edge_color='gray')
    #
    # plt.savefig('PD-alphaG.png')

    print(result)
    return result
    # plt.show()

# PDST_cplex(datos)
