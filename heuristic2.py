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

np.random.seed(1)
orig = [50, 50]
dest = orig

nG = 10
datos = Data([], m=nG, r=2, modo=4, tmax=120, alpha = True,
             init=True,
             show=True,
             seed=2)
datos.generar_grafos()
grafos = datos.mostrar_datos()

T_index = range(datos.m + 2)
T_index_prima = range(1, datos.m+1)
T_index_primaprima = range(datos.m+1)

vD = 2

vC = 1

centroides = {}

for g in T_index_prima:
    centroides[g] = np.mean(grafos[g-1].V, axis = 0)

centros = []
centros.append(orig)
for g in centroides.values():
    centros.append(g)
centros.append(dest)

elipses = []

radio = 1

for c in centros:
    P = np.identity(2)
    q = -2*np.array(c)
    r = c[0]**2 + c[1]**2 - radio**2
    elipse = Elipse(P, q, r)
    elipses.append(elipse)

elipse = Data(elipses, m = 6, r = 2, modo  = 4, alpha = True, tmax = 1200, init = 0, prepro = 0, refor = 0, show = True, seed = 2)

path, path_P, obj  = MTZ(elipse)
print(path)
z = af.path2matrix(path)

xL, xR, obj = af.XPPNZ(datos, z, orig, dest, 0)


# print(xL)
xL_dict = {}

for g in T_index:
    xL_dict[g] = [xL[(g, 0)], xL[(g, 1)]]

# xL_dict = dict(sorted(xL_dict.items(), key=lambda x: x[0]))

radio = 0.5
elipses = []
for c in xL_dict.values():
    P = np.identity(2)
    q = -2*np.array(c)
    r = c[0]**2 + c[1]**2 - radio**2
    elipse = Elipse(P, q, r)
    elipses.append(elipse)

elipse = Data(elipses, m = 6, r = 2, modo  = 4, alpha = True, tmax = 1200, init = 0, prepro = 0, refor = 0, show = True, seed = 2)


objective = []
objective.append(10000)

for i in range(5):

    grafos = [grafos[a-1] for a in path[1:-1]]
    elipses = [elipses[a] for a in path]

    u_dict = {}
    zgij_dict = {}
    v_dict = {}
    xL_dict = {}
    obj_dict = {}

    path, path_P, obj = MTZ(elipse)

    z = af.path2matrix(path)

    # grafos = [grafos[a-1] for a in path[1:-1]]
    # elipses = [elipses[a] for a in path]


    # for g in range(1, nG+1):
    #     print('PROBLEMA del Grafo: ' + str(path[g]))
    #     vals_u, vals_zgij, vals_v, obj = af.XPPND(datos, elipses[g], grafos[g-1], elipses[g+1])
    #     for key, value in vals_u.items():
    #         u_dict[(g, key)] = value
    #     for key, value in vals_zgij.items():
    #         zgij_dict[(g, key[0], key[1])] = value
    #     for key, value in vals_v.items():
    #         v_dict[(g, key)] = value

        # xL_dict[0] = orig
        # xL_dict[g] = [vals_xL[0], vals_xL[1]]
        # xL_dict[g+1] = dest

    obj_dict[g] = obj

    print(u_dict)
    print(zgij_dict)

    # radio = 1
    # elipses = []
    # for p in path_P:
    #     P = np.identity(2)
    #     q = -2*np.array(c)
    #     r = c[0]**2 + c[1]**2 - radio**2
    #     elipse = Elipse(P, q, r)
    #     elipses.append(elipse)
    #
    # elipse = Data(elipses, m = 6, r = 2, modo = 4, alpha = True, tmax = 1200, init = 0, prepro = 0, refor = 0, show = True, seed = 2)

    datos = Data(grafos, m=nG, r=3, modo=4, tmax=120, alpha = True,
                 init=True,
                 show=True,
                 seed=2)

    # xL, xR, path_P, obj = af.XPPNZZ(datos, u_dict, zgij_dict, v_dict, z, orig, dest, i+1)
    xL, xR, obj = af.XPPNZ(datos, z, orig, dest, i+1, elipses)

    # print(xL)
    xL_dict = {}

    for g in T_index:
        xL_dict[g] = [xL[(g, 0)], xL[(g, 1)]]

    # xL_dict = dict(sorted(xL_dict.items(), key=lambda x: x[0]))

    radio = 0.5
    elipses = []
    for c in xL_dict.values():
        P = np.identity(2)
        q = -2*np.array(c)
        r = c[0]**2 + c[1]**2 - radio**2
        elipse = Elipse(P, q, r)
        elipses.append(elipse)

    elipse = Data(elipses, m = 6, r = 2, modo  = 4, alpha = True, tmax = 1200, init = 0, prepro = 0, refor = 0, show = True, seed = 2)

    # path, path_P, obj = MTZ(elipse)

#     xL_val = np.zeros((len(path)+1, 2))
#     xR_val = np.zeros((len(path)+1, 2))
#
#     #
#     for index in xL:
#         xL_val[index] = xL[index]
#         xR_val[index] = xR[index]
#
# #   xL_val = path_P
#
#     radio = 0.05
#     elipses = []
#     for c in xL_val:
#         P = np.identity(2)
#         q = -2*np.array(c)
#         r = c[0]**2 + c[1]**2 - radio**2
#         elipse = Elipse(P, q, r)
#         elipses.append(elipse)
#
#     elipse = Data(elipses, m = 6, r = 2, modo = 4, alpha = True, tmax = 1200, init = 0, prepro = 0, refor = 0, show = True, seed = 2)
    if obj >= min(objective):
        break
    else:
        objective.append(obj)


print(objective)
#
# fig, ax = plt.subplots()
# plt.axis([0, 100, 0, 100])
#
# for g in T_index_prima:
#     grafo = grafos[g-1]
#     nx.draw(grafo.G, grafo.pos, node_size=20,
#             node_color='black', alpha=0.3, edge_color='gray')
#     ax.annotate(g, xy = (centroides[g][0], centroides[g][1]))
#
#
# for c in centros:
#     plt.plot(c[0], c[1], 'ko', markersize=1, color='cyan')
#
# plt.show()
