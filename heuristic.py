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

np.random.seed(2)
orig = [50, 50]
dest = orig

nG = 10
datos = Data([], m=nG, r=3, modo=4, tmax=120, alpha = True,
             init=True,
             show=True,
             seed=2)
datos.generar_grafos()

grafos = datos.mostrar_datos()

objective = []
# Paso 1: Calculo el punto de minima distancia de cada grafo al origen  (Idea para crear otro tipo de clustering por grafos)
puntos = [af.dist_grafo(orig, grafo)[1] for grafo in grafos]

# Paso 2: Calculo el punto medio entre los puntos calculados y el origen
centros = []
centros.append(orig)
for p in puntos:
    centros.append([(orig[0] + p[0])/2, (orig[1] + p[1])/2])
centros.append(dest)

# Paso 3: Calculo un TSPN con bolas cuyos centros son los puntos calculados en el Paso 2
elipses = []

radio = 1

for c in centros:
    P = np.identity(2)
    q = -2*np.array(c)
    r = c[0]**2 + c[1]**2 - radio**2
    elipse = Elipse(P, q, r)
    elipses.append(elipse)


elipse = Data(elipses, m = 6, r = 2, modo = 4, alpha = True, tmax = 1200, init = 0, prepro = 0, refor = 0, show = True, seed = 2)

path, path_P, obj  = MTZ(elipse)
print(path)

z = af.path2matrix(path)
af.XPPNZ(datos, z, orig, dest, elipses, 1)


# Paso 4: Opcion A. Resolver un problema de TSP para cada grafo donde xL esta en el entorno asociado a g y xR esta en el entorno asociado al grafo del siguiente tour
u_dict = {}
zgij_dict = {}
v_dict = {}

for g in range(1, nG+1):
    print('PROBLEMA del Grafo: ' + str(path[g]))
    vals_u, vals_zgij, vals_v  = af.XPPND(datos, elipses[g], grafos[g-1], elipses[g+1])
    for key, value in vals_u.items():
        u_dict[(g, key)] = value
    for key, value in vals_zgij.items():
        zgij_dict[(g, key[0], key[1])] = value
    for key, value in vals_v.items():
        v_dict[(g, key)] = value

xL, xR, path, obj = af.XPPNM(datos, u_dict, zgij_dict, v_dict, orig, dest, elipses, -1)

for i in range(5):
    elipses = [elipses[a] for a in path]

    u_dict = {}
    zgij_dict = {}
    v_dict = {}
    obj_dict = {}

    for g in range(1, nG+1):
        print('PROBLEMA del Grafo: ' + str(path[g]))
        vals_u, vals_zgij, vals_v, obj = af.XPPND(datos, elipses[g-1], grafos[g-1], elipses[g])
        for key, value in vals_u.items():
            u_dict[(g, key)] = value
        for key, value in vals_zgij.items():
            zgij_dict[(g, key[0], key[1])] = value
        for key, value in vals_v.items():
            v_dict[(g, key)] = value

        obj_dict[g] = obj

    # obj_dict = dict(sorted(obj_dict.items(), key=lambda x: x[1], reverse=True))
    # objective.append(sum(list(obj_dict.values())))

    # path = []
    # path.append(0)
    # for i in list(obj_dict.keys()):
    #     path.append(i)
    # #path = list(obj_dict.keys())
    # path.append(nG+1)

    # for e in elipses:
    #     e.radio = 5
    # print(obj_dict)
    #
    xL, xR, path, obj = af.XPPNM(datos, u_dict, zgij_dict, v_dict, orig, dest, elipses, i)

    # path, path_P, obj = MTZ(elipse)

    xL_val = np.zeros((len(path)+1, 2))
    xR_val = np.zeros((len(path)+1, 2))

    objective.append(obj)
    #
    for index in xL:
        xL_val[index] = xL[index]
        xR_val[index] = xR[index]

#   xL_val = path_P

    elipses = []
    for c in xL_val:
        P = np.identity(2)
        q = -2*np.array(c)
        r = c[0]**2 + c[1]**2 - radio**2
        elipse = Elipse(P, q, r)
        elipses.append(elipse)

    elipse = Data(elipses, m = 6, r = 2, modo = 4, alpha = True, tmax = 1200, init = 0, prepro = 0, refor = 0, show = True, seed = 2)
print(objective)
# print(u_list)
# print(v_list)

# Paso 5: Fijando estas variables binarias calculo xL y calculo xR libremente y el tour
# nG = 6
# datos = Data([], m=nG, r=3, modo=4, tmax=120, alpha = True,
#              init=True,
#              show=True,
#              seed=2)
# datos.generar_grafos()

# grafos = datos.mostrar_datos()


# ind = 0
# path_C = []
# paths_D = []
#
# for p in path:
#     path_C.append([xL_val[p, 0], xL_val[p, 1]])
#     path_C.append([xR_val[p, 0], xR_val[p, 1]])
#
#
# for p in path[1:]:
#     #    if ind < nG:
#     if ind < nG:
#         path_D = []
#         path_D.append([xL_val[p, 0], xL_val[p, 1]])
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
#         path_D.append([xR_val[p, 0], xR_val[p, 1]])
#     paths_D.append(path_D)

# path_C.append([xLt[nG+1, 0].X, xLt[nG+1, 1].X])
#
#
# for g, i in rhogi.keys():
#     plt.plot(Rgi[g, i, 0].X, Rgi[g, i, 1].X, 'ko', markersize=1, color='cyan')
#     plt.plot(Lgi[g, i, 0].X, Lgi[g, i, 1].X, 'ko', markersize=1, color='cyan')
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
# for g in range(nG):
#     grafo = grafos[g]
#     nx.draw(grafo.G, grafo.pos, node_size=20,
#             node_color='black', alpha=0.3, edge_color='gray')
#
# plt.savefig('heuristic-alphaG.png')
#
# plt.show()

#XPPNM(datos, u_list, zgij_list, v_list)







# fig, ax = plt.subplots()
# ax.set_aspect('equal')
# plt.axis([0, 100, 0, 100])
#
# ind = 0
# for g in grafos:
#     nx.draw(g.G, g.pos, node_size=20, node_color='black', alpha=0.3, edge_color='gray')
#     ax.annotate(ind, xy = g.V[0])
#     ind += 1
#
# for p in puntos:
#     plt.plot(p[0], p[1], 'ko', alpha = 1, markersize=2, color='black')
#
# for q in centros:
#     plt.plot(q[0], q[1], 'ko', alpha = 1, markersize = 2, color = 'green')
#
# for elipse in elipses:
#     ax.add_artist(elipse.artist)
#
# plt.plot(orig[0], orig[1], 'ko', alpha = 1, color = 'black')
#
# for p in xL_val:
#     plt.plot(p[0], p[1], 'ko', alpha = 1, color = 'black')

# for p in xR_list:
#     plt.plot(p[0], p[1], 'ko', alpha = 1, color = 'black')


# polygon = Polygon(path_P, fill=False, linestyle='-', alpha=1)
# ax.add_artist(polygon)

# plt.show()
