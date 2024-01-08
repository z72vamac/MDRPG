"""Tenemos un conjunto E de entornos y un conjunto de poligonales P de las que queremos recorrer un porcentaje alfa p . Buscamos
   un tour de mínima distancia que alterne poligonal-entorno y que visite todas las poligonales"""


# Incluimos primero los paquetes
import gurobipy as gp
from gurobipy import GRB
import numpy as np
from itertools import product
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
from data import *
from entorno import *
import copy
import estimacion_M as eM
import vns
from tsp import tsp
from auxiliar_functions import path2matrix

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

n = 2
# origin = [0, 0]
# dest = [0, 0]
vD = 2
vC = 1
# nE = 2
# np.random.seed(15)
np.random.seed(10)
# origin = np.random.uniform(0, 100, 2).tolist()

origin = [50, 50]

dest = origin
#
# datos1 = Data([], nE, 3, 1, None, False, True, 0)
# datos1.generar_muestra()
# E = datos1.data

nP = 10
datos2 = Data([], nP, 2, 5, None, False, True, 0)
datos2.generar_muestra()
P = datos2.data
path = tsp(P)
matriz = path2matrix(path)

# E_index = range(nE)
T_index = range(nP+2)
T_index_prima = range(1, nP+1)
T_index_primaprima = range(0, nP+1)
N = range(n)

# Creacion del modelo
MODEL = gp.Model('klmedian-continuous')

# Variables de los puntos por donde salen y se recogen de la poligonal
# xl es el punto por el que se sale a la poligonal
# xr es el punto por el que entra de la poligonal

x_index = []

for t in T_index:
    for dim in N:
        x_index.append((t, dim))

xl = MODEL.addVars(x_index, vtype=GRB.CONTINUOUS, name='x1')
xr = MODEL.addVars(x_index, vtype=GRB.CONTINUOUS, name='x2')

# Variables de asignacion de la poligonal p a la etapa t
y_index = []

for p in T_index_prima:
    for t in T_index_prima:
        y_index.append((p, t))

y = MODEL.addVars(y_index, vtype=GRB.BINARY, name='y1')

# Denotamos por dlr la distancia del punto de donde sale del camion al punto de la poligonal
dlr_index = []

for t in T_index:
    dlr_index.append(t)

dlr = MODEL.addVars(T_index_prima, vtype=GRB.CONTINUOUS, lb=0.0, name='dlr')
diflr = MODEL.addVars(T_index_prima, n, vtype=GRB.CONTINUOUS, lb=0.0, name='dlr')

drl = MODEL.addVars(T_index_primaprima, vtype=GRB.CONTINUOUS, lb=0.0, name='drl')
difrl = MODEL.addVars(T_index_primaprima, n, vtype=GRB.CONTINUOUS, lb=0.0, name='drl')

dlp = MODEL.addVars(T_index_prima, T_index_prima, vtype=GRB.CONTINUOUS, lb=0.0, name='difz')
diflp = MODEL.addVars(T_index_prima, T_index_prima,  2, vtype=GRB.CONTINUOUS, lb=0.0, name='difz')
dpr = MODEL.addVars(T_index_prima, T_index_prima, vtype=GRB.CONTINUOUS, lb=0.0, name='difz')
difpr = MODEL.addVars(T_index_prima, T_index_prima, 2,
                      vtype=GRB.CONTINUOUS, lb=0.0, name='difz')

plp = MODEL.addVars(T_index_prima, T_index_prima, vtype=GRB.CONTINUOUS,
                    lb=0.0, name='difz')
ppr = MODEL.addVars(T_index_prima, T_index_prima, vtype=GRB.CONTINUOUS,
                    lb=0.0, name='difz')

MODEL.update()

start = True
if start:
    for p, t in y_index:
        y[p, t].start = matriz[p-1, t-1]

MODEL.addConstrs((diflp[p, t, dim] >=  xl[t, dim] - P[p-1].V[dim]) for p, t, dim in diflp.keys())
MODEL.addConstrs((diflp[p, t, dim] >= -xl[t, dim] + P[p-1].V[dim]) for p, t, dim in diflp.keys())

MODEL.addConstrs(diflp[p, t, 0]*diflp[p, t, 0] + diflp[p, t, 1]*diflp[p, t, 1] <= dlp[p, t]*dlp[p, t] for p, t in dlp.keys())

SmallM = 0
BigM = max([np.linalg.norm(np.array(origin) - np.array(P[p].V)) for p in range(nP)])
BigM = 10000

MODEL.addConstrs(plp[p, t] >= SmallM * y[p, t] for p, t in plp.keys())
MODEL.addConstrs(plp[p, t] >= dlp[p, t] - BigM * (1 - y[p, t]) for p, t in plp.keys())

MODEL.addConstrs((difpr[p, t, dim] >=  P[p-1].V[dim] - xr[t, dim]) for p, t, dim in difpr.keys())
MODEL.addConstrs((difpr[p, t, dim] >= -P[p-1].V[dim] + xr[t, dim]) for p, t, dim in difpr.keys())

MODEL.addConstrs(difpr[p, t, 0]*difpr[p, t, 0] + difpr[p, t, 1]*difpr[p, t, 1] <= dpr[p, t]*dpr[p, t] for p, t in dpr.keys())

SmallM = 0
#BigM = 10000
MODEL.addConstrs(ppr[p, t] >= SmallM * y[p, t] for p, t in ppr.keys())
MODEL.addConstrs(ppr[p, t] >= dpr[p, t] - BigM * (1 - y[p, t]) for p, t in ppr.keys())

MODEL.addConstrs((difrl[t, dim] >=   xr[t, dim] - xl[t + 1, dim]) for t in T_index_primaprima for dim in range(2))
MODEL.addConstrs((difrl[t, dim] >=   -xr[t, dim] + xl[t + 1, dim]) for t in T_index_primaprima for dim in range(2))
MODEL.addConstrs(difrl[t, 0]*difrl[t, 0] + difrl[t, 1] * difrl[t, 1] <= drl[t] * drl[t] for t in T_index_primaprima)

MODEL.addConstrs((diflr[t, dim] >=   xl[t, dim] - xr[t, dim]) for t, dim in diflr.keys())
MODEL.addConstrs((diflr[t, dim] >= - xl[t, dim] + xr[t, dim]) for t, dim in diflr.keys())
MODEL.addConstrs(diflr[t, 0]*diflr[t, 0] + diflr[t, 1] * diflr[t, 1] <= dlr[t] * dlr[t] for t in dlr.keys())

MODEL.addConstrs((ppr[p, t] + plp[p, t])/vD <= dlr[t]/vC for p, t in ppr.keys())
MODEL.addConstrs((dlr[t]/vD <= 20) for t in dlr.keys())

MODEL.addConstrs(xl[0, dim] == origin[dim] for dim in N)
MODEL.addConstrs(xr[0, dim] == origin[dim] for dim in N)
#
MODEL.addConstrs(xl[nP+1, dim] == dest[dim] for dim in N)
MODEL.addConstrs(xr[nP+1, dim] == dest[dim] for dim in N)

# MODEL.addConstrs(v.sum('*', e, '*') == 1 for e in E_index)
# MODEL.addConstrs(z.sum(e, '*', '*') == 1 for e in E_index)
MODEL.addConstrs(y.sum(p, '*') == 1 for p in T_index_prima)
MODEL.addConstrs(y.sum('*', t) == 1 for t in T_index_prima)

MODEL.update()

# Funcion objetivo
# + gp.quicksum(0.5*plp[index] for index in plp.keys()) + gp.quicksum(0.5*ppr[index] for index in ppr.keys())
# objective = gp.quicksum(drl[index] for index in drl.keys()) + gp.quicksum(dlr[index] for index in dlr.keys()) + gp.quicksum(plp[index] for index in dlp.keys()) + gp.quicksum(ppr[index] for index in dpr.keys())

objective = gp.quicksum(dlr[index] for index in dlr.keys()) + gp.quicksum(drl[index] for index in drl.keys())

MODEL.setObjective(objective, GRB.MINIMIZE)

MODEL.update()

MODEL.setParam('TimeLimit', 600)
MODEL.write('aver.lp')

# Optimizamos
MODEL.optimize()
# MODEL.computeIIS()
# MODEL.write('infactible.ilp')

valsy = MODEL.getAttr('x', y)

selected_y = gp.tuplelist(e for e in valsy if valsy[e] > 0)
# print(selected_y)
#

# print(xl)
# print(xr)

fig, ax = plt.subplots()

min_x = []
max_x = []
min_y = []
max_y = []

# for p in T_index_prima:
#     dato = P[p-1]
    # min_x.append(min(P[0] for P in dato.V))
    # max_x.append(max(P[0] for P in dato.V))
    # min_y.append(min(P[1] for P in dato.V))
    # max_y.append(max(P[1] for P in dato.V))
    # plt.plot(dato.V)

#
# for e in E_index:
#     dato = E[e]
#     min_x.append(dato.centro[0] - dato.width)
#     max_x.append(dato.centro[0] + dato.width)
#     min_y.append(dato.centro[1] - dato.height)
#     max_y.append(dato.centro[1] + dato.height)
#     ax.add_artist(dato.artist)
#
#     ax.autoscale_view()
#     ax.axis([min(min_x) - 1, max(max_x) + 1,
#              min(min_y) - 1, max(max_y) + 1])
#     ax.set_aspect('equal')

# colores = []
# for t in range(2*nP + 4):
#     rgb = np.random.rand(3,)
#     colores.append(rgb)

# len(T_index)
path = []

# Representacion de los centros
# for e, p, t in z.keys():
#     if z[e, p, t].X == 1:
#         path.append([x1[t, 0].X, x1[t, 1].X])
#         path.append([a[p, 0, 0].X, a[p, 0, 1].X])
#         path.append([a[p, 1, 0].X, a[p, 1, 1].X])
#         path.append([x2[t, 0].X, x2[t, 1].X])

path_camion = []
for t in T_index:
    path_camion.append([xl[t, 0].X, xl[t, 1].X])
    path_camion.append([xr[t, 0].X, xr[t, 1].X])

paths = []

for t in T_index_prima:
    for p in T_index_prima:
        if y[p, t].X == 1:
            path_dron = []
            path_dron.append([xl[t, 0].X, xl[t, 1].X])
            path_dron.append(P[p-1].V)
            path_dron.append([xr[t, 0].X, xr[t, 1].X])
            ax.add_patch(Polygon(path_dron, fill=False, closed = False, linestyle='-', alpha=0.5, color='black'))
            for p in path_dron:
                plt.plot(p[0], p[1], 'ko', markersize=5, alpha=1, color='black')

for p in path_camion:
    plt.plot(p[0], p[1], 'ko', markersize=3, alpha=1, color='red')
#    ax.add_artist(Circle((p[0], p[1]), radius = 1.0))


# polygon_dron = Polygon(path_dron, fill=False, animated=True, linestyle='-', alpha=1, color = 'black')
# ax.add_patch(polygon_dron)

polygon_camion = Polygon(path_camion, fill=False, linestyle=':', alpha=1, color='red')
ax.add_patch(polygon_camion)


plt.title(str(nP) + "coordinate ")
# string = 'imagenes/' + str(m) + "-XPPN - MTZ - Mode " + \
#     str(modo) + " - Radii - " + str(datos.r) + ".png"
plt.savefig('poikonen.png')
#plt.show()
#plt.close()
plt.show()
