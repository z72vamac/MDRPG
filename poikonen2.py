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
matriz = af.path2matrix(path)

# E_index = range(nE)
T_index = range(nP+2)
T_index_prima = range(1, nP+1)
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

xl = MODEL.addVars(x_index, vtype=GRB.CONTINUOUS, name='xL')
xr = MODEL.addVars(x_index, vtype=GRB.CONTINUOUS, name='xR')

# Variables de asignacion de la poligonal p a la etapa t
z_index = []

for v in T_index:
    for w in T_index:
        if v != w:
            z_index.append((v, w))

z = MODEL.addVars(z_index, vtype=GRB.BINARY, name='z')

# Denotamos por dlr la distancia del punto de donde sale del camion al punto de la poligonal
dlr = MODEL.addVars(T_index_prima, vtype=GRB.CONTINUOUS, lb=0.0, name='dlr')
diflr = MODEL.addVars(T_index_prima, n, vtype=GRB.CONTINUOUS, lb=0.0, name='diflr')

drl = MODEL.addVars(z_index, vtype=GRB.CONTINUOUS, lb=0.0, name='drl')
difrl = MODEL.addVars(z_index, n, vtype=GRB.CONTINUOUS, lb=0.0, name='difrl')

dlp = MODEL.addVars(T_index_prima, vtype=GRB.CONTINUOUS, lb=0.0, name='dlp')
diflp = MODEL.addVars(T_index_prima,  2, vtype=GRB.CONTINUOUS, lb=0.0, name='diflp')
dpr = MODEL.addVars(T_index_prima, vtype=GRB.CONTINUOUS, lb=0.0, name='dpr')
difpr = MODEL.addVars(T_index_prima, 2, vtype=GRB.CONTINUOUS, lb=0.0, name='difpr')

prl = MODEL.addVars(z_index, vtype=GRB.CONTINUOUS, lb=0.0, name='prl')

s = MODEL.addVars(T_index, vtype=GRB.CONTINUOUS, lb=0, name='sgi')

MODEL.update()

start = False
if start:
    for v, w in z_index:
        if w < 10 and v < 10:
            z[v, w].start = matriz[v, w]

# Distancia del dron al punto de servicio
MODEL.addConstrs((diflp[p, dim] >=  xl[p, dim] - P[p-1].V[dim]) for p, dim in diflp.keys())
MODEL.addConstrs((diflp[p, dim] >= -xl[p, dim] + P[p-1].V[dim]) for p, dim in diflp.keys())

MODEL.addConstrs(diflp[p, 0]*diflp[p, 0] + diflp[p, 1]*diflp[p, 1] <= dlp[p]*dlp[p] for p in dlp.keys())

SmallM = 0
BigM = max([np.linalg.norm(np.array(origin) - np.array(P[p].V)) for p in range(nP)])
#BigM = 10000

# Distancia del dron del punto de servicio al punto de recogida
MODEL.addConstrs((difpr[p, dim] >=  xr[p, dim] - P[p-1].V[dim]) for p, dim in difpr.keys())
MODEL.addConstrs((difpr[p, dim] >= -xr[p, dim] + P[p-1].V[dim]) for p, dim in difpr.keys())

MODEL.addConstrs(difpr[p, 0]*difpr[p, 0] + difpr[p, 1]*difpr[p, 1] <= dpr[p]*dpr[p] for p in dpr.keys())


# Distancia del camion del punto de lanzamiento al de recogida
MODEL.addConstrs((diflr[p, dim] >=   xl[p, dim] - xr[p, dim]) for p, dim in diflr.keys())
MODEL.addConstrs((diflr[p, dim] >= - xl[p, dim] + xr[p, dim]) for p, dim in diflr.keys())
MODEL.addConstrs(diflr[p, 0]*diflr[p, 0] + diflr[p, 1] * diflr[p, 1] <= dlr[p] * dlr[p] for p in dlr.keys())

# Distancia del camino del punto de recogida al punto de lanzamiento al siguiente punto de servicio
MODEL.addConstrs((difrl[v, w, dim] >=   xr[v, dim] - xl[w, dim]) for v, w, dim in difrl.keys())
MODEL.addConstrs((difrl[v, w, dim] >=  -xr[v, dim] + xl[w, dim]) for v, w, dim in difrl.keys())
MODEL.addConstrs((difrl[v, w, 0]*difrl[v, w, 0] + difrl[v, w, 1]*difrl[v, w, 1] <=  drl[v, w]*drl[v, w]) for v, w in drl.keys())


MODEL.addConstrs(prl[v, w] >= SmallM * z[v, w] for v, w in prl.keys())
MODEL.addConstrs(prl[v, w] >= drl[v, w] - BigM * (1 - z[v, w]) for v, w in prl.keys())



# SmallM = 0
#BigM = 10000


MODEL.addConstrs((dlp[p] + dpr[p])/vD <= dlr[p]/vC for p in dlr.keys())
MODEL.addConstrs((dlr[p]/vD <= 20) for p in dlp.keys())

MODEL.addConstrs(xl[0, dim] == origin[dim] for dim in N)
MODEL.addConstrs(xr[0, dim] == origin[dim] for dim in N)
#
MODEL.addConstrs(xl[nP+1, dim] == dest[dim] for dim in N)
MODEL.addConstrs(xr[nP+1, dim] == dest[dim] for dim in N)

# MODEL.addConstrs(v.sum('*', e, '*') == 1 for e in E_index)
# MODEL.addConstrs(z.sum(e, '*', '*') == 1 for e in E_index)

# MODEL.addConstr(z[0, 4] == 1)
# MODEL.addConstr(z[1, 9] == 1)
# MODEL.addConstr(z[2, 6] == 1)
# MODEL.addConstr(z[3, 1] == 1)
# MODEL.addConstr(z[4, 7] == 1)
# MODEL.addConstr(z[5, 3] == 1)
# MODEL.addConstr(z[6, 11] == 1)
# MODEL.addConstr(z[7, 5] == 1)
# MODEL.addConstr(z[8, 10] == 1)
# MODEL.addConstr(z[9, 8] == 1)
# MODEL.addConstr(z[10, 2] == 1)


MODEL.addConstr(gp.quicksum(z[v, 0] for v in T_index_prima) == 0)

MODEL.addConstr(gp.quicksum(z[nP+1, w] for w in T_index_prima) == 0)

MODEL.addConstrs(gp.quicksum(z[v , w] for w in T_index if w != v) == 1 for v in T_index_prima)
MODEL.addConstrs(gp.quicksum(z[w , v] for w in T_index if w != v) == 1 for v in T_index_prima)

# # Eliminación de subtours
for v in T_index_prima:
    for w in T_index_prima:
        if v != w:
            MODEL.addConstr(len(T_index) - 1 >= (s[v] - s[w]) + len(T_index) * z[v, w])
#
# # for g in range(nG):
# #     MODEL.addConstr(sgi[g, grafos[g].aristas[0]] == 0)
#
for v in T_index_prima:
    MODEL.addConstr(s[v] >= 0)
    MODEL.addConstr(s[v] <= len(T_index) - 1)

MODEL.addConstr(s[0] == 0)
MODEL.addConstr(s[nP + 1] == nP+1)

MODEL.update()

# Funcion objetivo
# + gp.quicksum(0.5*plp[index] for index in plp.keys()) + gp.quicksum(0.5*ppr[index] for index in ppr.keys())
# objective = gp.quicksum(drl[index] for index in drl.keys()) + gp.quicksum(dlr[index] for index in dlr.keys()) + gp.quicksum(plp[index] for index in dlp.keys()) + gp.quicksum(ppr[index] for index in dpr.keys())

objective = gp.quicksum(dlr[index] for index in dlr.keys()) + gp.quicksum(prl[index] for index in drl.keys()) + gp.quicksum(dlp[index] for index in dlp.keys()) + gp.quicksum(dpr[index] for index in dpr.keys())

MODEL.setObjective(objective, GRB.MINIMIZE)

MODEL.update()

MODEL.setParam('TimeLimit', 600)
MODEL.write('aver.lp')

# Optimizamos
MODEL.optimize()
# MODEL.computeIIS()
# MODEL.write('infactible.ilp')

valsy = MODEL.getAttr('x', z)

selected_y = gp.tuplelist(e for e in valsy if valsy[e] > 0)
print(selected_y)
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

path = af.subtour(selected_y)
print(path)

# path = path[0:-1]

# len(T_index)

# Representacion de los centros
# for e, p, t in z.keys():
#     if z[e, p, t].X == 1:
#         path.append([x1[t, 0].X, x1[t, 1].X])
#         path.append([a[p, 0, 0].X, a[p, 0, 1].X])
#         path.append([a[p, 1, 0].X, a[p, 1, 1].X])
#         path.append([x2[t, 0].X, x2[t, 1].X])

path_camion = []
for p in path:
    path_camion.append([xl[p, 0].X, xl[p, 1].X])
    path_camion.append([xr[p, 0].X, xr[p, 1].X])


# paths = []
#
for p in path:
    path_dron = []
    path_dron.append([xl[p, 0].X, xl[p, 1].X])
    path_dron.append(P[p-1].V)
    path_dron.append([xr[p, 0].X, xr[p, 1].X])
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
plt.savefig('poikonen2.png')
#plt.show()
#plt.close()
plt.show()
