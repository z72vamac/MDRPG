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
origin = [0, 0]
dest = [0, 0]
vD = 30
vC = 10
# nE = 2
# np.random.seed(15)
np.random.seed(10)
#
# datos1 = Data([], nE, 3, 1, None, False, True, 0)
# datos1.generar_muestra()
# E = datos1.data

nP = 3
datos2 = Data([], nP, 2, 3, None, False, True, 0)
datos2.generar_muestra()
P = datos2.data

# E_index = range(nE)
P_index = range(nP)
T_index = range(0, nP+1)
T_index_prima = range(0, nP)
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

for p in P_index:
    for t in T_index:
        y_index.append((p, t))

y = MODEL.addVars(y_index, vtype=GRB.BINARY, name='y1')

# Denotamos por dlr la distancia del punto de donde sale del camion al punto de la poligonal
dlr_index = []

for t in T_index:
    dlr_index.append(t)

dlr = MODEL.addVars(T_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dlr')
diflr = MODEL.addVars(T_index, n, vtype=GRB.CONTINUOUS, lb=0.0, name='dlr')

drl = MODEL.addVars(T_index, vtype=GRB.CONTINUOUS, lb=0.0, name='drl')
difrl = MODEL.addVars(T_index, n, vtype=GRB.CONTINUOUS, lb=0.0, name='drl')

dlp = MODEL.addVars(T_index, P_index, vtype=GRB.CONTINUOUS,
                    lb=0.0, name='difz')
diflp = MODEL.addVars(T_index, P_index,  2,
                      vtype=GRB.CONTINUOUS, lb=0.0, name='difz')
dpr = MODEL.addVars(P_index, T_index, vtype=GRB.CONTINUOUS,
                    lb=0.0, name='difz')
difpr = MODEL.addVars(P_index, T_index, 2,
                      vtype=GRB.CONTINUOUS, lb=0.0, name='difz')

plp = MODEL.addVars(T_index, P_index, vtype=GRB.CONTINUOUS,
                    lb=0.0, name='difz')
ppr = MODEL.addVars(P_index, T_index, vtype=GRB.CONTINUOUS,
                    lb=0.0, name='difz')

# Variables que representan los puntos de la poligonal por donde se entra y se sale
a_index = []

for p in P_index:
    for dim in range(n):
        a_index.append((p, dim))

a = MODEL.addVars(a_index, 2, vtype=GRB.CONTINUOUS, name='a')

# Variables que parametrizan la poligonal
mu_index = []
landa_index = []
sublanda_index = []
u_index = []
s_index = []

for p in P_index:
    comp = P[p]
    u_index.append(p)
    # landa de la variable de entrada en la poligonal c
    landa_index.append((p, 0))
    # landa de la variable de salida en la poligonal c
    landa_index.append((p, 1))
    for segm in range(comp.num_segmentos):
        s_index.append((p, 0, segm))
        s_index.append((p, 1, segm))
    for punto in range(comp.num_puntos):
        sublanda_index.append((p, 0, punto))
        sublanda_index.append((p, 1, punto))

landa = MODEL.addVars(landa_index, vtype=GRB.CONTINUOUS, name='landa')
sublanda = MODEL.addVars(
    sublanda_index, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='sublanda')
s = MODEL.addVars(s_index, vtype=GRB.BINARY, name='s')
u = MODEL.addVars(u_index, vtype=GRB.BINARY, name='u')
minc = MODEL.addVars(u_index, vtype=GRB.CONTINUOUS, lb=0.0, name='min')
maxc = MODEL.addVars(u_index, vtype=GRB.CONTINUOUS, lb=0.0, name='max')

MODEL.update()

# Linearizacion del producto
for p, t in product(P_index, T_index_prima):
    MODEL.addConstrs(diflr[t, dim] >= xl[t+1, dim] - xr[t, dim] for dim in range(n))
    MODEL.addConstrs(diflr[t, dim] >= -xl[t+1, dim] + xr[t, dim] for dim in range(n))
    MODEL.addConstr(diflr[t, 0]*diflr[t, 0] + diflr[t, 1] * diflr[t, 1] <= dlr[t] * dlr[t])

    MODEL.addConstrs(difrl[t, dim] >= xl[t, dim] - xr[t, dim] for dim in range(n))
    MODEL.addConstrs(difrl[t, dim] >= -xl[t, dim] + xr[t, dim] for dim in range(n))
    MODEL.addConstr(difrl[t, 0]*difrl[t, 0] + difrl[t, 1] * difrl[t, 1] <= drl[t] * drl[t])

    MODEL.addConstrs(diflp[t, p, dim] >= xl[t, dim] - a[p, u, dim] for dim in range(n) for u in range(2))
    MODEL.addConstrs(diflp[t, p, dim] >= -xl[t, dim] + a[p, u, dim] for dim in range(n) for u in range(2))
    MODEL.addConstr(diflp[t, p, 0]*diflp[t, p, 0] + diflp[t, p, 1] * diflp[t, p, 1] <= dlp[t, p]*dlp[t, p])

    MODEL.addConstrs(difpr[p, t, dim] >= xr[t, dim] - a[p, u, dim] for dim in range(n) for u in range(2))
    MODEL.addConstrs(difpr[p, t, dim] >= -xr[t, dim] + a[p, u, dim] for dim in range(n) for u in range(2))
    MODEL.addConstr(difpr[p, t, 0]*difpr[p, t, 0] + difpr[p, t, 1] * difpr[p, t, 1] <= dpr[p, t]*dpr[p, t])

    SmallM = 0
    BigM = 10000
    comp = P[p]
    MODEL.addConstr(plp[t, p] >= SmallM * y[p, t])
    MODEL.addConstr(plp[t, p] >= dlp[t, p] - BigM * (1 - y[p, t]), name='liny1-const-2')
    MODEL.addConstr(ppr[p, t] >= SmallM * y[p, t])
    MODEL.addConstr(ppr[p, t] >= dpr[p, t] - BigM * (1 - y[p, t]), name='liny1-const-2')

    MODEL.addConstr((ppr[p, t] + plp[t, p] + comp.alpha * comp.longitud)/vD <= drl[t]/vC)
    MODEL.addConstr((drl[t]/vD <= 2))

for p in P_index:
    comp = P[p]
    if type(comp) is Poligonal:
        MODEL.addConstr(landa[p, 0] - landa[p, 1] == maxc[p] - minc[p], name='u0')
        # si u = 0, entonces landa0 >= landa1
        MODEL.addConstr(maxc[p] + minc[p] == comp.alpha * comp.num_segmentos, name='u1')
        MODEL.addConstr(maxc[p] <= comp.num_segmentos * (1 - u[p]), name='u2')
        MODEL.addConstr(minc[p] <= comp.num_segmentos * u[p], name='u3')
        for j in range(2):
            for punto in range(1, comp.num_puntos):
                MODEL.addConstr(landa[p, j] - punto >= sublanda[p, j, punto] - comp.num_puntos * (1 - s[p, j, punto - 1]))
                MODEL.addConstr(landa[p, j] - punto <= sublanda[p, j, punto] + comp.num_puntos * (1 - s[p, j, punto - 1]))
            MODEL.addConstr(sublanda[p, j, 0] <= s[p, j, 0])
            MODEL.addConstr(sublanda[p, j, comp.num_puntos - 1] <= s[p, j, comp.num_puntos - 2])
            for punto in range(1, comp.num_puntos - 1):
                MODEL.addConstr(sublanda[p, j, punto] <= s[p, j, punto - 1] + s[p, j, punto])
            MODEL.addConstr(s.sum(p, j, '*') == 1)
            MODEL.addConstr(sublanda.sum(p, j, '*') == 1)
            for dim in range(2):
                MODEL.addConstr(a[p, j, dim] == gp.quicksum(sublanda[p, j, punto] * comp.V[punto][dim] for punto in range(comp.num_puntos)), name='seg1')


# MODEL.addConstrs(xl[0, dim] == origin[dim] for dim in N)
# MODEL.addConstrs(xr[0, dim] == origin[dim] for dim in N)
#
# MODEL.addConstrs(xl[nP, dim] == dest[dim] for dim in N)
# MODEL.addConstrs(xr[nP, dim] == dest[dim] for dim in N)

# MODEL.addConstrs(v.sum('*', e, '*') == 1 for e in E_index)
# MODEL.addConstrs(z.sum(e, '*', '*') == 1 for e in E_index)
MODEL.addConstrs(y.sum(p, '*') == 1 for p in P_index)
MODEL.addConstrs(y.sum('*', t) == 1 for t in T_index_prima)

MODEL.update()

# Funcion objetivo
# + gp.quicksum(0.5*plp[index] for index in plp.keys()) + gp.quicksum(0.5*ppr[index] for index in ppr.keys())
objective = gp.quicksum(drl[index] for index in drl.keys()) + gp.quicksum(dlr[index] for index in dlr.keys())

MODEL.setObjective(objective, GRB.MINIMIZE)

MODEL.update()

MODEL.setParam('TimeLimit', 600)
MODEL.write('aver.lp')

# Optimizamos
MODEL.optimize()
# MODEL.computeIIS()
# MODEL.write('infactible.ilp')

valsy = MODEL.getAttr('x', y)
# valsz = MODEL.getAttr('x', z)
# valsv = MODEL.getAttr('x', v)
#
selected_y = gp.tuplelist(e for e in valsy if valsy[e] > 0)
# selected_z = gp.tuplelist(e for e in valsz if valsz[e] > 0)
# selected_v = gp.tuplelist(e for e in valsv if valsv[e] > 0)
#
# print(selected_z)
print(selected_y)
#
fig, ax = plt.subplots()

min_x = []
max_x = []
min_y = []
max_y = []

for p in P_index:
    dato = P[p]
    # min_x.append(min(P[0] for P in dato.V))
    # max_x.append(max(P[0] for P in dato.V))
    # min_y.append(min(P[1] for P in dato.V))
    # max_y.append(max(P[1] for P in dato.V))
    ax.add_artist(dato.artist)


vals_landa = MODEL.getAttr('x', landa)
vals_landa

vals_u = MODEL.getAttr('x', u)
vals_u
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
path_camion.append(origin)
for t in T_index:
    path_camion.append([xl[t, 0].X, xl[t, 1].X])
    path_camion.append([xr[t, 0].X, xr[t, 1].X])
path_camion.append(dest)

paths = []

for t in T_index:
    for p in P_index:
        if y[p, t].X == 1:
            path_dron = []
            path_dron.append([xl[t, 0].X, xl[t, 1].X])
            if u[p].X > 0.5:
                path_dron.append([a[p, 0, 0].X, a[p, 0, 1].X])
                for i in range(P[p].num_puntos):
                    if i >= round(vals_landa[(p, 0)], 3) and i <= round(vals_landa[(p, 1)]-1, 3):
                        print(i)
                        path_dron.append(P[p].V[i-1])
                path_dron.append([a[p, 1, 0].X, a[p, 1, 1].X])
            else:
                path_dron.append([a[p, 1, 0].X, a[p, 1, 1].X])
                for i in np.arange(0, P[p].num_puntos):
                    if i <= round(vals_landa[(p, 0)], 3) and i >= round(vals_landa[(p, 1)], 3):
                        print(i)
                        path_dron.append(P[p].V[i])
                path_dron.append([a[p, 0, 0].X, a[p, 0, 1].X])

            path_dron.append([xr[t, 0].X, xr[t, 1].X])
            ax.add_patch(Polygon(path_dron, fill=False, animated=True,
                         linestyle='-', alpha=0.5, color='black'))
            for p in path_dron:
                plt.plot(p[0], p[1], 'ko', markersize=1,
                         alpha=1, color='black')

for p in path_camion:
    plt.plot(p[0], p[1], 'ko', markersize=3, alpha=1, color='red')
#    ax.add_artist(Circle((p[0], p[1]), radius = 1.0))


# polygon_dron = Polygon(path_dron, fill=False, animated=True, linestyle='-', alpha=1, color = 'black')
# ax.add_patch(polygon_dron)

polygon_camion = Polygon(path_camion, fill=False,
                         animated=True, linestyle=':', alpha=1, color='red')
ax.add_patch(polygon_camion)


plt.title(str(nP) + "coordinate ")
# string = 'imagenes/' + str(m) + "-XPPN - MTZ - Mode " + \
#     str(modo) + " - Radii - " + str(datos.r) + ".png"
plt.savefig('imagenes\modelo7.png')
#plt.show()
#plt.close()
plt.show()
