
# Incluimos primero los paquetes
import gurobipy as gp
from gurobipy import GRB
import numpy as np
from itertools import product
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from data import *
from entorno import *
import estimacion_M as eM

# Definicion de los datos
""" M: numero de poligonales
    N: dimension del problema
    data: conjunto de poligonales a agrupar
    J: numero de entornos donde pueden ir los centroides
    C(j): centro del entorno j
    R(j): radio del entorno j
    i: indice de las poligonales
    j: indice de los entornos

"""


# Datos iniciales
# np.random.seed(1)
# random.seed(1)
M, N = 6, 2
J = M
L = 2

vP = 1
vD = 10

datos = Data([], m = M,
                 r = 3,
                 modo = 3,
                 tmax = 120,
                 init = True,
                 show = True,
                 seed = 2)

datos.generar_muestra()
data = datos.mostrar_datos()

entornos = Data([], m = J,
                    r = 2,
                    modo = 1,
                    tmax = 120,
                    init = True,
                    show = True,
                    seed = 2)

entornos.generar_muestra()
entorno = entornos.mostrar_datos()


# Creacion del modelo
MODEL = gp.Model('klmedian-continuous')

# Subindices de la variable continua xkl
x_index = []

for j in range(J):
    for n in range(N):
        x_index.append((j, n))

x = MODEL.addVars(x_index, 2, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'x')

# Subindices de la variable binaria yij
y_index = []

for i, j in product(range(M), range(J)):
    y_index.append((i, j))

y = MODEL.addVars(y_index, vtype = GRB.BINARY, name = 'y')


# Subindices de las variables continuas p, d
d_index = []

for i, j in product(range(M), range(J)):
    d_index.append((i, j))

d = MODEL.addVars(d_index, 2, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'd')
p = MODEL.addVars(d_index, 2, vtype = GRB.CONTINUOUS, lb = 0.0, name = 'p')

# Variables auxiliares para representar restricciones SOC de distancia entre ai y xj. Primeras dos variables correspondientes a
# la i j, tercera componente es el punto de referencia si ai1 o ai2, cuarta componente es la dimension.

dif = MODEL.addVars(d_index, 2, 2, vtype = GRB.CONTINUOUS, name = 'dif')

# Variables auxiliares para representar restricciones SOC de distancia entre xj y cj
dif2 = MODEL.addVars(x_index, 2, vtype = GRB.CONTINUOUS, name = 'dif2')

# Variables que representan los puntos de la poligonal por donde se entra y se sale
a_index = []

for i in range(M):
    a_index.append(i)

a = MODEL.addVars(a_index, 2, 2, vtype = GRB.CONTINUOUS, name = 'a')

# Variables que parametrizan la poligonal
mu_index = []
landa_index = []
sublanda_index = []
u_index = []
s_index = []

for i in range(M):
    comp = data[i]
    u_index.append(i)
    # landa de la variable de entrada en la poligonal c
    landa_index.append((i, 0))
    # landa de la variable de salida en la poligonal c
    landa_index.append((i, 1))
    for segm in range(comp.num_segmentos):
        s_index.append((i, 0, segm))
        s_index.append((i, 1, segm))
    for punto in range(comp.num_puntos):
        sublanda_index.append((i, 0, punto))
        sublanda_index.append((i, 1, punto))

landa = MODEL.addVars(landa_index, vtype=GRB.CONTINUOUS, name='landa')
sublanda = MODEL.addVars(sublanda_index, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='sublanda')
s = MODEL.addVars(s_index, vtype=GRB.BINARY, name='s')
u = MODEL.addVars(u_index, vtype=GRB.BINARY, name='u')
minc = MODEL.addVars(u_index, vtype=GRB.CONTINUOUS, lb=0.0, name='min')
maxc = MODEL.addVars(u_index, vtype=GRB.CONTINUOUS, lb=0.0, name='max')

MODEL.update()

# initval = False
# if initval:
#     ys = algorithm(A, K, 5000, 60)
#
#     for i, j in ys:
#         y[i, j].start = 1

# Restricciones

for i, j, l in d.keys():

    # Restricciones de la diferencia
    MODEL.addConstr( a[i, l, 0] - x[j, l, 0] <= dif[i, j, l, 0], name = 'dif-const-l0')
    MODEL.addConstr(-a[i, l, 0] + x[j, l, 0] <= dif[i, j, l, 0], name = 'dif-const+l0')
    MODEL.addConstr( a[i, l, 1] - x[j, l, 1] <= dif[i, j, l, 1], name = 'dif-const-l1')
    MODEL.addConstr(-a[i, l, 1] + x[j, l, 1] <= dif[i, j, l, 1], name = 'dif-const+l1')
    MODEL.addConstr(dif[i, j, l, 0] * dif[i, j, l, 0] + dif[i, j, l, 1] * dif[i, j, l, 1] <= d[i, j, l] * d[i, j, l], name = 'd-const0')

    # Restricciones de linearizar
    ent = entorno[j]
    pol = data[i]

    smallM = eM.estima_SmallM_local(ent, pol)
    bigM = eM.estima_BigM_local(ent, pol)

    MODEL.addConstr(p[i, j, l] >= smallM*y[i, j], name = 'lin-const-0')
    MODEL.addConstr(p[i, j, l] >= d[i, j, l] - bigM*(1 - y[i, j]), name = 'lin-const-2')


for j, l in product(range(J), range(2)):

    # Inicialización de los datos
    comp = entorno[j]
    C = comp.centro
    R = comp.radio

    # Restricciones de la diferencia
    MODEL.addConstr( C[0] - x[j, l, 0] <= dif2[j, l, 0], name = 'dif2-const-0')
    MODEL.addConstr(-C[0] + x[j, l, 0] <= dif2[j, l, 0], name = 'dif2-const+0')
    MODEL.addConstr( C[1] - x[j, l, 1] <= dif2[j, l, 1], name = 'dif2-const-1')
    MODEL.addConstr(-C[1] + x[j, l, 1] <= dif2[j, l, 1], name = 'dif2-const+1')
    MODEL.addConstr(dif2[j, l, 0] * dif2[j, l, 0] + dif2[j, l, 1] * dif2[j, l, 1] <= R*R, name = 'd2-const')

# Restriccion de que hay que recorrer un porcentaje de la poligonal
for i in range(M):
    comp = data[i]
    if type(comp) is Poligonal:
        MODEL.addConstr(landa[i, 0] - landa[i, 1] == maxc[i] - minc[i], name='u0')
        # si u = 0, entonces landa0 >= landa1
        MODEL.addConstr(maxc[i] + minc[i] == comp.alpha * comp.num_segmentos, name='u1')
        MODEL.addConstr(maxc[i] <= comp.num_segmentos * (1 - u[i]), name='u2')
        MODEL.addConstr(minc[i] <= comp.num_segmentos * u[i], name='u3')
        for j in range(2):
            for punto in range(1, comp.num_puntos):
                MODEL.addConstr(landa[i, j] - punto >= sublanda[i, j, punto] - comp.num_puntos * (1 - s[i, j, punto - 1]))
                MODEL.addConstr(landa[i, j] - punto <= sublanda[i, j, punto] + comp.num_puntos * (1 - s[i, j, punto - 1]))
            MODEL.addConstr(sublanda[i, j, 0] <= s[i, j, 0])
            MODEL.addConstr(sublanda[i, j, comp.num_puntos - 1] <= s[i, j, comp.num_puntos - 2])
            for punto in range(1, comp.num_puntos - 1):
                MODEL.addConstr(sublanda[i, j, punto] <= s[i, j, punto - 1] + s[i, j, punto])
            MODEL.addConstr(s.sum(i, j, '*') == 1)
            MODEL.addConstr(sublanda.sum(i, j, '*') == 1)
            for dim in range(2):
                MODEL.addConstr(a[i, j, dim] == gp.quicksum(sublanda[i, j, punto] * comp.V[punto][dim] for punto in range(comp.num_puntos)), name='seg1')

# Restricciones de asignacion
# A cada cluster tienes que asignarle al menos un elemento
# MODEL.addConstrs(y.quicksum('*', k) == 1 for k in range(K))

# Cada elemento se asocia solo a un cluster
MODEL.addConstrs(y.sum(i, '*') == 1 for i in range(M))
#MODEL.addConstrs(y.sum('*', j) == 1 for j in range(J))

#MODEL.addConstr(gp.quicksum(z[j] for j in range(J)) >= 1)

MODEL.update()

# Funcion objetivo
objective = gp.quicksum(p[u]/vD for u in p.keys())

MODEL.setObjective(objective, GRB.MINIMIZE)

MODEL.update()

MODEL.setParam('TimeLimit', 600)
MODEL.write('aver.lp')

# Optimizamos
MODEL.optimize()

# MODEL.computeIIS()
# MODEL.write('aver.ilp')

# Representamos la solucion
vals = MODEL.getAttr('x', y)
y_selected = gp.tuplelist(u for u in vals.keys() if vals[u] > 0.5)

# vals = MODEL.getAttr('x', z)
# z_selected = gp.tuplelist(u for u in vals.keys() if vals[u] > 0.5)

# Configuracion de los ejes
fig, ax = plt.subplots()
ax.autoscale_view()
ax.axis([0, 100, 0, 100])

colores = []
for j in range(J):
    rgb = np.random.rand(3,)
    colores.append(rgb)

# Representacion de los centros
for j in range(J):
    xs0 = [x[j, 0, 0].X, x[j, 0, 1].X]
    xs1 = [x[j, 1, 0].X, x[j, 1, 1].X]
    ax.add_artist(Circle((xs0[0], xs0[1]), radius = 1.0, color = colores[j]))
    ax.add_artist(Circle((xs1[0], xs1[1]), radius = 1.0, color = colores[j]))

# Representacion de los puntos
for i, j in y_selected:
    ax.add_artist(Circle((a[i, 0, 0].X, a[i, 0, 1].X), radius = 0.5, color = colores[j]))
    ax.add_artist(Circle((a[i, 1, 0].X, a[i, 1, 1].X), radius = 0.5, color = colores[j]))


# Representacion de las poligonales
for i in range(M):
    dato = data[i]
    ax.add_artist(dato.artist)


# Representacion de los entornos
for j in range(J):
    dato = entorno[j]
    ax.add_artist(dato.artist)


plt.show()
