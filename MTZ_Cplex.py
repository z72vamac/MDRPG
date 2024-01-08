#!/usr/bin/python

# Copyright 2019, Gurobi Optimization, LLC

# Solve the classic diet model, showing how to add constraints
# to an existing model.

import docplex.mp.model as cpx
from docplex.mp.model import Model
from docplex.util.environment import get_environment
import numpy as np
from entorno import Poligono, Elipse, Poligonal
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle
import estimacion_M as eM
from data import *
from itertools import combinations
import auxiliar_functions as af
import time
# Con 15 todavia no la encuentra
# np.random.seed(112)
# np.random.seed(11112)
# 168.60477938110762
np.random.seed(11)
datos = Data([], m = 8,
                 r = 1,
                 modo = 1,
                 tmax = 120,
                 init = False,
                 show = True,
                 seed = 2)
datos.generar_muestra()


def MTZ(datos):

    data = datos.data
    m = len(data)

    # modo = datos.modo

    # Model
    M = Model(name='MTZ_MODEL') #, **kwargs)

    # Generando variables continuas de entrada y salida
    x_index = []
    d_index = []
    difc_index = []


    for c in range(m):
        # comp = data[c]
        # if type(comp) is Poligonal:
        x_index.append((c, 0))
        x_index.append((c, 1))
        d_index.append(c)

    x = M.continuous_var_dict(x_index, name='x')

    z_index = []

    for c1 in d_index:
        for c2 in d_index:
            if c1 != c2:
                z_index.append((c1, c2))

    p = M.continuous_var_dict(z_index, lb=0, name='p')

    z = M.binary_var_dict(z_index, name='z')
    d = M.continuous_var_dict(z_index, lb=0, name='dcc')

    # if datos.init:
    #     d = M.addVars(z_index, vtype = GRB.CONTINUOUS, lb = ds, name = 'dcc')
    dif_index = []

    for c1 in d_index:
        for c2 in d_index:
         for i in range(2):
          if c1 != c2:
            dif_index.append((c1, c2, i))

    dif = M.continuous_var_dict(dif_index, lb=0, name='dif')

    dc = M.continuous_var_dict(d_index, lb=0, name='dc')

    for c in range(m):
     for i in range(2):
       difc_index.append((c,i))

    difc = M.continuous_var_dict(difc_index, lb=0, name='difc')

    # Generando los mus de la envolvente convexa, los landas de la poligonal y las
    # variables binarias que indican qu segmento se elige

    mu_index = []
    landa_index = []
    sublanda_index = []
    u_index = []
    s_index = []

    for c in d_index:
        comp = data[c]
        if type(comp) is Poligono:
            for mu in range(comp.num_puntos):
                mu_index.append((c, mu))
        if type(comp) is Poligonal:
            u_index.append(c)
            # landa de la variable de entrada en la poligonal c
            landa_index.append((c, 0))
            # landa de la variable de salida en la poligonal c
            landa_index.append((c, 1))
            for segm in range(comp.num_segmentos):
                s_index.append((c, 0, segm))
                s_index.append((c, 1, segm))
            for punto in range(comp.num_puntos):
                sublanda_index.append((c, 0, punto))
                sublanda_index.append((c, 1, punto))

    mu = M.continuous_var_dict(mu_index, lb=0.0, ub=1.0, name='mu')
    landa = M.continuous_var_dict(landa_index, name='landa')
    sublanda = M.continuous_var_dict(sublanda_index, lb=0.0, ub=1.0, name='sublanda')
    s = M.binary_var_dict(s_index, name='s')
    u = M.binary_var_dict(u_index, name='u')
    minc = M.continuous_var_dict(u_index, lb=0.0, name='min')
    maxc = M.continuous_var_dict(u_index, lb=0.0, name='max')

    sc = M.continuous_var_dict(d_index, lb=0, ub=m - 1, name='sc')



    #M.update()

    #
    #     for c in sc.keys():
    #         sc[c].start = scs[c]
#            d[i, j].start = ds[i, j]

    # for c in s.keys():
    #     s[c].start = ss[c]
    #
    # for c in u.keys():
    #     u[c].start = us[c]
#        for c in dc.keys():
#            dc[c].start = dcs[c]

    # Constraints
    for c1, c2 in z.keys():
        comp1 = data[c1]
        comp2 = data[c2]
        BigM = eM.estima_BigM_local(comp1, comp2)
        SmallM = eM.estima_SmallM_local(comp1, comp2)
        # M.addConstr(d[c1, c2] <= BigM)
        # M.addConstr(d[c1, c2]>= SmallM)
        # SmallM, x_0, x_1 = af.min_dist(comp1, comp2)
        M.add_constraint(p[c1, c2] >= SmallM * z[c1, c2] , ctname='p2')
        M.add_constraint(p[c1, c2] >= d[c1, c2] - BigM * (1 - z[c1, c2]), ctname='p3')
        #
        # M.addConstr(p[c1, c2] >= SmallM*z[c1, c2])
        # M.addConstr(p[c1, c2] <= d[c1, c2] + z[c1, c2]*SmallM - SmallM)
        # M.addConstr(p[c1, c2] <= z[c1, c2]*BigM)

        for dim in range(2):
            M.add_constraint(dif[c1, c2, dim] >=  x[c2, dim] - x[c1, dim])
            M.add_constraint(dif[c1, c2, dim] >= -x[c2, dim] + x[c1, dim])


        M.add_constraint(dif[c1, c2, 0] * dif[c1, c2, 0] +
                    dif[c1, c2, 1] * dif[c1, c2, 1] <= d[c1, c2] * d[c1, c2])


    # Restricciones para formar un tour
    M.add_constraint(M.sum(z[v, 0] for v in range(1, m-1)) == 0)
    M.add_constraint(M.sum(z[m-1, w] for w in range(1, m-1)) == 0)

    for v in range(m):
        M.add_constraint(M.sum(z[v , w] for w in range(m) if w != v) == 1)
        M.add_constraint(M.sum(z[w , v] for w in range(m) if w != v) == 1)

    for v in range(1, m-1):
        for w in range(1, m-1):
            if v != w:
                M.add_constraint(m - 1 >= (sc[v] - sc[w]) + m * z[v, w])

    # for v in range(1, nG+1):
    #     M.add_constraint(s[v] - s[0] + (nG+1 - 2)*z[0, v] <= len(T_index) - 1)
    #
    # for v in range(1, nG+1):
    #     M.add_constraint(s[0] - s[v] + (nG+1 - 1)*z[v, 0] <= 0)

    # for v in range(1, nG+1):
    #     M.add_constraint(-z[0,v] - s[v] + (nG+1-3)*z[v,0] <= -2, name="LiftedLB(%s)"%v)
    #     M.add_constraint(-z[v,0] + s[v] + (nG+1-3)*z[0,v] <= nG+1-2, name="LiftedUB(%s)"%v)

    #     M.add_constraint(s[v] >= 1)
    #     for v in range(1, m-1):
    #     M.add_constraint(s[v] <= len(T_index) - 1)
    #
    # M.add_constraint(s[0] == 0)
    # M.add_constraint(s[m-1] == m-1)


    # Restriccion 10
    M.add_constraint(sc[0] == 0, ctname = 'rest10')
    M.add_constraint(sc[m-1] == m-1)

    for c in range(1, m-1):
        M.add_constraint(sc[c] >= 1, ctname = 'rest11')

    for c in range(m):
        for dim in range(2):
            M.add_constraint(difc[c, dim] >=  x[c, dim] - data[c].centro[dim])
            M.add_constraint(difc[c, dim] >= -x[c, dim] + data[c].centro[dim])
            M.add_constraint(difc[c, 0]*difc[c, 0] + difc[c, 1]*difc[c, 1] <= data[c].radio*data[c].radio)

    # M.update()

    fc = np.ones(m) * 1

    # for c in range(m):
    #     if fc[c] >= 1.0 and type(data[c]) is not Poligonal:
    #         for j in range(2):
    #             M.add_constraint(x[c, 0, j] == x[c, 1, j])

    objective = M.sum(p[c1, c2] for c1 in range(m) for c2 in range(m) if c1 != c2)

    # objective = gp.quicksum(p[index] for index in p.keys())

    M.minimize(objective)

    # M.update()

    #M.add_constraint(M.getObjective() <= obj_heur)
    #M.add_constraint(M.getObjective() >= suma)


    # if datos.tmax is not None:
    #     M.Params.timeLimit = datos.tmax
    #
    # if not(datos.show):
    #     M.setParam('OutputFlag', 0)
    #
    # M.Params.Threads = 8
    # M.Params.tuneTimeLimit = 600
    # M.Params.timeLimit = 20
    # M.tune()
    # M.Params.Method = 4
    # M.Params.Heuristic = 0
    # M.Params.IntFeasTol = 1e-1
    # M.Params.NonConvex = 2
    # M.Params.FeasibilityTol = 1e-6
    #
    # M.update()
    #
    # M.write('MTZ.lp')

    # Solve
    M.solve()

    print(M.z.solution_value())
    print(M.solution.get_objective_value())
    print()
    # vals = M.getAttr('x', z)
    # selected = gp.tuplelist((i, j)
    #                         for i, j in vals.keys() if vals[i, j] > 0.5)
    #
    # path = af.subtour(selected)
    #
    # print('Camino = ' + str(path))
    #
    # path_P = []
    #
    # for p in path:
    #     comp = data[p]
    #     path_P.append([x[(p, 0)].X, x[(p, 1)].X])

    #return path, path_P, M.ObjVal

MTZ(datos)
