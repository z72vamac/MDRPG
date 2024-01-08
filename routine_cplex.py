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
from PDMTZ_Cplex import PDMTZ_cplex
from PDSEC_Cplex import PDSEC_cplex
from PDST_Cplex import PDST_cplex
import csv
import pandas as pd
import pickle as pickle

instancias_grid = pickle.load(open("instancias_grid.pickle", "rb"))
instancias_deulonay = pickle.load(open("instancias_deulonay.pickle", "rb"))

dataframe = pd.DataFrame(columns=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form'])

dataframe_h = pd.DataFrame(columns=['Obj', 'Time', 'Type'])

for i in range(5):
    datos = instancias_grid[i]

    # print()
    # print('--------------------------------------------')
    # print('AMDRPG-Stages: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Grid"))
    # print('--------------------------------------------')
    # print()

    # sol_Stages = PDST_cplex(datos)

    # print()
    # print('--------------------------------------------')
    # print('AMDRPG-MTZ: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Grid"))
    # print('--------------------------------------------')
    # print()

    # sol_MTZ = PDMTZ_cplex(datos)

    print()
    print('--------------------------------------------')
    print('AMDRPG-SEC: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Grid"))
    print('--------------------------------------------')
    print()

    sol_SEC = PDSEC_cplex(datos)

    # dataframe = dataframe.append(pd.Series([sol_Stages[0], sol_Stages[1], sol_Stages[2],sol_Stages[3], sol_Stages[4], sol_Stages[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)

    # dataframe = dataframe.append(pd.Series([sol_MTZ[0], sol_MTZ[1], sol_MTZ[2],sol_MTZ[3], sol_MTZ[4], sol_MTZ[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)

    dataframe = dataframe.append(pd.Series([sol_SEC[0], sol_SEC[1], sol_SEC[2],sol_SEC[3], sol_SEC[4], sol_SEC[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)
    dataframe.to_csv('AMDRPG_results_cplex' + '.csv', header = True, mode = 'w')

    # print()
    # print('--------------------------------------------')
    # print('AMDRPG-Heuristic: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Grid"))
    # print('--------------------------------------------')
    # print()
    #
    # sol_h= heuristic_cplex(datos)
    #
    # dataframe_h = dataframe_h.append(pd.Series([sol_h[0], sol_h[1], sol_h[2]], index=['Obj', 'Time', 'Type']), ignore_index=True)
    #
    # dataframe_h.to_csv('Heuristic_results_cplex' + '.csv', header = True, mode = 'w')

    datos2 = instancias_deulonay[i]

    # print()
    # print('--------------------------------------------')
    # print('AMDRPG-Stages: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Delaunay"))
    # print('--------------------------------------------')
    # print()

    # sol_Stages = PDST_cplex(datos2)

    # print()
    # print('--------------------------------------------')
    # print('AMDRPG-MTZ: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Delaunay"))
    # print('--------------------------------------------')
    # print()

    # sol_MTZ = PDMTZ_cplex(datos2)

    print()
    print('--------------------------------------------')
    print('AMDRPG-SEC: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Delaunay"))
    print('--------------------------------------------')
    print()

    sol_SEC = PDSEC_cplex(datos2)

    # dataframe = dataframe.append(pd.Series([sol_Stages[0], sol_Stages[1], sol_Stages[2],sol_Stages[3], sol_Stages[4], sol_Stages[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)
    # dataframe.to_csv('AMDRPG_results_cplex' + '.csv', header = True, mode = 'w')

    # dataframe = dataframe.append(pd.Series([sol_MTZ[0], sol_MTZ[1], sol_MTZ[2],sol_MTZ[3], sol_MTZ[4], sol_MTZ[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)

    dataframe = dataframe.append(pd.Series([sol_SEC[0], sol_SEC[1], sol_SEC[2],sol_SEC[3], sol_SEC[4], sol_SEC[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)
    dataframe.to_csv('AMDRPG_results_cplex' + '.csv', header = True, mode = 'w')


    # print()
    # print('--------------------------------------------')
    # print('AMDRPG-Heuristic: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Delaunay"))
    # print('--------------------------------------------')
    # print()
    #
    # sol_h= heuristic_cplex(datos2)
    #
    # dataframe_h = dataframe_h.append(pd.Series([sol_h[0], sol_h[1], sol_h[2]], index=['Obj', 'Time', 'Type']), ignore_index=True)

    # dataframe_h.to_csv('Heuristic_results_cplex' + '.csv', header = True, mode = 'w')


# grid= False
#
# for i in range(5):
#     datos = Data([], m = nG, r = 1, grid = grid, tmax =tmax, alpha = alpha, init = True, show = True, elim=False, prepro=False, refor=False, pol=False, porc=1, seed=0+i)
#
#     datos.generar_grafos(lista)
#     sol = gurobi(datos)
#     solutions.append(sol)
