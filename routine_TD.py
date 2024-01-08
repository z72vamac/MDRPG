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
from TDMTZ import TDMTZ
from TDSEC import TDSEC
from TDST import TDST
from tsp_heuristic import heuristic
import csv
import pandas as pd
import pickle as pickle

instancias_grid = pickle.load(open("instancias_grid.pickle", "rb"))
instancias_deulonay = pickle.load(open("instancias_deulonay.pickle", "rb"))
instancias_ciclos = pickle.load(open("instancias_ciclos.pickle", "rb"))

dataframe = pd.DataFrame(columns=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form'])

# dataframe_h = pd.DataFrame(columns=['Obj', 'Time', 'Type'])

for i in range(5):
    datos = instancias_grid[i]
    ciclo = instancias_ciclos[i]

    print()
    print('--------------------------------------------')
    print('PMDRPG-Stages: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Grid"))
    print('--------------------------------------------')
    print()

    sol_Stages = TDST(datos, ciclo)

    print()
    print('--------------------------------------------')
    print('PMDRPG-MTZ: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Grid"))
    print('--------------------------------------------')
    print()

    sol_MTZ = TDMTZ(datos, ciclo)

    print()
    print('--------------------------------------------')
    print('PMDRPG-SEC: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Grid"))
    print('--------------------------------------------')
    print()

    sol_SEC = TDSEC(datos, ciclo)

    dataframe = dataframe.append(pd.Series([sol_Stages[0], sol_Stages[1], sol_Stages[2],sol_Stages[3], sol_Stages[4], sol_Stages[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)

    dataframe = dataframe.append(pd.Series([sol_MTZ[0], sol_MTZ[1], sol_MTZ[2],sol_MTZ[3], sol_MTZ[4], sol_MTZ[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)

    dataframe = dataframe.append(pd.Series([sol_SEC[0], sol_SEC[1], sol_SEC[2],sol_SEC[3], sol_SEC[4], sol_SEC[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)
    dataframe.to_csv('PMDRPG_results' + '.csv', header = True, mode = 'w')

    # print()
    # print('--------------------------------------------')
    # print('AMDRPG-Heuristic: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Grid"))
    # print('--------------------------------------------')
    # print()
    #
    # sol_h= heuristic(datos)
    #
    # dataframe_h = dataframe_h.append(pd.Series([sol_h[0], sol_h[1], sol_h[2]], index=['Obj', 'Time', 'Type']), ignore_index=True)
    #
    # dataframe_h.to_csv('Heuristic_results' + '.csv', header = True, mode = 'w')

    datos2 = instancias_deulonay[i]

    print()
    print('--------------------------------------------')
    print('PMDRPG-Stages: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Deulonay"))
    print('--------------------------------------------')
    print()

    sol_Stages = TDST(datos2, ciclo)

    print()
    print('--------------------------------------------')
    print('PMDRPG-MTZ: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Deulonay"))
    print('--------------------------------------------')
    print()

    sol_MTZ = TDMTZ(datos2, ciclo)

    print()
    print('--------------------------------------------')
    print('PMDRPG-SEC: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Deulonay"))
    print('--------------------------------------------')
    print()

    sol_SEC = TDSEC(datos2, ciclo)

    dataframe = dataframe.append(pd.Series([sol_Stages[0], sol_Stages[1], sol_Stages[2],sol_Stages[3], sol_Stages[4], sol_Stages[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)

    dataframe = dataframe.append(pd.Series([sol_MTZ[0], sol_MTZ[1], sol_MTZ[2],sol_MTZ[3], sol_MTZ[4], sol_MTZ[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)

    dataframe = dataframe.append(pd.Series([sol_SEC[0], sol_SEC[1], sol_SEC[2],sol_SEC[3], sol_SEC[4], sol_SEC[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)
    dataframe.to_csv('PMDRPG_results' + '.csv', header = True, mode = 'w')


    # print()
    # print('--------------------------------------------')
    # print('AMDRPG-Heuristic: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Delaunay"))
    # print('--------------------------------------------')
    # print()
    #
    # sol_h= heuristic(datos2)
    #
    # dataframe_h = dataframe_h.append(pd.Series([sol_h[0], sol_h[1], sol_h[2]], index=['Obj', 'Time', 'Type']), ignore_index=True)
    #
    # dataframe_h.to_csv('Heuristic_results' + '.csv', header = True, mode = 'w')
