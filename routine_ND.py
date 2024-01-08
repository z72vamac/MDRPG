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
from NDMTZ import NDMTZ
from NDST import NDST
import csv
import pandas as pd
import pickle as pickle

instancias = pickle.load(open("instancias_ND.pickle", "rb"))

dataframe = pd.DataFrame(columns=['Instance', 'List', 'GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Pct', 'Model', 'Mode'])

def print_instance(i, j, alpha, grid):
    string1 = '#########################################################################'
    string = '# ND -- Instance #' + str(i) + '-- List: ' + str(j)
    string2 = '#########################################################################'

    if alpha:
        string += '-- Edge % '
    else:
        string += '-- Graph %'

    if grid:
        string1 += '\n'
        string += '-- Graph structure: Grid #\n'
        string2 += '\n\n'
    else:
        string1 += '####\n'
        string += '-- Graph structure: Delauney #\n'
        string2 += '####\n\n'


    total_string = string1 + string + string2
    print(total_string)

it = 0
# dataframe_h = pd.DataFrame(columns=['Obj', 'Time', 'Type'])
# ha hecho hasta el 3-0-1
for i, j, alpha, grid in instancias:

    print_instance(i, j, alpha, grid)

    datos = instancias[i, j, alpha, grid]

    for mode in range(1, 4):

        it += 1

        grafo_data = Data([], m = 1, grid = True, tmax=600, alpha = True,
                        init = True,
                        show = True,
                        seed = 2)
        grafo_data.generar_grafo_personalizado(mode = mode)

        print('Mode = ' + str(mode))
        print('Stages formulation')

        sol_Stages = NDST(datos, grafo_data)

        lista = [i, j] + sol_Stages
        serie = pd.Series(lista, index = dataframe.columns)
        print(serie)

        dataframe = dataframe.append(pd.Series([sol_Stages[0], sol_Stages[1], sol_Stages[2], sol_Stages[3], sol_Stages[4], sol_Stages[5], sol_Stages[6], sol_Stages[7]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Pct', 'Model', 'Mode']), ignore_index=True)

        dataframe.to_csv('./results/ND_results2' + '.csv', header = True, mode = 'w')

        print('Mode = ' + str(mode))
        print('MTZ formulation')

        sol_MTZ = NDMTZ(datos, grafo_data)
        dataframe = dataframe.append(pd.Series([sol_MTZ[0], sol_MTZ[1], sol_MTZ[2], sol_MTZ[3], sol_MTZ[4], sol_MTZ[5], sol_MTZ[6], sol_MTZ[7]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Pct', 'Model', 'Mode']), ignore_index=True)

        lista = [i, j] + sol_Stages
        serie = pd.Series(lista, index = dataframe.columns)
        print(serie)

        dataframe.to_csv('./results/ND_results2' + '.csv', header = True, mode = 'w')


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
