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

it = -1
dataframe = pd.DataFrame(columns=['Instance', 'List', 'GAP', 'Time', 'Nodes', 'ObjVal', 'ObjVal_h', 'Time_h', 'Structure', 'Pct', 'Model', 'Mode'])

for i, j, alpha, grid in instancias:

    print_instance(i, j, alpha, grid)

    datos = instancias[i, j, alpha, grid]

    for mode in range(1, 4):

        it += 1

        if it == 42:

            grafo_data = Data([], m = 1, grid = True, tmax=600, alpha = True,
                            init = True,
                            show = True,
                            seed = 2)
            grafo_data.generar_grafo_personalizado(mode = mode)

            print('Mode = ' + str(mode))
            print('MTZ formulation')

            sol_MTZ = NDMTZ(datos, grafo_data)

            print(sol_MTZ)
            lista = [i, j] + sol_MTZ
            serie = pd.Series(lista, index = dataframe.columns)
            print(serie)

            dataframe = dataframe.append(pd.Series(serie, index=['Instance', 'List', 'GAP', 'Time', 'Nodes', 'ObjVal', 'ObjVal_h', 'Time_h', 'Structure', 'Pct', 'Model', 'Mode']), ignore_index=True)
            dataframe.to_csv('results_init_ND6.csv')


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
