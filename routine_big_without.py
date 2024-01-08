#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gurobipy as gp
import pdb
#from gurobipy import GRB
#import numpy as np
#from itertools import product
#import random
#import matplotlib.pyplot as plt
#from matplotlib.patches import Circle, Polygon
#from matplotlib.collections import PatchCollection
#import matplotlib.lines as mlines
#from data import *
#from entorno import *
#import copy
#import estimacion_M as eM
#import networkx as nx
#import auxiliar_functions as af
from PDMTZ import PDMTZ
#from PDSEC import PDSEC
#from PDST import PDST
#from PDMTZ_heuristic import PDMTZ_heuristic
#import csv
import pandas as pd
import pickle as pickle

instancias = pickle.load(open("instancias.pickle", "rb"))
# instancias_deulonay = pickle.load(open("instancias_deulonay.pickle", "rb"))

# dataframe_h = pd.DataFrame(columns=['Obj', 'Time', 'Type'])
# print(instancias)
def print_instance(i, j, alpha, grid):
    string1 = '#########################################################################'
    string = '# AMDRPG-MTZ -- Instance #' + str(i) + '-- List: ' + str(j)
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


dataframe = pd.DataFrame(columns=['Instance', 'List', 'Gap', 'Time', 'Nodes', 'ObjVal', 'Structure', 'Pct', 'Model'])

for i, j, alpha, grid in instancias.keys():
    # if i == 2 and j == 1 and alpha == True and grid == False:
        print_instance(i, j, alpha, grid)


        datos = instancias[i, j, alpha, grid]
        datos.init = False
        sol_MTZ = PDMTZ(datos)

        lista = [i, j] + sol_MTZ
        serie = pd.Series(lista, index = dataframe.columns)
        print(serie)

        dataframe = dataframe.append(pd.Series(serie, index=['Instance', 'List', 'Gap', 'Time', 'Nodes', 'ObjVal', 'Structure', 'Pct', 'Model']), ignore_index=True)
        dataframe.to_csv('./results_noinit.csv')


        # dataframe = dataframe.append(pd.Series(serie, index=['Instance', 'List', 'Gap', 'Time', 'Nodes', 'ObjVal', 'Structure', 'Pct', 'Model']), ignore_index=True)
        # dataframe.to_csv('results.csv', mode = 'w')


    # print()
    # print('--------------------------------------------')
    # print('AMDRPG-Stages: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Grid"))
    # print('--------------------------------------------')
    # print()

    # sol_Stages = PDST(datos)

    # print()
    # print('--------------------------------------------')
    # print('AMDRPG-MTZ: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Grid"))
    # print('--------------------------------------------')
    # print()

    # sol_MTZ = PDMTZ(datos)

    # print()
    # print('--------------------------------------------')
    # print('AMDRPG-SEC: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Grid"))
    # print('--------------------------------------------')
    # print()

    # sol_SEC = PDSEC(datos)

    # dataframe = dataframe.append(pd.Series([sol_Stages[0], sol_Stages[1], sol_Stages[2],sol_Stages[3], sol_Stages[4], sol_Stages[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)

    # dataframe = dataframe.append(pd.Series([sol_MTZ[0], sol_MTZ[1], sol_MTZ[2],sol_MTZ[3], sol_MTZ[4], sol_MTZ[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)

    # dataframe = dataframe.append(pd.Series([sol_SEC[0], sol_SEC[1], sol_SEC[2],sol_SEC[3], sol_SEC[4], sol_SEC[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)
    # dataframe.to_csv('AMDRPG_results' + '.csv', header = True, mode = 'w')

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
    #
    # datos2 = instancias_deulonay[i]

    # print()
    # print('--------------------------------------------')
    # print('AMDRPG-Stages: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Delaunay"))
    # print('--------------------------------------------')
    # print()

    # sol_Stages = PDST(datos2)

    # print()
    # print('--------------------------------------------')
    # print('AMDRPG-MTZ: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Delaunay"))
    # print('--------------------------------------------')
    # print()

    # sol_MTZ = PDMTZ(datos2)

    # print()
    # print('--------------------------------------------')
    # print('AMDRPG-SEC: Iteration: {i} - Graphs: {j}'.format(i = i, j = "Delaunay"))
    # print('--------------------------------------------')
    # print()

    # sol_SEC = PDSEC(datos2)

    # dataframe = dataframe.append(pd.Series([sol_Stages[0], sol_Stages[1], sol_Stages[2],sol_Stages[3], sol_Stages[4], sol_Stages[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)
    # dataframe.to_csv('AMDRPG_results' + '.csv', header = True, mode = 'w')

    # dataframe = dataframe.append(pd.Series([sol_MTZ[0], sol_MTZ[1], sol_MTZ[2],sol_MTZ[3], sol_MTZ[4], sol_MTZ[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)

    # dataframe = dataframe.append(pd.Series([sol_SEC[0], sol_SEC[1], sol_SEC[2],sol_SEC[3], sol_SEC[4], sol_SEC[5]], index=['GAP', 'Time', 'Nodes', 'Obj', 'Type', 'Form']), ignore_index=True)
    # dataframe.to_csv('AMDRPG_results' + '.csv', header = True, mode = 'w')


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
