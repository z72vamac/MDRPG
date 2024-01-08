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


dataframe = pd.DataFrame(columns=['Instance', 'List', 'Gap', 'Time', 'Nodes', 'ObjVal', 'ObjVal_h', 'Time_h', 'Structure', 'Pct', 'Model'])

it = 0
its = [73]
for i, j, alpha, grid in instancias.keys():
    # if i == 2 and j == 1 and alpha == True and grid == False:
    print_instance(i, j, alpha, grid)

    if it in its:
        datos = instancias[i, j, alpha, grid]
        sol_MTZ = PDMTZ(datos)

        lista = [i, j] + sol_MTZ
        serie = pd.Series(lista, index = dataframe.columns)
        print(serie)

        dataframe = dataframe.append(pd.Series(serie, index=['Instance', 'List', 'Gap', 'Time', 'Nodes', 'ObjVal', 'ObjVal_h', 'Time_h', 'Structure', 'Pct', 'Model']), ignore_index=True)
        dataframe.to_csv('results_heuristic2.csv')

    it += 1
