############################################################################

# Created by: Prof. Valdecy Pereira, D.Sc.
# UFF - Universidade Federal Fluminense (Brazil)
# email:  valdecy.pereira@gmail.com
# Course: Metaheuristics
# Lesson: Variable Neighborhood Search

# Citation:
# PEREIRA, V. (2018). Project: Metaheuristic-Local_Search-Variable_Neighborhood_Search, File: Python-MH-Local Search-Variable Neighborhood Search.py, GitHub repository: <https://github.com/Valdecy/Metaheuristic-Local_Search-Variable_Neighborhood_Search>

############################################################################

# Required Libraries
import pandas as pd
import numpy as np
import random
import copy
from matplotlib import pyplot as plt

# Function: Tour Distance
def distance_calc(X, data, city_tour):
    "X = ([(x1, x2, city_tour, i, componente)], data) donde i = 0 si solo hay una componente, y 1 si hay dos"

    restaruno = lambda x: x-1
    distance = 0
    for c in range(len(city_tour[0])-1):
        m = c + 1
        loc1 = city_tour[0][c]
        loc2 = city_tour[0][m]
        punto1 = X[loc1]
        punto2 = X[loc2]
        comp = data[int(punto1[4])]
        if punto1[4] == punto2[4]:
            distance = distance + comp.alpha*comp.longitud
        else:
            distance = distance + np.linalg.norm(punto1[0:2] - punto2[0:2])
    return distance

# # Function: Initial Seed
# def seed_function(X, data):
#     seed = [[],float("inf")]
#     sequence = random.sample(list(range(1,X.shape[0]+1)), X.shape[0])
#     sequence.append(sequence[0])
#     seed[0] = sequence
#     seed[1] = distance_calc(X, data, seed)
#     print(seed[0])
#     return seed

def seed_function(X, data):
    seed = [[],float("inf")]
    sequence = list(range(X.shape[0]))
    sequence.append(sequence[0])
    seed[0] = sequence
    seed[1] = distance_calc(X, data, seed)
    return seed

# Function: Build Distance Matrix
def buid_distance_matrix(coordinates):
    X = np.zeros((coordinates.shape[0], coordinates.shape[0]))
    for i in range(0, X.shape[0]):
        for j in range(0, X.shape[1]):
            if (i != j):
                x = coordinates[i,:]
                y = coordinates[j,:]
                X[i,j] = np.linalg.norm(x - y)
    return X

def buscar(punto, puntos):
    m = len(puntos)
    for ind1 in range(m):
        for ind2 in range(2):
            if puntos[ind1][ind2] == punto:
                return (ind1, ind2)
    return -1

# Function: Stochastic 2_opt
def stochastic_2_opt(X, data, puntos, city_tour):
    best_route = copy.deepcopy(city_tour)
    i, j  = random.sample(range(len(city_tour[0])-1), 2)
    minc = min(i, j)
    maxc = max(i, j)
    i = minc
    j = maxc
#    print(str(i)+' ' +str(j))
    if buscar(city_tour[0][i], puntos) != -1:
        ind1, ind2 = buscar(city_tour[0][i], puntos)
        i1 = i
        i2 = [c for c in range(len(city_tour[0])-1) if city_tour[0][c] == puntos[ind1][1-ind2]][0]
        i = min(i1, i2)


    if buscar(city_tour[0][j], puntos) != -1:
        ind1, ind2 = buscar(city_tour[0][j], puntos)
        j1 = j
        j2 = [c for c in range(len(city_tour[0])-1) if city_tour[0][c] == puntos[ind1][1-ind2]][0]
        j = max(j1, j2)

    best_route[0][i:j+1] = list(reversed(best_route[0][i:j+1]))
#    print(str(i)+' ' +str(j))
#    print(best_route[0])
    best_route[0][-1]  = best_route[0][0]
    best_route[1] = distance_calc(X, data, best_route)
    return best_route

# Function: Local Search
def local_search(X, data, puntos, city_tour, max_attempts = 50, neighbourhood_size = 5):
    count = 0
    solution = copy.deepcopy(city_tour)
    while (count < max_attempts):
        for i in range(0, neighbourhood_size):
            candidate = stochastic_2_opt(X, data, puntos, city_tour = solution)
        if candidate[1] < solution[1]:
            solution  = copy.deepcopy(candidate)
            count = 0
        else:
            count = count + 1
    return solution

# Function: Variable Neighborhood Search
def variable_neighborhood_search(X, data, puntos, city_tour, max_attempts = 20, neighbourhood_size = 5, iterations = 50):
    count = 0
    solution = copy.deepcopy(city_tour)
    best_solution = copy.deepcopy(city_tour)
    while (count < iterations):
        for i in range(0, neighbourhood_size):
            for j in range(0, neighbourhood_size):
                solution = stochastic_2_opt(X, data, puntos, city_tour = best_solution)
            solution = local_search(X, data, puntos, city_tour = solution, max_attempts = max_attempts, neighbourhood_size = neighbourhood_size )
            if (solution[1] < best_solution[1]):
                best_solution = copy.deepcopy(solution)
                break
        count = count + 1
        if count % 5 == 0:
            print("Iteration = ", count, "-> Distance ", best_solution[1])
    return best_solution

######################## Part 1 - Usage ####################################

# Load File - A Distance Matrix (17 cities,  optimal = 1922.33)
#X = pd.read_csv('Python-MH-Local Search-Variable Neighborhood Search-Dataset-01.txt', sep = '\t')
#X = X.values
#
## Start a Random Seed
#seed = seed_function(X)
#
## Call the Function
#lsvns = variable_neighborhood_search(X, city_tour = seed, max_attempts = 25, neighbourhood_size = 5, iterations = 1000)
#
## Plot Solution. Red Point = Initial city; Orange Point = Second City # The generated coordinates (2D projection) are aproximated, depending on the data, the optimum tour may present crosses
#plot_tour_distance_matrix(X, lsvns)
#
######################### Part 2 - Usage ####################################
#
## Load File - Coordinates (Berlin 52,  optimal = 7544.37)
#Y = pd.read_csv('Python-MH-Local Search-Variable Neighborhood Search-Dataset-02.txt', sep = '\t')
#Y = Y.values
#
## Build the Distance Matrix
#X = buid_distance_matrix(Y)
#
## Start a Random Seed
#seed = seed_function(X)
#
## Call the Function
#lsvns = variable_neighborhood_search(X, city_tour = seed, max_attempts = 25, neighbourhood_size = 5, iterations = 1000)
#
## Plot Solution. Red Point = Initial city; Orange Point = Second City
#plot_tour_coordinates(Y, lsvns)
