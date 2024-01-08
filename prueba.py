import networkx as nx
import matplotlib.pyplot as plt

g = nx.random_tree(8, seed = 1)

nx.draw(g)

g.edges()
