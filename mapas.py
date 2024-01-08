import networkx as nx
import osmnx as ox
import requests
import matplotlib.cm as cm
import matplotlib.colors as colors
# %matplotlib inline

ox.config(use_cache=True, log_console=True)
ox.__version__

G = ox.graph_from_place('Sevilla, Spain', network_type='drive')
ox.plot_graph(G)
