# Get summary statistics for .dot files (network graphs generated using graphviz)
# Code mnodified from https://programminghistorian.org/en/lessons/exploring-and-analyzing-network-data-with-python#metrics-available-in-networkx
import argparse
import os
import glob
import csv

import networkx as nx
import pandas as pd
import numpy as np

from operator import itemgetter

def generate_statistics(name, dot_file):
    # Read dot file
    G = nx.Graph(nx.drawing.nx_pydot.read_dot(dot_file))

    # Density:  how interconnected a graph is in terms of a ratio of actual over possible connections
    density = nx.density(G)

    # Transitivity: ratio of all triangles over all possible triangles
    triadic_closure = nx.transitivity(G)

    # Degree: Number of edges of a node
    degree_dict = dict(G.degree(G.nodes))

    # Centrality: More important (interconnected) nodes
    sorted_degree = sorted(degree_dict.items(), key=itemgetter(1), reverse=True)
    most_connected = sorted_degree[0][1]
    less_connected = sorted_degree[-1][1]

    # Number of clusters (nodes)
    total_nodes  = len(G)

    # Number of edges ot total og all edge weights
    total_edges = G.size()

    # return name, total_nodes, total_edges, density, triadic_closure, most_connected, less_connected
    return G

def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c)

files = glob.glob('/home/laura/Projects/ovrf-review/data/dot_plots_all/*.dot')

# with open('summary_stats.csv', 'w') as out_f:
#     w = csv.writer(out_f)
#     w.writerow(['family', 'total_nodes', 'total_edges', 'density', 'triadic_closure', 'most_connected', 'less_connected'])
#
#     for file in files:
#         base=os.path.basename(file)
#         name = os.path.splitext(base)[0]
#         row = generate_statistics(name, file)
#         w.writerow(row)


# file = '/home/laura/Projects/ovrf-review/data/adenoviridae/plot_opt_adeno_2.dot'
# G = nx.Graph(nx.drawing.nx_pydot.read_dot(file))

for file in files:
    base=os.path.basename(file)
    name = os.path.splitext(base)[0]
    G = generate_statistics(name, file)
    H = nx.connected_components(G)
    nx.draw(G)
