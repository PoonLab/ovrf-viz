# Visualize overlapping reading frames
from graphviz import Digraph
import argparse
import math
import csv

# import numpy as np
# from matplotlib.path import Path
# from matplotlib.patches import PathPatch
# import matplotlib.pyplot as plt
# import matplotlib.cm as cm



import numpy as np
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def get_args(parser):
    parser.add_argument(
        'file',
        help='Path to file containing cluster output in csv format'
    )
    parser.add_argument(
        '--outfile', default = None, help = 'Path to dot file'
    )
    parser.add_argument(
        '--edge_count', default = None, help = 'Minimum number of genomes required to draw an edge'
    )

    return parser.parse_args()


class Protein:
    """
    Creates protein objects with information
    """

    def __init__(self, name, genome, location, start, end, cluster):
        """
        Stores relevant information about the protein and its location on the genome
        """
        self.name = name
        self.genome = genome
        self.location = location
        self.cluster = cluster
        self.adjacent_proteins = []
        self.start = start
        self.end = end
        self.overlaps = []

    def __repr__(self):
        return self.name


class Genome:
    """
    Creates list of proteins in the genome, sorted by position
    """

    def __init__(self, accno, all_proteins):
        self.accno = accno
        self.proteins = self.get_my_proteins(all_proteins)
        self.get_adjacent_proteins()  # Get adjacent proteins for every protein in the genome

    def __repr__(self):
        return self.accno

    def get_my_proteins(self, all_proteins):
        """
        Create a list with all proteins in the genome
        """

        proteins = [protein for protein in all_proteins if protein.genome == self.accno]
        return sorted(proteins, key=lambda x: x.start)

    def get_adjacent_proteins(self):
        """
        For every protein in the genome, store adjacent proteins
        """

        for i in range(len(self.proteins)):
            current_prot = self.proteins[i]
            for j in range(i+1, len(self.proteins)):
                other_prot = self.proteins[j]
                current_prot.adjacent_proteins.append(other_prot)  # Store following protein
                # Check if proteins overlap
                # Check if the start site is between the range of the other protein
                if other_prot.start <= current_prot.start <= other_prot.end:
                    current_prot.overlaps.append(other_prot)
                # Check if the end site is between the range of the other protein
                elif other_prot.start <= current_prot.end <= other_prot.end:
                    current_prot.overlaps.append(other_prot)
                else:  # If proteins don't overlap
                    break


class Cluster:
    """
    Creates a list of proteins in the cluster
    """

    def __init__(self, cluster, all_proteins):
        self.cluster = cluster
        self.proteins = self.get_my_proteins(all_proteins)
        self.edges = self.get_edge_list()

    def __repr__(self):
        return self.cluster

    def get_my_proteins(self, all_proteins):
        """
        Create a list with all proteins in the cluster
        """
        return [protein for protein in all_proteins if protein.cluster == self.cluster]

    def get_edge_list(self):
        """
        Get all edges for proteins in the cluster
        """
        edges = []
        for prot in self.proteins:
            edges.extend(prot.adjacent_proteins)
        return edges

    def compare_clusters(self, other_cluster):
        """
        Count the number of times that a edge is created between a protein in current cluster and protein in other cluster
        """
        edges_count = 0
        for protein in other_cluster.proteins:  # Loop trough proteins in the other cluster
            if protein in self.edges:
                edges_count += self.edges.count(protein)

        return edges_count

    def get_all_edges(self, clusters):
        """
        Get edges between current cluster and every other cluster on the data set
        :return: dictionary keyed by cluster name with the count of edges between them
        """
        total_edges = {}
        for other_cluster in clusters:
            if other_cluster.cluster != self.cluster:
                total_edges[other_cluster] = self.compare_clusters(other_cluster)

        self.cluster_edges = total_edges
        return total_edges



def get_info(handle):
    """
    Creates a list of protein objects from cluster analysis and a list of genomes
    ex: "24","Bat adenovirus 2","NC_015932","fiber",TRUE,"26819:28493",23
    input labels: 'desc', 'accession', 'gene.name', 'is.forward', 'coords'
    """
    with open(handle) as file:
        records = csv.DictReader(file)
        protein_list = []
        genomes_names = []
        clusters = []
        for row in records:
            # Find start and end position for each protein
            loc = row['coords'].strip().split(';')
            # Create as many proteins as splicing fragments and store them on the protein list
            for splice in loc:
                start, end = splice.split(':')
                #new_prot = Protein(name, genome, location, int(start), int(end), cluster)
                new_prot = Protein(row['gene.name'], row['accession'], row['coords'], int(start), int(end), row['clusters'])
                protein_list.append(new_prot)

            # Store genomes
            if row['accession'] not in genomes_names:
                genomes_names.append(row['accession'])
            # Store clusters
            if row['clusters'] not in clusters:
                clusters.append(row['clusters'])

    return protein_list, genomes_names, clusters


def label(xy, width, height, text):
    """
    Create labels with protein names
    """
    x = (xy[0]+width/2)
    y = (xy[1]+height/2)
    plt.text(x, y, text, ha="center", family='sans-serif', size=5.5)


def genome_plot(colors, genome_list):
    """
    Use matplotlib to create a plot of each genome, where each protein is colored according to the cluster they belong to
    """

    fig= plt.figure() # Figure object
    for i in range(len(genome_list)):
        gen = genome_list[i]
        prots = gen.proteins
        ax = fig.add_subplot(len(genome_list), 1, i+1)

        patches = []
        # Each protein is represented as a rectangle in the plot
        for prot in prots:
            start = prot.start
            end = prot.end
            bottom = 1
            top = 100
            width = end - start
            height = top - bottom
            rect1 = mpatches.Rectangle((start,bottom), width, height,  color = colors[int(prot.cluster)-1])
            patches.append(rect1)
            #label((start,bottom), width, height, text=prot.name)  # Optional to print name of the protein

        collection = PatchCollection(patches, match_original=True)
        ax.add_collection(collection)
        # TO DO: Find a better way to plot the name of the genome to have them aligned
        ax.text(-3500, 1, gen, fontsize=6, ha='center', va='center')
        plt.axis('auto')  # Automatically adjust the size of the ax to the size of the actual plot
        plt.axis('off')  # Don't show the axis label

    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description = "Create DOT file for cluster analysis"
    )
    args = get_args(parser)
    print(args)
    file = args.file
    if args.edge_count != None:
        min_edge = int(args.edge_count)
    else:
        min_edge = 1

    # Set information
    protein_list, genomes_names, clusters_names = get_info(file)
    genome_list = [Genome(name, protein_list) for name in genomes_names]
    cluster_list = [Cluster(name, protein_list) for name in clusters_names]

    # Color pallette used in the Adenoviridae clustering method in R
    colors = ["#F8766D", "#F17D51", "#E98429", "#DF8B00", "#D49200",
     "#C89800", "#BA9E00", "#AAA300", "#97A800", "#82AD00", "#67B100",
     "#3FB500", "#00B929", "#00BC4F", "#00BE6B", "#00BF82", "#00C097",
      "#00C1AA", "#00C0BC", "#00BECD", "#00BBDC", "#00B7E9", "#00B1F4",
       "#00AAFE", "#30A2FF", "#7299FF", "#988FFF", "#B584FF", "#CC7AFF", "#DE70F9",
        "#EC68EE", "#F663E1", "#FD61D2", "#FF61C1", "#FF64AF", "#FF699B", "#FD6F85"]

    # Get edges of cluster  for all clusters
    for cluster in cluster_list:
        cluster.get_all_edges(cluster_list)

    # Create plot
    dot = Digraph(comment='Cluster plot')
    dot.graph_attr['outputorder'] = 'endgesfirst'

    for cluster in cluster_list:
        #print(cluster, len(cluster.proteins))
        size = len(cluster.proteins)
        dot.node(cluster.cluster, label=None, fixedsize="true", width=str(math.sqrt(size)/3), height=str(math.sqrt(size)/3),
                 fontsize=str(85), style='filled', color=colors[int(cluster.cluster)-1], fontname = 'Courier-Bold')
        for other_cluster, count in cluster.cluster_edges.items():
            edges = []
            pair = (cluster, other_cluster)
            if (pair not in edges and count > min_edge):

                edges.append(pair)
                #Create edge
                if count >= 10:  # Color the most frequent edges
                    dot.edge(cluster.cluster, other_cluster.cluster, label=None, penwidth=str(count),
                     color='grey60', arrowsize = str(0.01), len = str(10))
                else:  # Light color for the least frequent edges
                    dot.edge(cluster.cluster, other_cluster.cluster, label=None, penwidth=str(count),
                     color='grey83', arrowsize = str(0.01))

    dot.render(filename="{}.dot".format(args.outfile))

    # Create genome plot
    genome_plot(colors, genome_list)

if __name__ =='__main__':
    main()
