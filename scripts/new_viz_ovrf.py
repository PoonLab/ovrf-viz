## Visualize overlapping reading frames
from graphviz import Digraph
import argparse
import math
import csv

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from matplotlib.path import Path
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

from wordcloud import WordCloud
from collections import Counter

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
        self.adjacent_clust, self.overlapping_clust = self.get_cluster_count()


    def __repr__(self):
        return self.cluster

    def get_my_proteins(self, all_proteins):
        """
        Create a list with all proteins in the cluster
        """
        return [protein for protein in all_proteins if protein.cluster == self.cluster]

    def get_connected_clusters(self):
        """
        Get proteins from other clusters associated with proteins in my cluster
        """
        adjacent = []  # Adjacent edge formation
        overlapping = []  # Overlapping edge formation
        for prot in self.proteins:
            overlapping.extend(prot.overlaps)
            adjacent.extend(set(prot.adjacent_proteins) -  set(prot.overlaps))

        adjacent_clust = [protein.cluster for protein in adjacent]
        overlapping_clust = [protein.cluster for protein in overlapping]

        return adjacent_clust, overlapping_clust


    def find_self_edges(self):
        print("Running the function")
        adjacent = []  # Adjacent edge formation
        overlapping = []  # Overlapping edge formation
        for prot in self.proteins:
            overlapping.extend(prot.overlaps)
            adjacent.extend(set(prot.adjacent_proteins) -  set(prot.overlaps))
        genomes = []
        for protein in overlapping:

            if protein.genome not in genomes:
                genomes.append(protein.genome)


    def get_cluster_count(self):
        """
        Creates dictionaries counting the ammount of times a cluster has an edge with current cluster for adjacent and overlapping cases
        """
        adjacent_clust, overlapping_clust = self.get_connected_clusters()

        # Create dictionary for adjacent clusters
        adj_freq = {}
        for cluster in adjacent_clust:
            adj_freq[cluster] = adjacent_clust.count(cluster)

        overlap_freq = {}
        for cluster in overlapping_clust:
            overlap_freq[cluster] = overlapping_clust.count(cluster)

        return adj_freq, overlap_freq


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
        ax.text(-700, 1, gen, fontsize=6, ha='center', va='center')
        plt.axis('auto')  # Automatically adjust the size of the ax to the size of the actual plot
        plt.axis('off')  # Don't show the axis label

    plt.show()

def wordcloud_plot(cluster_list):
    """
    Create a wordcloud plot for all proteins in each cluster
    """
    number_of_subplots = len(cluster_list)
    # TO DO: nrows and ncols have to be manually edited in order to properly distribute the plots across the Figure
    number = math.sqrt(number_of_subplots)

    # Calculate number of rows and cols
    if number % 1 == 0:
        nrows = number
        ncols = nrows
    else:
        r = round(number)
        nrows = math.trunc(number)
        ncols = nrows +1
        if r != nrows:
            nrows = r
            ncols = nrows

    # Make plot
    fig= plt.figure()
    for i in range(number_of_subplots):
        cluster = cluster_list[i]
        s = i+1
        ax = fig.add_subplot(nrows, ncols, i+1)
        ax.set_title("Cluster #" + str(cluster))
        protein_list = [str(protein) for protein in cluster.proteins]
        wordcloud_dict = Counter(protein_list)
        wordcloud = WordCloud(background_color = "white").generate_from_frequencies(wordcloud_dict)
        plt.imshow(wordcloud)
        plt.axis("off")  # No axis

    fig.tight_layout(pad=1) # Increase spacing between figures
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description = "Create DOT file for cluster analysis"
    )
    args = get_args(parser)
    print(args)
    file = args.file

    # Set the minimim ammount of connections requiered to draw an edge
    if args.edge_count != None:
        min_edge = int(args.edge_count)
    else:
        min_edge = 1

    # Set information
    protein_list, genomes_names, clusters_names = get_info(file)
    genome_list = [Genome(name, protein_list) for name in genomes_names]
    cluster_list = [Cluster(name, protein_list) for name in clusters_names]

    # Define colors for plots
    pal = sns.color_palette(palette="husl", n_colors=len(cluster_list))
    colors = pal.as_hex()

    # Create plot
    dot = Digraph(comment='Cluster plot')
    dot.graph_attr['outputorder'] = 'endgesfirst'

    for cluster in cluster_list:
        cluster_size = len(cluster.proteins)
        node_size = math.sqrt(cluster_size)/3

        # Create a node
        dot.node(cluster.cluster, label=None, fixedsize="true", width=str(node_size), height=str(node_size),
                 fontsize=str(85), style='filled', color=colors[int(cluster.cluster)-1], fontname = 'Courier-Bold')


        # Create adges for adjacent proteins
        for adj_cluster, count in cluster.adjacent_clust.items():
            if count >= min_edge:
                dot.edge(cluster.cluster, adj_cluster, label = None, penwidth = str(count),
                color = "grey76", arrowsize = str(0.01), len = str(10))

        # print("Overlapping", cluster.overlapping_clust)
        # overlap_count = 0
        # Create edges for overlapping proteins
        for overlap_cluster, count in cluster.overlapping_clust.items():
            if count >= min_edge:
                dot.edge(cluster.cluster, overlap_cluster, label = None, penwidth=str(count),
                color="#143D59", arrowsize = str(0.01), len = str(10))
                if cluster.cluster == overlap_cluster:
                    print(f"Cluster: {cluster.cluster}, Numer of self edges: {count}")



    dot.render(filename="{}.dot".format(args.outfile))

    # Create wordcloud plot
    wordcloud_plot(cluster_list)
    # Create genome plot
    genome_plot(colors, genome_list)

if __name__ =='__main__':
    main()
