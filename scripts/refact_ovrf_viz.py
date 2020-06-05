# Visualize overlapping reading frames
from graphviz import Digraph

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
        total_edges = {}
        for other_cluster in clusters:
            if other_cluster.cluster != self.cluster:
                total_edges[other_cluster] = self.compare_clusters(other_cluster)

        self.cluster_edges = total_edges
        return total_edges



def get_info(handle):
    """
    Creates a list of protein objects from cluster analysis and a list of genomes
    ex: "NC_015932.E1A.489:1000;1084:1305",1
    """
    with open(handle) as file:
        _ = next(file)
        protein_list = []
        genomes_names = []
        clusters = []
        for line in file:

            # Extract info from the line
            line = line.strip().split(',')
            if len(line) != 2:  # Some proteins have a comma inside their info
                protein_info = line[0] + line[1]
            else:
                protein_info = line[0]
            info = protein_info.strip('\"').split('.')
            cluster = line[-1]  # Cluster is the last element in the line
            genome = info[0]
            name = info[1]
            location = info[-1]

            # Find start and end positions
            loc = location.strip().split(';')
            # Create as many proteins as splicing fragments and store them on the protein list
            for splice in loc:
                start, end = splice.split(':')
                new_prot = Protein(name, genome, location, int(start), int(end), cluster)
                protein_list.append(new_prot)

            # Store genomes
            if genome not in genomes_names:
                genomes_names.append(genome)
            if cluster not in clusters:
                clusters.append(cluster)

    return protein_list, genomes_names, clusters

handle2 = "../data/clusters_hc_k32.csv"
protein_list, genomes_names, clusters_names = get_info(handle2)
genome_list = [Genome(name, protein_list) for name in genomes_names]
cluster_list = [Cluster(name, protein_list) for name in clusters_names]


# Get edges of cluster  for all clusters
for cluster in cluster_list:
    cluster.get_all_edges(cluster_list)

# Create plot
dot = Digraph("New_Dot")
for cluster in cluster_list:
    for other_cluster, count in cluster.cluster_edges.items():
        edges = []
        pair = (cluster, other_cluster)
        if (pair not in edges and count != 0):
            edges.append(pair)
            dot.edge(cluster.cluster, other_cluster.cluster, label=None, penwidth=str(count/10))

print(dot.source)
dot.view()
