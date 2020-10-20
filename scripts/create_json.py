## Create JSON files for visualization of overlapping reading frames
import csv
import argparse
import json
import seaborn as sns

from collections import Counter

def get_args(parser):
    parser.add_argument(
        'file',
        help='Path to file containing cluster output in csv format'
    )
    parser.add_argument(
        '--outfile', default = None, help = 'Path to dot file'
        # With the structure: "","desc","accession","gene.name","is.forward","coords","clusters"
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


def main():
    parser = argparse.ArgumentParser(
        description = "Create DOT file for cluster analysis"
    )
    args = get_args(parser)
    print(args)
    file = args.file

    # Set information
    protein_list, genomes_names, clusters_names = get_info(file)
    genome_list = [Genome(name, protein_list) for name in genomes_names]
    cluster_list = [Cluster(name, protein_list) for name in clusters_names]

    # Define node colors
    pal = sns.color_palette(palette="husl", n_colors=len(cluster_list))
    colors = pal.as_hex()

    json_dict = {"nodes":[], "links":[], "ovrfs":[]}
    for cluster_obj in cluster_list:
        # Find most common protein name
        comp_name = Counter([str(protein) for protein in cluster_obj.proteins]).most_common(1)[0][0]
        prot_name = comp_name.replace(' protein', '')
        #print(prot_name)

        cluster = str(cluster_obj)
        start = f"start{str(cluster)}"
        end =  f"end{str(cluster)}"
        cluster_size = len(cluster_obj.proteins)
        json_dict["nodes"].append({
                                    "id": start, 
                                    "group": int(cluster), 
                                    "size": cluster_size, 
                                    "color": colors[int(cluster)-1]
                                })

        json_dict["nodes"].append({
                                    "id": end, 
                                    "group": int(cluster), 
                                    "size": cluster_size, 
                                    "color": colors[int(cluster)-1]
                                })
        
        # Self edges (between start and end nodes from the same cluster)
        json_dict["links"].append({
                                    "source": start, 
                                    "target": end, 
                                    "count": cluster_size,
                                    "className": "self",
                                    "protName":prot_name,
                                    "color": colors[int(cluster)-1]
                                    })                         

        # Create edges for adjacent proteins
        for adj_cluster_obj, count in cluster_obj.adjacent_clust.items():
            adj_cluster = str(adj_cluster_obj)
            json_dict["links"].append({
                                        "source": end, 
                                        "target": f"start{str(adj_cluster)}", 
                                        "count": count, 
                                        "className": "adj"
                                    })

        # Create edged for overlapping proteins
        for overlap_cluster_obj, count in cluster_obj.overlapping_clust.items():
            overlap_cluster = str(overlap_cluster_obj)     
            json_dict["ovrfs"].append({
                                        "source": end, 
                                        "target": f"start{str(overlap_cluster)}", 
                                        "count": count, 
                                        "className": "overlap"
                                    })

    with open(f"{args.outfile}.json", "w") as outfile:
        json.dump(json_dict, outfile, indent=4)
        


if __name__ =='__main__':
    main()
