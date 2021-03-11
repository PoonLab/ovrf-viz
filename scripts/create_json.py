## Create JSON files for visualization of overlapping reading frames
import csv
import argparse
import json
import seaborn as sns

from collections import Counter
from new_viz_ovrf import Protein, Genome, Cluster

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

    #print(json_dict)
    # with open(f"{args.outfile}.json", "w") as outfile:
    #     json.dump(json_dict, outfile, indent=4)

    def encode_genome_list(z):
        if isinstance(z, list):
            return z
        elif isinstance(z, Genome):
            return { "accno": z.accno, "proteins": z.proteins }
        elif isinstance(z, Protein):
            return { "name": z.name, "start": z.start, "end": z.end, "color": colors[int(z.cluster)-1] }
        else:
            type_name = z.__class__.__name__
            raise TypeError(f"Object of type '{type_name}' is not JSON serializable") 

    with open(f"genome_{args.outfile}.json", "w") as outfile:
        json.dump(genome_list, outfile, default=encode_genome_list, indent=4)

if __name__ =='__main__':
    main()
