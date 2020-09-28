# Create visualization for OVRFs

from graphviz import Digraph
from itertools import tee

handle = "../data/test_prot.csv"

def create_prot_dict(handle):
    """
    Creates a list of dictionaries with protein info sorted by the start position of the protein
    """
    with open(handle) as file:
        _ = next(file)
        protein_list =[]
        for line in file:
            name, start, end = line.strip().split(',')
            protein_dict = {'name':name ,'start': int(start),'end': int(end)}
            protein_list.append(protein_dict)

    prot_sorted = sorted(protein_list, key = lambda i: i['start'])
    return prot_sorted

def get_overlaping_proteins(protein_list):
    """
    Check which proteins overlap
    """
    updated_list = protein_list
    for i in range(len(updated_list)):
        current_prot = updated_list[i]
        current_prot.update({'overlap':[]})
        for j in range(i+1, len(updated_list)):
            other_prot = updated_list[j]
            print(current_prot['name'], other_prot['name'])
            # Check if the start site is between the range of the other protein
            if other_prot['start'] <= current_prot['start'] <= other_prot['end']:
                current_prot['overlap'].append(other_prot['name'])
            # Check if the end site is between the range of the other protein
            elif other_prot['start'] <= current_prot['end'] <= other_prot['end']:
                current_prot['overlap'].append(other_prot['name'])

    return updated_list

def pairwise(iterable):
    """
    Retrieve pairs on an iterable object
    """
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def create_diagraph(prtoteins):
    """
    Creates a diagram of proteins considering their location and overlap
    """

    dot = Digraph('Protein_diagraph')
    for prot in proteins:
        dot.node(prot['name'])

    for i in range(len(proteins)):
        current_prot = proteins[i]
        current_overlaps = current_prot['overlap']
        no_overlap = False
        for j in range(i+1,len(proteins)):
            follow_prot = proteins[j]
            follow_name = follow_prot['name']
            if  follow_name in current_overlaps:
                dot.edge(current_prot['name'], follow_name)
                print("There is an overlap")
            else:
                dot.edge(current_prot['name'], follow_name)
                print("No overlap, I will jump")
                no_overlap = True
                break
        if no_overlap == True:
            continue

    return dot


protein_list = create_prot_dict(handle)
proteins = get_overlaping_proteins(protein_list)
dot = create_diagraph(proteins)
print(dot.source)
#dot.render()
dot.view()
