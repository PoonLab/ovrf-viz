# Gather information about Baltimore Classification for every species
from ete3 import NCBITaxa
from Bio import Entrez
from csv import DictReader, DictWriter

#lineage = ncbi.get_lineage(taxid)
#ncbi.get_taxid_translator(lineage)
#out (taxid = 438782) = {1: 'root', 10239: 'Viruses', 251095: 'Nanoviridae', 251096: 'Babuvirus', 438782: 'Abaca bunchy top virus'}

Entrez.email = 'lmuoz@uwo.ca'


accn = ['NC_029899', 'NC_029898']
genome = ['Bat mastadenovirus WIV18', 'Bat adenovirus TJM']
ncbi = NCBITaxa()
for name in genome:

    handle = Entrez.esearch(db='taxonomy', term = name)
    response = Entrez.read(handle)
    taxid = response['IdList'][0]

    print(taxid)

    handle = Entrez.efetch(db='taxonomy', id = taxid)
    response = Entrez.read(handle)
    lineage = ncbi.get_lineage(taxid)
    name = ncbi.get_taxid_translator(lineage)

    print(name)