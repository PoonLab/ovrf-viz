"""
The subtitles in the .tbl file are not consistently virus family names.
Use the NCBI taxonomy database to query the virus name and retrieve the
family annotation.
"""

from scraper import *
from csv import DictReader, DictWriter
from time import sleep
import sys
from Bio import Entrez
import xml.etree.ElementTree as ET

tbl = '../data/taxid10239.tbl'
taxid = parse_table(tbl)

classif = [
    ('viria', 'realm'),
    ('vira', 'subrealm'),
    ('viriae', 'kingdom'),
    ('virites', 'subkingdom'),
    ('viricota', 'phylum'),
    ('viricotina', 'subphylum'),
    ('viricetes', 'class'),
    ('viricetidae', 'subclass'),
    ('virales', 'order'),
    ('virineae', 'suborder'),
    ('viridae', 'family'),
    ('virinae', 'subfamily')
]

#reader = DictReader(open('../data/viruses.csv'))

# check current file for restart position
reader = DictReader(open('../data/viruses.csv'))
for row in reader:
    pass  # scan through to last row

last_accn = row['Accession']
if last_accn.startswith('['):
    # parse as list
    last_accn = last_accn.strip('"[]"').split(',')[0].strip("'")


handle = open('../data/viruses-fixed.csv', 'w')
writer = DictWriter(handle,
                    delimiter=',',
                    extrasaction='ignore',
                    fieldnames=[
                        'Family', 'Genome', 'Source information',
                        'Accession', 'Date completed', 'Date updated',
                        'Genome length', 'Number of proteins', 'Host'
                    ])
writer.writeheader()


restart = False

for row in taxid:
    family = row['Family']
    if not family.endswith('viridae') and not family.endswith('litidae'):
        accn = row['Accession']
        if type(accn) is list:
            accn = accn[0]

        if accn == last_accn:
            restart = True
            continue

        if not restart:
            continue

        print(row)

        #gid = retrieve_gid(accn)
        handle = Entrez.esearch(db='taxonomy', term=row['Genome'])
        response = Entrez.read(handle)
        taxid = response['IdList'][0]

        sleep(1)

        handle = Entrez.efetch(db='taxonomy', id=taxid)
        response = Entrez.read(handle)
        taxonomy = response[0]['Lineage'].split('; ')

        #record = retrieve_record('{}.1'.format(accn))
        family = None
        for term in taxonomy:
            if term.endswith('viridae'):
                family = term

        if family is None:
            family = 'unclassified {}'.format(taxonomy[-1])

        print('{}\n'.format(family))
        row['Family'] = family
        sleep(2)

    if restart:
        writer.writerow(row)
        handle.flush()  # FIXME: not supported in Python 3.5?
