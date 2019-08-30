from scraper import *
from csv import DictReader, DictWriter
from time import sleep
import sys

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
if type(last_accn) is list:
    last_accn = last_accn[0]


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
        print(row)

        accn = row['Accession']
        if type(accn) is list:
            accn = accn[0]

        if accn == last_accn:
            restart = True
            continue

        if not restart:
            continue

        gid = retrieve_gid(accn)
        record = retrieve_record(gid)
        sleep(1)

        taxonomy = record.annotations['taxonomy']
        print(taxonomy)

        family = None
        for term in taxonomy:
            if term.endswith('viridae'):
                family = term

        if family is None:
            family = 'unclassified {}'.format(taxonomy[-1])

        print('{}\n'.format(family))
        row['Family'] = family
        sleep(1)

    if restart:
        writer.writerow(row)
        handle.flush()
