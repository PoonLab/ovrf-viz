from scraper import *
from csv import DictWriter
import sys

tbl = '../data/taxid10239.tbl'
taxid = parse_table(tbl)


#reader = DictReader(open('../data/viruses.csv'))

# overwrite file!!!
writer = DictWriter(open('../data/viruses.csv', 'w'),
                    delimiter=',',
                    extrasaction='ignore',
                    fieldnames=[
                        'Family', 'Genome', 'Source information',
                        'Accession', 'Date completed', 'Date updated',
                        'Genome length', 'Number of proteins', 'Host'
                    ])
writer.writeheader()


sys.exit()

for row in taxid:
    family = row['Family']

    accn = row['Accession']
    if type(accn) is list:
        accn = accn[0]
    gid = retrieve_gid(accn)
    record = retrieve_record(gid)

    taxonomy = record.annotations['taxonomy']

    nucleic = record.annotations['molecule_type']
    topology = record.annotations['topology']
    sleep(1)

    writer.writerow(row)
