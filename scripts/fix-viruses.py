from scraper import parse_table
from csv import DictReader, DictWriter

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

for row in taxid:
    writer.writerow(row)
