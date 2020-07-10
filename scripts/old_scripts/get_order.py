from Bio import Entrez
from csv import DictReader, DictWriter
from time import sleep

Entrez.email = 'apoon42@uwo.ca'


handle = open('../data/viruses-order.csv', 'w')
writer = DictWriter(handle,
                    delimiter=',',
                    extrasaction='ignore',
                    fieldnames=[
                        'Order', 'Family', 'Genome', 'Source information',
                        'Accession', 'Date completed', 'Date updated',
                        'Genome length', 'Number of proteins', 'Host'
                    ])
writer.writeheader()


reader = DictReader(open('../data/viruses.csv'))
for row in reader:
    handle = Entrez.esearch(db='taxonomy', term='"{}"'.format(row['Genome']))
    response = Entrez.read(handle)
    taxid = response['IdList'][0]

    sleep(1)

    handle = Entrez.efetch(db='taxonomy', id=taxid)
    response = Entrez.read(handle)
    taxonomy = response[0]['Lineage'].split('; ')

    row['Order'] = None
    for term in taxonomy:
        if term.endswith('virales'):
            row['Order'] = term

    print(row['Genome'], row['Order'])

    writer.writerow(row)
    handle.flush()
    sleep(1)
