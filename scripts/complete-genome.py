"""
Issue #4 - screen for complete genomes or segments
"""
from Bio import Entrez
from time import sleep
from csv import DictReader, DictWriter

Entrez.email = 'apoon42@uwo.ca'

reader = DictReader(open('../data/viruses.csv'))

handle = open('../data/complete.csv', 'w')
writer = DictWriter(handle, fieldnames=reader.fieldnames)
writer.writeheader()

def get_title(accn):
    response = Entrez.esummary(db='nucleotide', id='{}.1'.format(accn))
    return Entrez.read(response)

for row in reader:
    accn = row['Accession']
    if accn.startswith('['):
        # handle multi-accession record
        items = map(lambda x: x.strip("'"), accn.strip('[]').split(', '))
        accn = next(items)

    try:
        response = get_title(accn)
    except RuntimeError:
        # try again
        try:
            sleep(1)
            response = get_title(accn)
        except:
            print('failed to retrieve {}'.format(accn))
            continue

    title = response[0]['Title']
    is_complete = False
    if 'complete genome' in title or \
            'complete sequence' in title or \
            'genomic sequence' in title or \
            (title.startswith('Complete') and title.endswith('genome')):
        is_complete = True

    if is_complete:
        writer.writerow(row)
        handle.flush()
    else:
        print('Not complete: {} {}'.format(accn, title))

    sleep(0.5)
