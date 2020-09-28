from Bio import Entrez, SeqIO, SeqFeature
from csv import DictWriter
from time import sleep

Entrez.email = 'bioinfo@uwo.ca'
#accn = 'NC_015932'
tbl = '../data/taxid10239.tbl'
#tbl = 'next.tbl'

def parse_table(tbl):
    handle = open(tbl)
    _ = next(handle)
    header = next(handle).strip('\n').split('\t')
    keywords = list(map(lambda x: x.strip('"'), header))

    family = None
    row = None

    for line in handle:
        delimiter = '\t' if '\t' in line else '     '
        values = line.strip('\n').split(delimiter)

        if len(values) == 1:
            # title row
            if values[0].startswith(' '):
                # extra annotation, ignore
                continue

            if row is not None:
                row['Family'] = family
                yield row
                row = None

            family = values[0]
            continue

        values = line.strip(' \n').split(delimiter)
        if delimiter == '\t':
            if row is not None:
                # output previous entry
                row['Family'] = family
                yield row

            # prepare next entry
            row = dict(zip(keywords, values))
            if row['Accession'] == '-':
                row['Accession'] = []
                continue
        else:
            # segment
            try:
                genome, length, accno = values[:3]
            except:
                print(values)
                raise
            row['Accession'].append(accno)

    row['Family'] = family
    yield(row)  # last entry



def retrieve_gid(accn):
    handle = Entrez.esearch(db='nucleotide', term=accn, retmax=1)
    response = Entrez.read(handle)
    if response['Count'] == '0':
        # retry query
        handle = Entrez.esearch(db='nucleotide', term=accn, retmax=1)
        response = Entrez.read(handle)
        if response['Count'] == '0':
            return None

    return response['IdList'][0]

def retrieve_record(gid):
    handle = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id=gid)
    gb = SeqIO.parse(handle, format='genbank')
    return next(gb)

def retrieve_CDS(record):
    """
    Analyze features in Genbank record to extract (1) the number of coding
    regions (CDS),
    :param record:
    :return:
    """
    cds = [feat for feat in record.features if feat.type=='CDS']
    for cd in cds:
        q = cd.qualifiers
        parts = []
        for part in cd.location.parts:
            parts.append((part.start, part.end))
        locus = q.get('locus_tag', '')
        product = q.get('product', [''])
        aaseq = q.get('translation', [''])
        yield locus, product, cd.strand, parts, q['codon_start'], aaseq


def main():
    # handle = open('../data/viruses.csv', 'w')
    # outfile = DictWriter(handle,
    #                     delimiter=',',
    #                     extrasaction='ignore',
    #                     fieldnames=[
    #                         'Family', 'Genome', 'Source information',
    #                         'Accession', 'Date completed', 'Date updated',
    #                         'Genome length', 'Number of proteins', 'Host'
    #                     ])
    # outfile.writeheader()

    orffile = open('orfs_with_nt.csv', 'w')
    orffile.write('accno,product,strand,coords,aaseq\n')

    for row in parse_table(tbl):
        if row['Number of segments'] == '-':
            # unsegmented genome
            accnos = [row['Accession']]
        else:
            accnos = row['Accession']

        # outfile.writerow(row)
        # handle.flush()

        for accn in accnos:
            gid = retrieve_gid(accn)
            if gid is None:
                print('Warning, failed to retrieve gid for {}'.format(accn))
                continue

            print(row['Genome'], gid)  # track progress
            sleep(1)  # avoid spamming the server

            record = retrieve_record(gid)

            for locus, product, strand, parts, start, aaseq in retrieve_CDS(record):
                parts_str = ';'.join('{}:{}'.format(p[0].position, p[1].position) for p in parts)
                orffile.write('{},"{}",{},{},{}\n'.format(accn, product[0], strand, parts_str, aaseq[0]))

            sleep(1)



if __name__ == "__main__":
    main()