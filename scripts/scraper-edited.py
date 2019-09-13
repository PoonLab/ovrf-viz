from Bio import Entrez, SeqIO, SeqFeature
from csv import DictWriter
from time import sleep

Entrez.email = 'bioinfo@uwo.ca'
#'bioinfo@uwo.ca'
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

def retrieve_annotations(record):
    an = record.annotations
    an_dict = {'Topology':an['topology'], 'Taxonomy':an['taxonomy'], 'Molecule type': an['molecule_type'] }
    return an_dict

def main():

    #orffile = open('orfs.csv', 'w')
    #orffile.write('accno,product,strand,coords,aaseq\n')

    handle2 = open('../data/virus_with_info.csv', 'w')
    virus_info_file = DictWriter(handle2,
                                     delimiter = ',',
                                     extrasaction='ignore',
                                     fieldnames = [
                                         'Genome', 'Accession', 'Source information',
                                         'Topology', 'Taxonomy', 'Molecule type' ,
                                         'Genome length', 'Number of proteins', 'Host',
                                         'Date completed', 'Date updated'
                                     ])
    virus_info_file.writeheader()

    for row in parse_table(tbl):
        if row['Number of segments'] == '-':
            # unsegmented genome
            accnos = [row['Accession']]
        else:
            accnos = row['Accession']

        gid = retrieve_gid(accnos)
        record = retrieve_record(gid)
        annotation = retrieve_annotations(record)
        final_row = {**row, **annotation}
        virus_info_file.writerow(final_row)
        handle2.flush()

        print(row['Genome'], gid)  # track progress
        sleep(1)

if __name__ == "__main__":
    main()
