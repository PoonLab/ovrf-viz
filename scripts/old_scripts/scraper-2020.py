from Bio import Entrez, SeqIO, SeqFeature
from csv import DictWriter
from time import sleep

Entrez.email = 'lmuoz@uwo.ca'
#'bioinfo@uwo.ca'
#accn = 'NC_015932'
#tbl = '/home/lmunoz/Projects/ovrf-review/dataset_2020/taxid10239.nbr'
tbl = '/home/laura/Downloads/taxid.nbr'

def parse_table(tbl):
    """
    Yield every row on the table parsed for downstream analysis 
    """
    handle = open(tbl)
    _ = next(handle)
    header = next(handle).strip('\n').split('\t')
    keywords = list(map(lambda x: x.strip('"'), header))
    print(keywords)

    row = None
    accno = ''

    for line in handle:
        delimiter = '\t' if '\t' in line else '     '
        values = line.strip('\n').split(delimiter)

        #Prepare next entry
        row = dict(zip(keywords, values))
        yield(row)


def retrieve_gid(accn):
    """
    Get genome ID for an accession number
    """
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
    try:
        handle = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id=gid)
        gb = SeqIO.parse(handle, format='genbank')
        return next(gb)
    except:
        print("Failed, trying again")
        return retrieve_record(gid)

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
    an_dict = {'Length': len(record.seq) , 'Topology':an['topology'],
                'Taxonomy':an['taxonomy'], 'Molecule': an['molecule_type'],
                'Proteins': len(record.features)
                }
    return an_dict

def main():

    #orffile = open('orfs.csv', 'w')
    #orffile.write('accno,product,strand,coords,aaseq\n')
    handle = open('virus_with_info.csv', 'w')
    outfile = DictWriter(handle, delimiter=',',
                                 extrasaction='ignore',
                                 fieldnames=[
                                     'Representative', 'Neighbor', 'Host', 'Selected lineage',
                                     'Taxonomy name', 'Length', 'Proteins',
                                     'Topology', 'Taxonomy', 'Molecule', 'Segment name"i'
                                 ])
    outfile.writeheader()

    orffile = open('orfs.csv', 'w')
    orffile.write('accno,product,strand,coords,start_codon\n')

    refsecs = []
    for row in parse_table(tbl):
        accnos = [row['Representative']]

        for accn in accnos:
            if accn not in refsecs:
                gid = retrieve_gid(accnos)

                if gid is None:
                    print('Warning, failed to retrieve gid for {}'.format(accn))
                    continue

                print(row['Taxonomy name'], gid, accn)
                sleep(1)

                record = retrieve_record(gid)
                annotation = retrieve_annotations(record)
                final_row = {**row, **annotation}
                outfile.writerow(final_row)
                handle.flush()
                refsecs.append(accn)

                for locus, product, strand, parts, start_codon, aaseq in retrieve_CDS(record):
                    parts_str = ';'.join('{}:{}'.format(p[0].position, p[1].position) for p in parts)
                    orffile.write('{},"{}",{},{},{}\n'.format(accn, product[0], strand, parts_str, start_codon))

                sleep(1)


if __name__ == "__main__":
    main()
