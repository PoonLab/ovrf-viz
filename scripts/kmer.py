from gotoh2 import *
import argparse
import csv
import math
import re
import sys

pat = re.compile('^"([^,]+),([^,]+),(.+),(1|-1),([:;0-9]+)"$')

parser = argparse.ArgumentParser(
    description="Calculate k-mer distance matrix for all protein sequences "
                "in a FASTA file."
)
parser.add_argument('-infile', type=argparse.FileType('r'),
                    help='input, FASTA file with protein sequences')
parser.add_argument('-outfile', type=argparse.FileType('w'),
                    help='output, distance matrix as CSV with accession and '
                         'gene name as column names.')
parser.add_argument('-header', type=argparse.FileType('w'),
                    help='output, protein header info as CSV')
parser.add_argument('-k', type=int, default=3, help='k-mer length (default 3)')
args = parser.parse_args()


def kmer(seq, k=3):
    d = {}
    for i in range(len(seq)-k):
        trip = seq[i:(i+k)]
        if trip not in d:
            d.update({trip: 0})
        d[trip] += 1
    return d


def kdist(k1, k2):
    d11, d12, d22 = 0, 0, 0
    for km in set(k1.keys()).union(set(k2.keys())):
        c1 = k1.get(km, 0)
        c2 = k2.get(km, 0)
        d11 += c1 * c1
        d12 += c1 * c2
        d22 += c2 * c2

    return d12 / math.sqrt(d11 * d22)

print('k={}'.format(args.k))

writer = csv.writer(args.header)
writer.writerow(['desc', 'accession', 'gene.name', 'is.forward', 'coords'])

# parse FASTA input
labels = []
seqs = []
for h, s in iter_fasta(args.infile):
    m = pat.findall(h)
    desc, accno, gene_name, directn, coords = m[0]
    writer.writerow([
        desc, accno, gene_name,
        'TRUE' if directn == '1' else 'FALSE', coords
    ])
    labels.append('{}.{}.{}'.format(accno, gene_name, coords))
    seqs.append(s)


# pre-calculate k-mer counts
kmers = {}
for i, seq in enumerate(seqs):
    kmers.update({labels[i]: kmer(seq, k=args.k)})

# calculate and write distance matrix
outstr = None
for i, l1 in enumerate(labels):
    if i > 0:
        sys.stdout.write('\b' * len(outstr))

    outstr = '{} / {}'.format(i, len(labels))
    sys.stdout.write(outstr)
    sys.stdout.flush()

    args.outfile.write('"{}"'.format(l1))
    for l2 in labels:
        if l1 == l2:
            # exact match
            args.outfile.write(',1.0')
        else:
            d = kdist(kmers[l1], kmers[l2])
            args.outfile.write(',{:1.9f}'.format(d))

    args.outfile.write('\n')

sys.stdout.write('\n')
