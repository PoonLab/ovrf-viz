from gotoh2 import *
import argparse
import csv
import math
import re
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate k-mer distance matrix for all protein sequences "
                    "in a FASTA file."
    )
    parser.add_argument('infile', type=argparse.FileType('r'),
                        help='input, FASTA file with protein sequences')
    parser.add_argument('outfile', type=argparse.FileType('w'),
                        help='output, distance matrix as CSV with accession and '
                             'gene name as column names.')
    parser.add_argument('-header', type=argparse.FileType('w'), default=None,
                        help='output (optional), protein header info as CSV')
    parser.add_argument('--kernel', action='store_true',
                        help='optional, use string kernel (d2s) instead of '
                        'default intersection distance.')
    return parser.parse_args()


def kmer(seq):
    """
    Calculate word counts for word lengths varying from k={1,2,3}
    :param seq:  str, protein sequence
    :return:  dict, counts keyed by word - absent entries imply zero
    """
    d = {}
    for i in range(len(seq)):
        for k in [1, 2, 3]:
            if (len(seq)-1 - i) < k:
                break
            word = seq[i:(i+k)]
            if word not in d:
                d.update({word: 0})
                d[word] += 1
    return d


def kdist(k1, k2):
    """
    string kernel / d2s distance
    :param k1:  kmer counts for sequence 1
    :param k2:  kmer counts for sequence 2
    :return:  pairwise distance
    """
    d11, d12, d22 = 0, 0, 0
    for km in set(k1.keys()).union(set(k2.keys())):
        c1 = k1.get(km, 0)
        c2 = k2.get(km, 0)
        d11 += c1 * c1
        d12 += c1 * c2
        d22 += c2 * c2

    return d12 / math.sqrt(d11 * d22)


def intersection(k1, k2):
    """
    Intersection distance
    :param k1:  kmer counts for sequence 1
    :param k2:  kmer counts for sequence 2
    :return:  float, pairwise distance
    """
    res = 0
    for km in set(k1.keys()).intersection(set(k2.keys())):
        res += 2 * min(k1[km], k2[km])
    return res / (sum(k1.values()) + sum(k2.values()))


def kmer_dist(infile, outfile, header, kernel):
    """
    :param infile:  open file stream to read FASTA
    :param outfile:  open file stream to write distances as CSV
    :param header:  open file stream to write header info as CSV
    :param kernel:  if True, calculate kernel distance instead of
                    intersection distance
    """

    # regular expression to parse sequence headers
    pat = re.compile('^"([^,]+),([^,]+),(.+),(1|-1),([:;0-9]+)"$')

    # prepare output file
    writer = None
    if header:
        writer = csv.writer(header)
        writer.writerow(['desc', 'accession', 'gene.name', 'is.forward', 'coords'])

    # parse FASTA input
    labels = []
    seqs = []
    for h, s in iter_fasta(infile):
        m = pat.findall(h)
        desc, accno, gene_name, directn, coords = m[0]
        if writer:
            writer.writerow([
                desc, accno, gene_name,
                'TRUE' if directn == '1' else 'FALSE', coords
            ])
        labels.append('{}.{}.{}'.format(accno, gene_name, coords))
        seqs.append(s)


    # pre-calculate k-mer counts
    kmers = {}
    for i, seq in enumerate(seqs):
        kmers.update({labels[i]: kmer(seq)})

    # calculate and write distance matrix
    outstr = None
    for i, l1 in enumerate(labels):
        if i > 0:
            sys.stdout.write('\b' * len(outstr))

        outstr = '{} / {}'.format(i, len(labels))
        sys.stdout.write(outstr)
        sys.stdout.flush()

        outfile.write('"{}"'.format(l1))
        for l2 in labels:
            if kernel:
                d = 1. if l1 == l2 else kdist(kmers[l1], kmers[l2])
            else:
                # default
                d = intersection(kmers[l1], kmers[l2])


            outfile.write(',{:1.9f}'.format(d))

        outfile.write('\n')

    sys.stdout.write('\n')


if __name__ == '__main__':
    args = parse_args()
    kmer_dist(infile=args.infile, outfile=args.outfile, header=args.header,
             kernel=args.kernel)
