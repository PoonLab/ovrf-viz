# A single script to get all protein fasta files from all families, measure kmer distance, and create clustering results

from glob import glob
from kmer import kmer_dist

import os
import sys
import argparse
import subprocess

def get_args(parser):
    parser.add_argument(
        'files', help = 'location of family folders from where to grab protein data'
        )
    parser.add_argument(
        'outpath',
        help = 'Path where outputs should be stored'
    )
    return parser.parse_args()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='From multi fasta with proteins for virus families,\
         create clusters'
    )

    args = get_args(parser)
    fa_files = glob(args.files)
    
    print(f'Number of files: {len(fa_files)}')
    for f_path in fa_files:
        f = open(f_path)
        name=os.path.basename(f_path)
        virus=name.split('.')[0]
        out = open(f'{args.outpath}{virus}_kmer.csv', 'w')
        header = open(f'{args.outpath}{virus}_header.txt', 'w')

        print(f'Processing virus: {virus}')
        
        kmer_dist(infile=f, outfile=out, header=header, kernel=False)