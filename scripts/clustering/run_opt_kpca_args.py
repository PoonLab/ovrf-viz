from glob import glob
from kmer import kmer_dist

import json
import os
import sys
import csv
import re
import argparse
import subprocess

def get_args(parser):
    parser.add_argument(
        'clust', help = 'location of family folders from where to grab protein data'
        )
    parser.add_argument(
        'headers', help = 'location of family folders from where to grab protein data'
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
    clust = glob(args.clust)
    head = glob(args.headers)
    
    print(f'Number of files: {len(clust)}')
    for f_path in clust:
        # c_file = open(f_path)
        name = os.path.basename(f_path).split('_')[0]
        print(name)
        matching = [h_file for h_file in head if name in h_file][0]
        subprocess.call(['Rscript', 'opt_kpca_args.R', f_path, matching, args.outpath])
