# Cluster homologous proteins from distance matrix created using kmer.py
from email import header
import numpy as np
import pandas as pd

from sklearn.manifold import TSNE

import matplotlib.pyplot as plt

import seaborn as sns

import argparse

def get_args(parser):
    # Required arguments
    parser.add_argument(
        'km',
        help='Path to the file containing csv file of kmer distance between proteins'
    )

    parser.add_argument(
        'lab',
        help='Path to header with labels of distance matrix'
    )    

    return parser.parse_args()


def main():

    parser = argparse.ArgumentParser(
        description='From list of accession numbers, retrieve amino acid sequences as fasta file'
    )

    args = get_args(parser)
    labels = pd.read_csv(args.lab)
    km = pd.read_csv(args.km, header=None, index_col=0)
    # print(km.dtypes) km[rows, columns]
    # subset=km.iloc[: ,1:]
    print(km.shape, labels.shape)
    print( km.head())
    # print(km.index)  # index is rownames, header is colnames
    accn = km.index.str.split(".").str[0]

    

    tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
    tsne_results = tsne.fit_transform(km)

    km['tsne-2d-one'] = tsne_results[:,0]
    km['tsne-2d-two'] = tsne_results[:,1]
    km['accn'] = accn
    
    plt.figure(figsize=(16,10))
    sns.scatterplot(
        x="tsne-2d-one", y="tsne-2d-two",
        hue="accn",  # Color by
        palette=sns.color_palette("hls", 10),
        data=km,
        legend="full",
        alpha=0.3
    )
    plt.show()

if __name__ == "__main__":
    main()