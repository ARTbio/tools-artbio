#!/usr/bin/env python3

import argparse
import scrublet as scr
import matplotlib.pyplot as plt
import scipy.io
import numpy as np
import pandas as pd

def Parser():
    the_parser = argparse.ArgumentParser(
        description="")
    the_parser.add_argument('--input_barcodes', action="store", type=str)
    the_parser.add_argument('--input_features', action="store", type=str)
    the_parser.add_argument('--input_matrix', action="store", type=str)
    the_parser.add_argument('--expected_doublet_rate', action="store", type=float)
    the_parser.add_argument('--doublet_threshold', action="store", type=float)
    the_parser.add_argument('--output_csv', action="store", type=str)
    the_parser.add_argument('--output_hist', action="store", type=str)
    the_parser.add_argument('--output_umap', action="store", type=str)
    args = the_parser.parse_args()
    return args


def __main__():
    args = Parser()
    # Import des donnees 10X
    counts_matrix = scipy.io.mmread(args.input_matrix).T.tocsc()
    genes = np.array(scr.load_genes(args.input_features, delimiter='\t', column=2))
    barcodes = pd.read_csv(args.input_barcodes, delimiter='\t', header=None)
    print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
    print('Number of genes in gene list: {}'.format(len(genes)))
    print('Number of barcodes: {}'.format(len(barcodes)))

    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=args.expected_doublet_rate)

    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=50)
    print('User threshold at doublet score :', args.doublet_threshold)
    scrub.call_doublets(threshold=args.doublet_threshold)

    barcodes['doublets']=pd.DataFrame(scrub.predicted_doublets_)
    barcodes['score']=pd.DataFrame(doublet_scores)
    barcodes.to_csv(args.output_csv)


    scrub.plot_histogram();
    plt.savefig(args.output_hist)
    print('Running UMAP...')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    print('Done.')
    scrub.plot_embedding('UMAP', order_points=True);
    plt.savefig(args.output_umap)

if __name__ == "__main__":
    __main__()
