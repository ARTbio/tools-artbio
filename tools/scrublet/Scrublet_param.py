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
    the_parser.add_argument('--input_dir', action="store", type=str,
                            help="Chemin vers le dossier ou se trouve les fichiers au format 10x")
    the_parser.add_argument('--expected_doublet_rate', action="store", type=float,
                            help="")
    the_parser.add_argument('--seuil', action="store", type=float,
                            help="Seuil qui va discriminer le statut de doublet des cellules selon leur score")
    the_parser.add_argument('--output_dir', action="store", type=str,
                            help="Chemin vers le dossier où les trois fichiers de sorties seront enregistrés")
    args = the_parser.parse_args()
    return args


def __main__():
    args = Parser()
    # Import des donnees 10X
    counts_matrix = scipy.io.mmread(args.input_dir + '/matrix.mtx').T.tocsc()
    genes = np.array(scr.load_genes(args.input_dir + '/features.tsv', delimiter='\t', column=2))
    barcodes = pd.read_csv(args.input_dir + 'barcodes.tsv', delimiter='\t', header=None)
    print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
    print('Number of genes in gene list: {}'.format(len(genes)))
    print('Number of barcodes: {}'.format(len(barcodes)))

    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=args.expected_doublet_rate)

    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=50)
    print('User threshold at doublet score :', args.seuil)
    scrub.call_doublets(threshold=args.seuil)

    barcodes['doublets']=pd.DataFrame(scrub.predicted_doublets_)
    barcodes['score']=pd.DataFrame(doublet_scores)
    barcodes.to_csv(args.output_dir + 'scrublet_doublet.csv')


    scrub.plot_histogram();
    plt.savefig(args.output_dir + 'scrublet_histogram.pdf')
    print('Running UMAP...')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    print('Done.')
    scrub.plot_embedding('UMAP', order_points=True);
    plt.savefig(args.output_dir + 'scrublet_umap.pdf')

if __name__ == "__main__":
    __main__()
