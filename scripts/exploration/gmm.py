#! /usr/bin/env python3
# -*- coding: utf-8

import os, argparse, pysam
import numpy as np
import toolsbox
from time import sleep

from sklearn.mixture import GaussianMixture as GMM
from scipy.spatial.distance import squareform, pdist

# from sklearn.metrics import silhouette_score
from sklearn.metrics import calinski_harabasz_score as ch_score
# from sklearn.metrics import davies_bouldin_score as db_score


def main():
    desc = "desc here"
    parser = argparse.ArgumentParser(
        prog="GMM",
        description=desc
        )
    parser.add_argument(
        "-i",
        "--input",
        metavar=".fna",
        type=str,
        required=True,
        help="Input contigs files in nucleic fasta format",
        nargs="+"
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="output_folder",
        type=str,
        default="results/GMM/default",
        help="path to output folder"
    )
    parser.add_argument(
        "--bam",
        metavar='file.bam',
        type=str,
        required=True
    )
    # parser.add_argument(
    #     '-m',
    #     "--method",
    #     metavar='[ball-hall, dunn, silhouette, ch-index, db-index]',
    #     type=str
    # )
    parser.add_argument(
        "--nopal",
        action="store_true",
        default=False,
        help="merge palindromes"
    )
    parser.add_argument(
        "-f",
        "--filter",
        metavar="1500",
        type=int,
        default=0
    )
    parser.add_argument(
        "-w",
        "--weighting",
        action="store_true",
        default=False,
        help="Weights the tetranucleotide frequency of contigs by their reads coverage."
    )
    args = parser.parse_args()
    bam = pysam.AlignmentFile(args.bam, "rb")
    kmer = 4
    nf_matrix, c_tables, nucl_list = toolsbox.load(args.input, kmer, bam, nopal=args.nopal, c_filter=args.filter, ponderation=args.weighting)
    matrix = squareform(pdist(np.array(toolsbox.recentring(nf_matrix))), "cityblock")
    nf_matrix = np.array(nf_matrix)
    print("Clustering")
    bic_scores= []
    lowest_bic, best_gmm = np.infty , None
    n_components_range = range(1, 25)
    cv_types = ['spherical', 'tied', 'diag', 'full']
    for c in cv_types:
        for nc in n_components_range:
            gmm = GMM(n_components=nc, covariance_type=c)
            gmm.fit(matrix)
            bic_scores.append(gmm.bic(matrix))
            if bic_scores[-1] < lowest_bic:
                best_gmm = gmm
                lowest_bic = bic_scores[-1]
            print(f"NC: {nc}, Cov_type {c}, BIC: {bic_scores[-1]}")
    def_cluster = best_gmm.predict(matrix)
    bins = []
    for i in def_cluster:
        if i not in bins:
            bins.append(i)
    print(f"Number of cluster found : {len(bins)}")
    toolsbox.save(args.input, def_cluster, c_tables, args.output)

if __name__ == "__main__":
    try:
        main()
    except Exception as E:
        for i in range(3):
            print("\a", end='\r')
            sleep(0.33)
        print(f"{E}")
    except KeyboardInterrupt:
        print("\nYou have kill me :'(  MUURRRRDDDDEEEERRRRR !!!!!\a")
    else:
        print("\a Done")
