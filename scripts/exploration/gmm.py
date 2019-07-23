#! /usr/bin/env python3
# -*- coding: utf-8

import os, argparse, pysam
import numpy as np
import toolsbox
from time import sleep
import traceback

from sklearn.mixture import GaussianMixture as GMM
from scipy.spatial.distance import squareform, pdist

from sklearn.metrics import silhouette_score
from sklearn.metrics import calinski_harabasz_score as ch_score
from sklearn.metrics import davies_bouldin_score as db_score


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
        # required=True
    )
    parser.add_argument(
        '-m',
        "--method",
        metavar='[bic, aic, ball-hall, dunn, silhouette, ch-index, db-index]',
        type=str,
        required=True
    )
    parser.add_argument(
        "--nopal",
        action="store_true",
        default=False,
        help="merge palindromes"
    )
    parser.add_argument(
        "-k",
        "--kmer",
        default=4,
        type=int
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
    parser.add_argument(
        "--cd",
        type=str,
        help="Input the .gff protfile",
        metavar="file.gff"
    )
    args = parser.parse_args()
    if args.bam:
        bam = pysam.AlignmentFile(args.bam, "rb")
    else:
        bam = None
    if args.method not in ["bic", "aic", "ball-hall", "dunn", "silhouette", "ch-index", "db-index"]:
        print(f"'{args.method}' is not available.")
        exit(2)
    kmer = args.kmer
    nf_matrix, c_tables, nucl_list = toolsbox.load(args.input, kmer, bam, nopal=args.nopal, c_filter=args.filter, ponderation=args.weighting, cd=args.cd)
    matrix = squareform(pdist(np.array(toolsbox.recentring(nf_matrix))), "cityblock")
    nf_matrix = np.array(nf_matrix)
    print("Clustering")
    scores = []
    best_gmm = None
    if args.method in ["bic", "aic", "ball-hall", "db-index"]:
        best_score = np.infty
    else:
        best_score = -np.infty
    n_components_range = range(1, 25)
    cv_types = ['spherical', 'tied', 'diag', 'full']
    for c in cv_types:
        for nc in n_components_range:
            gmm = GMM(n_components=nc, covariance_type=c)
            gmm.fit(matrix)
            if args.method == "bic":
                scores.append(gmm.bic(matrix))
            elif args.method == "aic":
                scores.append(gmm.aic(matrix))
            clusters = gmm.predict(matrix)
            clust = []
            for i in clusters:
                if i not in clust:
                    clust.append(i)
            if len(clust) > 1:
                if args.method == "ball-hall":
                    scores.append(toolsbox.ball_hall(gmm.predict(matrix), matrix))
                elif args.method == "dunn":
                    scores.append(toolsbox.dunn(gmm.predict(matrix), c_tables, matrix))
                elif args.method == "silhouette":
                    scores.append(silhouette_score(matrix, gmm.predict(matrix)))
                elif args.method == "ch-index":
                    scores.append(ch_score(matrix, gmm.predict(matrix)))
                else:
                    scores.append(db_score(matrix, gmm.predict(matrix)))
            else:
                if args.method in ["ball-hall", "db-index"]:
                    scores.append(np.infty)
                elif args.method in ["dunn", "silhouette", "ch-index"]:
                    scores.append(-np.infty)
            if args.method in ["bic", "aic", "ball-hall", "db-index"]:
                if scores[-1] < best_score:
                    best_gmm = gmm
                    best_score = scores[-1]
                    print("==> ", end='')
                else:
                    print("    ", end='')
            else:
                if scores[-1] > best_score:
                    best_gmm = gmm
                    best_score = scores[-1]
                    print("==> ", end='')
                else:
                    print("    ", end='')
            print(f"Clusters number: {nc}, Cov_type {c}, {args.method}: {scores[-1]}")
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
        print(traceback.format_exc())
    except KeyboardInterrupt:
        print("\nYou have kill me :'(  MUURRRRDDDDEEEERRRRR !!!!!\a")
    else:
        print("Done")
