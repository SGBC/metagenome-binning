#! /usr/bin/env python3
# -*- coding: utf-8

import os, argparse, pysam
import numpy as np
import toolsbox
from time import sleep

from sklearn.mixture import GaussianMixture as GMM
from sklearn.mixture import BayesianGaussianMixture as VBGMM
from scipy.spatial.distance import squareform, pdist

# from sklearn.metrics import silhouette_score
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
        default="results/VBGMM/default",
        help="path to output folder"
    )
    parser.add_argument(
        "--bam",
        metavar='file.bam',
        type=str,
        # required=True
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
    parser.add_argument(
        "-k",
        "--kmer",
        default=4,
        type=int
    )
    parser.add_argument(
        "-c",
        "--clusters",
        metavar="10",
        type=int
    )
    parser.add_argument(
        "--aic",
        action="store_true",
        default=False,
        help="Use aic measurement instead of bic."
    )
    parser.add_argument(
        "--dbscore",
        action="store_true",
        default=False,
        help="Use Davies Bouldin score measurement instead of Calinski Harabasz score."
    )
    parser.add_argument(
        "--cd",
        type=str,
        help="Input the .gff protfile"
    )
    args = parser.parse_args()
    if args.weighting and not args.bam:
        raise ValueError("--weighting is not available without --bam")
        exit(2)
    if args.bam:
        bam = pysam.AlignmentFile(args.bam, "rb")
    else:
        bam = None
    kmer = args.kmer
    nf_matrix, c_tables, nucl_list = toolsbox.load(args.input, kmer, bam, nopal=args.nopal, c_filter=args.filter, ponderation=args.weighting, cd=args.cd)
    matrix = squareform(pdist(np.array(toolsbox.recentring(nf_matrix))), "cityblock")
    nf_matrix = np.array(nf_matrix)
    print("Clustering")
    lowest_bic, best_nb = np.infty, 0
    if not args.clusters:
        print("Not args.clusters")
        bic_scores = []
        max_ch = -(np.infty)
        n_components_range = range(1, 51)
        cv_types = ['spherical', 'tied', 'diag'] #, 'full']
        for cv in cv_types:
            print(cv)
            for nc in n_components_range:
                gmm = GMM(n_components=nc, covariance_type=cv)
                gmm.fit(matrix)
                score_unit = ''
                if args.aic:
                    score_unit = 'AIC'
                    bic_scores.append(gmm.aic(matrix))
                else:
                    score_unit = 'BIC'
                    bic_scores.append(gmm.bic(matrix))
                if bic_scores[-1] < lowest_bic:
                    best_nb = nc
                    lowest_bic = bic_scores[-1]
                    print(f"==> NC: {nc}, {score_unit}: {bic_scores[-1]}")
                else:
                    print(f"    NC: {nc}, {score_unit}: {bic_scores[-1]}")
    else:
        best_nb = args.clusters
    if args.dbscore:
        max_ch = np.infty
    else:
        max_ch = 0
    weight_opti = 0
    for a in np.arange(0.001, 1.001, 0.01):
        def_cluster = VBGMM(best_nb, covariance_type='tied', weight_concentration_prior=a, weight_concentration_prior_type="dirichlet_distribution").fit_predict(matrix)
        bins = []
        score_unit = ''
        if args.dbscore:
            score_unit = "DB_score"
        else:
            score_unit = "CH_score"
        for i in def_cluster:
            if i not in bins:
                bins.append(i)
        if args.dbscore:
            if db_score(matrix, def_cluster) < max_ch:
                max_ch = db_score(matrix, def_cluster)
                weight_opti = a
            else:
                print("    ", end="")
        else:
            if ch_score(matrix, def_cluster) > max_ch:
                print("==> ", end="")
                max_ch = ch_score(matrix, def_cluster)
                weight_opti = a
            else:
                print("    ", end="")
        print(f"Number of cluster found : {len(bins)}, {round(a, 3)}, {score_unit} {ch_score(matrix, def_cluster)}")
    def_cluster = VBGMM(best_nb, covariance_type='full', weight_concentration_prior=weight_opti).fit_predict(matrix)
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
