#! /usr/bin/env python3
# -*- coding: utf-8

import os, argparse, pysam, traceback
import numpy as np
import toolsbox
from time import sleep

from sklearn.cluster import AffinityPropagation as AP
from scipy.spatial.distance import squareform, pdist

from sklearn.metrics import silhouette_score
from sklearn.metrics import calinski_harabasz_score as ch_score
from sklearn.metrics import davies_bouldin_score as db_score


def main():
    desc = "desc here"
    parser = argparse.ArgumentParser(
        prog="Affinity_propagation",
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
        default="results/AF/default",
        help="path to output folder"
    )
    parser.add_argument(
        "--bam",
        metavar='file.bam',
        type=str,
        # required=True
    )
    parser.add_argument(
        "-k",
        "--kmer",
        default=4,
        type=int
    )
    parser.add_argument(
        '-m',
        "--method",
        metavar='[ball-hall, dunn, silhouette, ch-index, db-index]',
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
        "-f",
        "--filter",
        metavar="1500",
        type=int,
        default=0
    )
    parser.add_argument(
        "--cd",
        type=str,
        help="Input the .gff protfile"
    )
    parser.add_argument(
        "-w",
        "--weighting",
        action="store_true",
        default=False,
        help="Weights the tetranucleotide frequency of contigs by their reads coverage."
    )
    args = parser.parse_args()
    if args.bam:
        bam = pysam.AlignmentFile(args.bam, "rb")
    else:
        bam = None
    if args.method not in ["ball-hall", "dunn", "silhouette", "ch-index", "db-index"]:
        print(f"{args.method} not available")
        exit(2)
    kmer = args.kmer
    nf_matrix, c_tables, nucl_list = toolsbox.load(args.input, kmer, bam, nopal=args.nopal, c_filter=args.filter, ponderation=False)
    matrix = squareform(pdist(np.array(toolsbox.recentring(nf_matrix))), "cityblock")
    nf_matrix = np.array(nf_matrix)
    vars_clust = []
    def_cluster = None
    print("Clustering")
    scores, thresholds_computed, clusters = [], [], []
    threshold = [i for i in np.arange(0.50, 1.0, 0.01)]
    for t in threshold:
        cluster = AP(damping=t).fit_predict(matrix)
        dif_clust = []
        for j in list(cluster):
            if j not in dif_clust:
                dif_clust.append(j)
            if -1 in dif_clust:
                dif_clust.remove(-1)
            if args.method in ['ball-hall', "db-index"]:
                score = np.infty
            else:
                score = -np.infty
            if len(dif_clust) > 1:
                if args.method == "ch-index":
                    score = ch_score(matrix, cluster)
                elif args.method == "dunn":
                    score = toolsbox.dunn(cluster, c_tables, matrix)
                elif args.method == "db-index":
                    score = db_score(matrix, cluster)
                elif args.method == "silhouette":
                    score = silhouette_score(matrix, cluster)
                else:
                    score = toolsbox.ball_hall(cluster, matrix)
        scores.append(score)
        thresholds_computed.append(t)
        clusters.append(cluster)
        prefix = " "*3
        if (args.method in ['ball-hall', "db-index"] and score == min(scores)) or (args.method in ["dunn", "silhouette", "ch-index"] and score == max(scores)):
            prefix = "==>"
        print(f"{prefix} {round(t, 3)}: \tnb clust: {len(dif_clust)}\t{args.method}: {round(score,1)}")
        if args.method in ['ball-hall', "db-index"]:
            def_cluster = clusters[scores.index(min(scores))]
        else:
            def_cluster = clusters[scores.index(max(scores))]
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
