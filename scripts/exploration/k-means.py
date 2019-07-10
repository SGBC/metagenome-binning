#! /usr/bin/env python3
# -*- coding: utf-8

# import
import os
import argparse

import toolsbox

import numpy as np
import pysam
from Bio import SeqIO
from sklearn.cluster import KMeans as KM
from scipy.spatial.distance import squareform, pdist
from sklearn.metrics import silhouette_score
from sklearn.metrics import calinski_harabasz_score as ch_score
from sklearn.metrics import davies_bouldin_score as db_score


def main():
    desc = "desc here"
    parser = argparse.ArgumentParser(
        prog="k-means",
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
        default="results/K-means",
        help="path to output folder"
    )
    parser.add_argument(
        "--bam",
        metavar='file.bam',
        type=str,
        required=True
    )
    parser.add_argument(
        '-m',
        "--method",
        metavar='[ball-hall, dunn, silhouette, ch-index, db-index]',
        type=str
    )
    parser.add_argument(
        "-c",
        "--clusters",
        metavar='3',
        type=int,
        help="Number of cluster made by algorithm if any method specified or the maximum of cluster tested if method is specified",
        default=20
    )
    parser.add_argument(
        "--nopal",
        action="store_true",
        default=False,
        help="merge palindromes"
    )
    args = parser.parse_args()
    bam = None
    if args.bam:
        bam = pysam.AlignmentFile(args.bam, "rb")
    kmer = 4
    nf_matrix, c_tables, nucl_list = toolsbox.load(args.input, kmer, bam, nopal=args.nopal)
    matrix = squareform(pdist(np.array(toolsbox.recentring(nf_matrix))), "cityblock")
    nf_matrix = np.array(nf_matrix)
    vars_clust = []
    def_cluster = None
    print("Clustering")
    nb_clusters = args.clusters
    if nb_clusters < 1:
        nb_clusters == 1
    if args.method:
        print("You choice the number of cluster, -m/--methode is not used")
        def_cluster = KM(n_clusters=args.clusters, n_jobs=os.cpu_count()).fit_predict(nf_matrix)
    else:
        if args.method not in ['ball-hall', 'dunn', "silhouette", "ch-index", "db-index"]:
            print(f"{args.method} was not reconized :/")
            exit(2)
        else:
            for i in range(2, nb_clusters+1):
                print(f"Try {i:<2} cluster(s) :", end='')
                cluster = KM(n_clusters=i, n_jobs=os.cpu_count()).fit_predict(matrix)
                var = 0
                if args.method == "ball-hall":
                    var = toolsbox.ball_hall(cluster, matrix)
                elif args.method == "dunn":
                    var = toolsbox.dunn(cluster, c_tables, matrix)
                elif args.method == "silhouette":
                    var = silhouette_score(matrix, cluster, metric="cityblock")
                elif args.method == "ch-index":
                    var = ch_score(matrix, cluster)
                else:  # DB-score
                    var = db_score(matrix, cluster)
                vars_clust.append((var, i, cluster))
                print(var)
            x = [i[0] for i in vars_clust]
            y = [i[1] for i in vars_clust]
            nb_clust = 0
            if args.method in ['dunn', "silhouette", "ch-index"]:
                nb_clust = y[x.index(max(x))]
            else:  # Ball-hall db-score
                nb_clust = y[x.index(min(x))]
            print(f"{nb_clust} clusters is optimal value for {args.method}")
            def_cluster = KM(n_clusters=nb_clust, n_jobs=os.cpu_count()).fit_predict(nf_matrix)
    output = f"{args.output}"
    os.makedirs(f"{output}", exist_ok=True)
    toolsbox.save(args.input, def_cluster, c_tables, output)


if __name__ == "__main__":
    main()
