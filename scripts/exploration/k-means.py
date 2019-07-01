#! /usr/bin/env python3
# -*- coding: utf-8

# import
import os
import argparse

import toolsbox

import numpy as np
from Bio import SeqIO
from sklearn.cluster import KMeans as KM
from scipy.spatial.distance import squareform, pdist


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
        type=str
    )
    args = parser.parse_args()
    bam = None
    if args.bam:
        bam = pysam.AlignmentFile(args.bam, "rb")
    kmer = 4
    nf_matrix, c_tables, nucl_list = toolsbox.load(args.input, kmer, bam)
    matrix = squareform(pdist(nf_matrix, "cityblock"))
    vars_clust = []
    def_cluster = None
    print("Clustering")
    for i in range(1, 25):
        print(f"Try {i:<2} cluster(s) :", end='')
        cluster = KM(n_clusters=i, n_jobs=os.cpu_count()).fit_predict(matrix)
        # print(cluster)
        var = toolsbox.ball_hall(cluster, nf_matrix)
        dunn = toolsbox.dunn(cluster, c_tables, matrix)
        print(f" Ball-Hall = {round(var,3):<8} Dunn = {round(dunn, 3)}")
        vars_clust.append((var, i, cluster))
    x = [i[0] for i in vars_clust]
    y = [i[1] for i in vars_clust]
    if toolsbox.plateau(x, y):
        nb_clust = toolsbox.plateau(x, y)
        def_cluster = vars_clust[x.index(nb_clust)][2]
    output = f"{args.output}"
    os.makedirs(f"{output}", exist_ok=True)
    # cluster = KM(n_clusters=args.clusters, n_jobs=os.cpu_count()).fit_predict(nf_matrix)
    save(args.input, def_cluster, c_tables, output)


if __name__ == "__main__":
    main()
