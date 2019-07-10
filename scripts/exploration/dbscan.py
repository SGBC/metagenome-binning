#! /usr/bin/env python3
# -*- coding: utf-8

import os, argparse, pysam
import numpy as np
import toolsbox

from sklearn.cluster import DBSCAN
from scipy.spatial.distance import squareform, pdist

from sklearn.metrics import silhouette_score
from sklearn.metrics import calinski_harabasz_score as ch_score
from sklearn.metrics import davies_bouldin_score as db_score


def main():
    desc = "desc here"
    parser = argparse.ArgumentParser(
        prog="DB-scan",
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
        default="results/dbscan/default",
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
    args = parser.parse_args()
    bam = pysam.AlignmentFile(args.bam, "rb")
    kmer = 4
    nf_matrix, c_tables, nucl_list = toolsbox.load(args.input, kmer, bam, nopal=args.nopal, c_filter=args.filter, ponderation=False)
    matrix = squareform(pdist(np.array(toolsbox.recentring(nf_matrix))), "cityblock")
    nf_matrix = np.array(nf_matrix)
    vars_clust = []
    def_cluster = None
    print("Clustering")
    ch_scores, thresholds, clusters = [], [], []
    bound_down, bound_up = 0.1, 15.1
    seuils = [0.5, 0.1]  #, 0.025]
    for s in seuils:
        for i in np.arange(bound_down, bound_up, s):
            if i not in thresholds:
                cluster = DBSCAN(eps=i, min_samples=3).fit_predict(matrix)
                dif_clust = []
                for j in list(cluster):
                    if j not in dif_clust:
                        dif_clust.append(j)
                    if -1 in dif_clust:
                        dif_clust.remove(-1)
                    ch = -42
                    db = -42
                    if len(dif_clust) > 1:
                        ch = ch_score(matrix, cluster)
                        # db = db_score(matrix, cluster)
                ch_scores.append(ch)
                thresholds.append(i)
                clusters.append(cluster)
                print(f"{round(i, 3)}: {list(cluster).count(-1)} \t{len(cluster)}\t{len(dif_clust)}\tch {round(ch,1)}")
        if s == seuils[-1]:
            def_cluster = clusters[ch_scores.index(max(ch_scores))]
            # print(def_cluster, c_tables)
            break
        else:
            bound_down = thresholds[ch_scores.index(max(ch_scores))] - 2*s
            bound_up = thresholds[ch_scores.index(max(ch_scores))] + 2*s
    toolsbox.save(args.input, def_cluster, c_tables, args.output)

if __name__ == "__main__":
    main()