#! /usr/bin/env python3
# -*- coding: utf-8

# import
import os
import argparse

import toolsbox


from sklearn.cluster import AgglomerativeClustering as AC
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform, pdist
from sklearn.decomposition import PCA, NMF
from multiprocessing import Pool
from numpy import arange
import pysam

from sklearn.metrics import silhouette_score
from sklearn.metrics import calinski_harabasz_score as ch_score
from sklearn.metrics import davies_bouldin_score as db_score

from plotly.offline import plot
import plotly.graph_objs as go


"""
<!> Warning: <!>
This code is an experimental prototype.
It probably contains several bugs and is optimized to melt fragile computers
with many useless lines.
Used this one at your own risk.
Maybe in the future this code will be redone and made simpler.
Kisses, the wicked programmer <3
"""


# loading files
def clustering(nf_matrix, c_tables, matrix, dist):
    cluster = AC(affinity="cityblock", compute_full_tree=True, linkage="average", n_clusters=None, distance_threshold=dist).fit_predict(nf_matrix)
    return (dist, cluster)


# main function
def main():
    desc = "desc here"
    parser = argparse.ArgumentParser(
        prog="hclust",
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
        "-k",
        "--kmer",
        default=4,
        type=int
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="output_folder",
        type=str,
        default="results/hclust",
        help="path to output folder"
    )
    # parser.add_argument(
    #     '--pca',
    #     metavar="0.95",
    #     type=float,
    #     help="% of variance conserved in set"
    # )
    # parser.add_argument(
    #     "--nmf",
    #     metavar="25",
    #     type=int
    # )
    parser.add_argument(
        "-c",
        "--clusters",
        metavar="10",
        type=int
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
        help="merge palyndromes"
    )
    args = parser.parse_args()
    if args.method not in ['ball-hall', 'dunn', "silhouette", "ch-index", "db-index"] and not args.clusters:
            print(f"{args.method} was not reconized :/")
            exit(2)

    kmer = args.kmer
    bam = None
    if args.bam:
        bam = pysam.AlignmentFile(args.bam, "rb")
    nf_matrix, c_tables, nucl_list = toolsbox.load(args.input, kmer, bam, nopal=args.nopal)
    print("Pairwise distances calculation")
    output = f"{args.output}"
    os.makedirs(f"{output}", exist_ok=True)
    # if args.pca:
    #     print("PCA activated")
    #     pca = PCA(n_components=args.pca)
    #     nf_matrix = pca.fit_transform(nf_matrix)
    # if args.nmf:
    #     print("NMF activated")
    #     nmf = NMF(n_components=args.nmf)
    #     nf_matrix = nmf.fit_transform(nf_matrix)
    matrix = squareform(pdist(nf_matrix, "cityblock"))
    vars_clust = {}
    def_cluster = None
    print("Clustering")
    if args.clusters:
        def_cluster = AC(affinity="cityblock", compute_full_tree=True, linkage="average", n_clusters=args.clusters).fit_predict(nf_matrix)
    else:
        if args.method:
            seuils = [2, 0.1, 0.01]
            bound_up, bound_down = 52, 1
            # print(f"Distance threshold : [{round(bound_down,3)};{round(bound_up,3)};{seuils[0]}]")
            clusters = {}
            to_write = []
            sdist = 0
            for s in range(len(seuils)):
                var_clust = []
                print(f"Distance threshold : [{round(bound_down,3)};{round(bound_up,3)};{seuils[s]}]")
                with Pool(processes=(os.cpu_count())-1) as pool:
                    queue = [(nf_matrix, c_tables, matrix, i) for i in arange(bound_down, bound_up, seuils[s]) if i not in clusters.keys()]
                    var_clust = pool.starmap(clustering, queue)
                # print(var_clust)
                var_clust.sort()
                for i in range(len(var_clust)):
                    dist, cluster = var_clust[i]
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
                    var_clust[i] = (var, dist, cluster)
                    vars_clust[dist] = cluster
                x = [i[0] for i in var_clust]
                y = [i[1] for i in var_clust]
                if args.method in ['dunn', "silhouette", "ch-index"]:
                    sdist = y[x.index(max(x))]
                else:  # Ball-hall db-score
                    sdist = y[x.index(min(x))]
                if s == len(seuils)-1:
                    def_cluster = vars_clust[sdist]
                    break
                else:
                    bound_down = sdist
                    bound_up = bound_down+seuils[s]

    toolsbox.save(args.input, def_cluster, c_tables, output)


if __name__ == "__main__":
    main()
