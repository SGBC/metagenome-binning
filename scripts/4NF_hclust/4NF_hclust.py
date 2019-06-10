#! /usr/bin/env python3
# -*- coding: utf-8

# import
import os
from itertools import product
import argparse

import numpy as np
from Bio import SeqIO
from sklearn.metrics import pairwise_distances
from sklearn.cluster import AgglomerativeClustering as AC
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

from plotly.offline import plot
import plotly.graph_objs as go


# loading files
def load(files, kmer, nucl_list):
    print("Loading contigs")
    l_matrix = []
    contigs_table = []
    count = 0
    for file in files:
        f = open(file)
        with f:
            fasta = SeqIO.parse(f, "fasta")
            for record in fasta:
                l_matrix.append(nf(record.seq, kmer, nucl_list))
                contigs_table.append(record.id)
                count += 1
                print(f"\rImported contigs : {count}", end="")
    print("\nLoading succesfully")
    n_matrix = np.array(l_matrix)
    return n_matrix, contigs_table


def save(files, clust, contigs, output):
    print("sorting contigs")
    bins = {}
    for file in files:
        f = open(file)
        with f:
            fasta = SeqIO.parse(f, "fasta")
            for record in fasta:
                bnumber = clust[contigs.index(record.id)]
                if bnumber in bins.keys():
                    bins[bnumber].append(record)
                else:
                    bins[bnumber] = [record]
    print("Writing data")
    for k in bins.keys():
        SeqIO.write(bins[k], f"{output}/{k}.faa", "fasta")


# 4NF computation
def nf(seq, kmer, nucl_list):  # Nucleotides frequence
    nf_dict = {}
    length = len(seq)-(kmer-1)
    buffer = str(seq[:(kmer-1)])
    for nucl in seq[(kmer-1):]:
        buffer += str(nucl).capitalize()
        if buffer in nf_dict.keys():
            nf_dict[buffer] += 1
        else:
            nf_dict[buffer] = 1
        buffer = str(buffer[1:])
    nf_list = []
    for nucl in nucl_list:
        if nucl in nf_dict.keys():
            nf_list.append(nf_dict[nucl]/length)
        else:
            nf_list.append(0)
    return tuple(nf_list)


# main function
def main():
    desc = "desc here"
    parser = argparse.ArgumentParser(
        prog="4NF_clustering",
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
        default="4NF_hclust",
        help="path to output folder"
    )
    parser.add_argument(
        '-c',
        '--cpus',
        metavar=os.cpu_count()-1,
        type=int,
        default=os.cpu_count()-1,
        help="number of cpus allowed"
    )
    args = parser.parse_args()
    nucl_list = ["".join(i) for i in product("ATCG", repeat=args.kmer)]
    nf_matrix, c_tables = load(args.input, args.kmer, nucl_list)
    print("Clustering")
    seuil, n_clust = [], []
    i = 0.0
    while i <= 2:
        i += 0.01
        print(f"\r {round(i,4)}/2", end="")
        output = f"{args.output}/{str(i)}"
        # os.makedirs(f"{output}", exist_ok=True)
        cluster = AC(affinity="cityblock", compute_full_tree=True, linkage="complete", distance_threshold=i, n_clusters=None).fit_predict(nf_matrix)
        seuil.append(i)
        n_clust.append(max(cluster))
        #save(args.input, cluster, c_tables, output)

    trace = go.Scatter(
        x = seuil,
        y = n_clust,
        mode="lines",
        name='Number of clusters found'
    )

    trace2 = go.Scatter(
        x = seuil,
        y = [10 for i in seuil],
        mode='lines',
        name="Real number of organisms"
    )
    layout = go.Layout(
        title=f"Number of cluster on the distance threshold value (Hclust method, linkage complete, affinity cityblock, compute fulltree = true, k={args.kmer})",
        xaxis=dict(title="distance threshold value"),
        yaxis=dict(title="Number of cluster")
    )
    fig = dict(data=[trace, trace2], layout=layout)
    plot(fig)

if __name__ == "__main__":
    main()
# empty line
