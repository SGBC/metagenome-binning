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
                l_matrix.append(nf(record.seq, 4, nucl_list))
                contigs_table.append(record.id)
                count += 1
                print(f"\rImported contigs : {count}", end="")
    print("\nLoading succesfully")
    # vector_matrix = []
    # coord_dict = {}
    # print("Creation of vector matrix")
    # nb_entry = len(l_matrix)**2
    # cnt_op = 0
    # for i in l_matrix:
    #     line = []
    #     for j in l_matrix:
    #         line.append([i, j])
    #         coord_dict[(len(vector_matrix), len(line)-1)] = (contigs_table[l_matrix.index(i)], contigs_table[l_matrix.index(j)])
    #         cnt_op +=1
    #         print(f"\rVector matrix : {round((cnt_op/nb_entry)*100, 2)}%", end="")
    n_matrix = np.array(l_matrix)
    #contigs_table = np.array(contigs_table)
    return n_matrix, contigs_table


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


# hclust


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
    os.makedirs(f"{args.output}", exist_ok=True)
    nucl_list = ["".join(i) for i in product("ATCG", repeat=args.kmer)]
    nf_matrix, c_tables = load(args.input, args.kmer, nucl_list)
    dist_matrix = pairwise_distances(nf_matrix, metric="cityblock", n_jobs=args.cpus)
    tree = linkage(squareform(dist_matrix), method="complete", metric="cityblock")
    print(tree)


if __name__ == "__main__":
    main()
# empty line
