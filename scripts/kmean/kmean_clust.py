#! /usr/bin/env python3
# -*- coding: utf-8

# import
import os
from itertools import product
import argparse

import numpy as np
from Bio import SeqIO
from sklearn.cluster import KMeans as KM


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
        default="KM_clust",
        help="path to output folder"
    )
    parser.add_argument(
        '-c',
        '--clusters',
        metavar=10,
        type=int,
        default=10,
        help="number of clusters"
    )
    args = parser.parse_args()
    kmer = args.kmer
    nucl_list = ["".join(i) for i in product("ATCG", repeat=kmer)]
    nf_matrix, c_tables = load(args.input, kmer, nucl_list)
    print("Clustering")
    output = f"{args.output}"
    os.makedirs(f"{output}", exist_ok=True)
    cluster = KM(n_clusters=args.clusters, n_jobs=os.cpu_count()).fit_predict(nf_matrix)
    save(args.input, cluster, c_tables, output)


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
        SeqIO.write(bins[k], f"{output}/{k}.fa", "fasta")


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

if __name__ == "__main__":
    main()
