#! /usr/bin/env python3
# -*- coding: utf-8

import pysam
from Bio import SeqIO
from itertools import product
import numpy as np
from plotly.offline import plot
from plotly.figure_factory import create_distplot as distplot
from Bio.SeqUtils import GC
from plotly.graph_objs import Heatmap, Layout, Figure
import os
import argparse
import datetime

def load(files, kmer, nucl_list):
    print("Loading contigs")
    l_matrix = []
    contigs_table, gc = [], []
    count = 0
    for file in files:
        f = open(file)
        with f:
            fasta = SeqIO.parse(f, "fasta")
            for record in fasta:
                l_matrix.append(nf(record.seq, kmer, nucl_list))
                contigs_table.append(record.id)
                gc.append(GC(record.seq))
                count += 1
                print(f"\rImported contigs : {count}", end="")
    print("\nLoading succesfully")
    # n_matrix = np.array(l_matrix)
    return l_matrix, contigs_table, gc


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
    return nf_list


def min_max(value, mini, maxi):
    return (value-mini)/(maxi-mini)


def recentring(matrix):
    maxi_matrix = list(matrix[0])
    mini_matrix = list(matrix[0])
    # print(mini_matrix, maxi_matrix)
    for liste in matrix:
        for i in range(len(liste)):
            mini_matrix[i] = min(liste[i], mini_matrix[i])
            maxi_matrix[i] = max(liste[i], maxi_matrix[i])
    recentred_matrix = []
    # print(mini_matrix, maxi_matrix)
    for l in matrix:
        r_l = [min_max(c, mini, maxi) for c, mini, maxi in zip(l, mini_matrix, maxi_matrix)]
        recentred_matrix.append(r_l)
    return recentred_matrix


def main():
    desc = "desc here"
    parser = argparse.ArgumentParser(
        prog="stat viewer",
        description=desc
        )
    parser.add_argument(
        "-i",
        "--input",
        metavar=".fna",
        type=str,
        required=True,
        help="Input bins files",
        nargs="+"
    )
    parser.add_argument(
        "-k",
        "--kmer",
        metavar='4',
        type=int,
        default=4
    )
    parser.add_argument(
        "-n",
        "--name",
        metavar="My metabat binning",
        type=str,
        required=True,
        help="The name of the analisys"
    )
    args = parser.parse_args()

    # files = ["samples/contigs/final.contigs.fa"]
    # folder = "samples/metabat/fasta_bins/"
    # folder = args.input
    # if folder[-1] == "/":
    #     folder = folder[:-1]
    # files = [f"{folder}/{i}" for i in os.listdir(folder)]
    files = args.input
    nucl_list = ["".join(i) for i in product("ATCG", repeat=args.kmer)]
    matrix, contig_name, gc = load(files, args.kmer, nucl_list)
    bam = pysam.AlignmentFile("samples/mapping/reads.bam", "rb")
    cov = []
    print("Looking for coverage:")
    for i in range(len(contig_name)):
        r = int(bam.count(contig_name[i]))
        print(f"\r{round((i/len(contig_name))*100,2):>6}% {contig_name[i]:<10}: {r}", end=" "*10)
        cov.append(r)
    print("\r100% SUCCES", " "*15)
    for l, c, g in zip(matrix, cov, gc):
        l.append(g)
        l.append(c)
    matrix = recentring(matrix)
    now = datetime.datetime.now()
    data = [Heatmap(z=matrix, y=contig_name, x=nucl_list+["GC", "COV"], colorscale='Viridis')]
    layout = Layout(title=f"{args.name} - {now.day}/{now.month}/{now.year} | {now.hour}h{now.minute}")
    fig = Figure(data, layout)
    plot(fig, filename=f"Bin_viewer-{args.name}.html")
    print("done")

if __name__ == "__main__":
    main()
