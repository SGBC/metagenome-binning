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


def nucleotides_frequences(seq, kmer, vars_list, nopal):
    seq = str(seq).upper()
    var_dict = {}
    length = len(seq)-(kmer-1)
    buffer = str(seq[:(kmer-1)])
    for n in seq[(kmer-1):]:
        buffer += str(n)
        if buffer in var_dict.keys():
            var_dict[buffer] += 1
        else:
            var_dict[buffer] = 1
        buffer = str(buffer[1:])
    matrix = []
    if nopal:
        nopal_vars = []
        for i in vars_list:
            if i not in nopal_vars and pal(i) not in nopal_vars:
                nopal_vars.append(i)
        for i in nopal_vars:
            score = 0
            if i in var_dict.keys():
                score += var_dict[i]
            if pal(i) in var_dict.keys():
                score += var_dict[pal(i)]
            # print(score)
            matrix.append(score/length)
    else:
        for i in vars_list:
            if i in var_dict.keys():
                matrix.append(var_dict[i]/length)
            else:
                matrix.append(0)
    return matrix


def load(files, kmer, nucl_list, pal=False, sep=False):
    print("Loading contigs")
    l_matrix = []
    contigs_table, gc = [], []
    count = 0
    for file in files:
        f = open(file)
        with f:
            fasta = SeqIO.parse(f, "fasta")
            for record in fasta:
                l_matrix.append(nucleotides_frequences(record.seq, kmer, nucl_list, pal))
                contigs_table.append(record.id)
                gc.append(GC(record.seq))
                count += 1
                print(f"\rImported contigs : {count}", end="")
        if sep:
            nb_var = len(nucl_list)
            if pal:
                nb_var=nb_var//2
            l_matrix.append([0 for i in range(nb_var)])
            gc.append(0)
            contigs_table.append("separation")
    print("\nLoading succesfully")
    # n_matrix = np.array(l_matrix)
    return l_matrix, contigs_table, gc


# def nf(seq, kmer, nucl_list):  # Nucleotides frequence
#     nf_dict = {}
#     length = len(seq)-(kmer-1)
#     buffer = str(seq[:(kmer-1)])
#     for nucl in seq[(kmer-1):]:
#         buffer += str(nucl).capitalize()
#         if buffer in nf_dict.keys():
#             nf_dict[buffer] += 1
#         else:
#             nf_dict[buffer] = 1
#         buffer = str(buffer[1:])
#     nf_list = []
#     for nucl in nucl_list:
#         if nucl in nf_dict.keys():
#             nf_list.append(nf_dict[nucl]/length)
#         else:
#             nf_list.append(0)
#     return nf_list


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
    for l in matrix:
        r_l = [min_max(c, mini, maxi) for c, mini, maxi in zip(l, mini_matrix, maxi_matrix)]
        recentred_matrix.append(r_l)
    return recentred_matrix


def pal(seq):
    pal_seq = [seq[i] for i in range(len(seq)-1, -1, -1)]
    p = {"A": "T", "T": "A", "G": "C", "C": "G"}
    pal_seq = "".join([p[i] for i in pal_seq])
    return pal_seq


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
    parser.add_argument(
        "--nopal",
        action="store_true",
        default=False,
        help="merge palyndromes"
    )
    args = parser.parse_args()

    files = args.input
    nucl_list = ["".join(i) for i in product("ATCG", repeat=args.kmer)]
    matrix, contig_name, gc = load(files, args.kmer, nucl_list, args.nopal)
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
    # print(matrix)
    matrix = recentring(matrix)
    now = datetime.datetime.now()
    data = [Heatmap(z=matrix, y=contig_name, x=nucl_list+["GC", "COV"], colorscale='Viridis')]
    layout = Layout(title=f"{args.name} - {now.day}/{now.month}/{now.year} | {now.hour}h{now.minute}")
    fig = Figure(data, layout)
    os.makedirs("graphs", exist_ok=True)
    plot(fig, filename=f"graphs/Bin_viewer-{args.name}.html")
    print("done")

if __name__ == "__main__":
    main()
