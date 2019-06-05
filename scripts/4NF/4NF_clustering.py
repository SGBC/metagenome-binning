#! /usr/bin/env python3
# -*- coding: utf-8

import argparse

import os
from itertools import product

from F_4NF import *


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
        default="4NF_output",
        help="path to output folder"
    )
    args = parser.parse_args()
    os.makedirs(f"{args.output}/matrice", exist_ok=True)
    nucl_list = ["".join(i) for i in product("ATCG", repeat=args.kmer)]
    contigs = load(args.input, nucl_list, args.kmer, args.output)
    dendrogram_tetra(contigs, cpu=os.cpu_count()-1, output=args.output)
    print(end)


if __name__ == "__main__":
    main()
