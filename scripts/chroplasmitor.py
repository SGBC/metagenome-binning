#! /usr/bin/python3
# -*-coding:utf8-*-


import argparse
import os

from Bio import SeqIO


def writter(path, records):
    with open(path, "w") as file:
        for record in records:
            file.write(f">{record.description}\n")
            for i in range(0, len(record.seq), 80):
                file.write(f"{record.seq[i:i+80]}\n")
            file.writelines("\n")


def extract(args):
    for genome in args.genomes:
        f = open(genome, "r")
        plasmid, chromo, mito = [], [], []
        with f:
            fasta_file = SeqIO.parse(f, "fasta")
            for record in fasta_file:
                if "plasmid" in record.description:
                    plasmid.append(record)
                elif "mitochondrion" in record.description:
                    mito.append(record)
                else:
                    chromo.append(record)
        filename = os.path.basename(genome)
        if chromo:
            if not os.path.isdir(f"{args.output}/chromosomes"):
                os.makedirs(f"{args.output}/chromosomes")
            chromo_path = f"{args.output}/chromosomes/{filename}"
            writter(chromo_path, chromo)
        if plasmid:
            if not os.path.isdir(f"{args.output}/plasmids"):
                os.makedirs(f"{args.output}/plasmids")
            plasmid_path = f"{args.output}/plasmids/{filename}"
            writter(plasmid_path, plasmid)
        if mito:
            if not os.path.isdir(f"{args.output}/mitochondrions"):
                os.makedirs(f"{args.output}/mitochondrions")
            plasmid_path = f"{args.output}/mitochondrions/{filename}"
            writter(plasmid_path, mito)


def main():
    desc = "CHROmosomes PLASmids MITochondrion extractOR - \
Extracts and separates chromosome, plamids and mitochondrion from an NCBI file"
    parser = argparse.ArgumentParser(
        prog="chroplasor",
        description=desc
        )
    parser.add_argument(
        "-g",
        "--genomes",
        metavar=".fna",
        type=str,
        required=True,
        help="Input genomes files in nucleic fasta format",
        nargs="+"
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="output_folder",
        type=str,
        required=True,
        help="path to output folder"
    )

    args = parser.parse_args()
    extract(args)


if __name__ == "__main__":
    main()
