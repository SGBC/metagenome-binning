#! /usr/bin/env python3
# -*-coding:utf8-*-

"""
Chroplasmitor: CHROmosomes PLASmids MITochondrions extractOR
Chroplasmitor is used to sort contigs based on their entry name in a fasta file.
If the word "plasmid" is in the name of the entry, that if will go into the plasmids folder.
If the word "mitochondrion" is in the name of the entry, that will go into the mitochondrial file.
The rest will go into the chromosomes file.
"""
import argparse
import os

from Bio import SeqIO


def writer(path, records):
    with open(path, "w") as file:
        SeqIO.write(records, path, "fasta")
        # for record in records:
            # file.write(f">{record.description}\n")
            # for i in range(0, len(record.seq), 80):
            #     file.write(f"{record.seq[i:i+80]}\n")
            # file.writelines("\n")


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
            writer(chromo_path, chromo)
        if plasmid:
            if not os.path.isdir(f"{args.output}/plasmids"):
                os.makedirs(f"{args.output}/plasmids")
            plasmid_path = f"{args.output}/plasmids/{filename}"
            writer(plasmid_path, plasmid)
        if mito:
            if not os.path.isdir(f"{args.output}/mitochondrions"):
                os.makedirs(f"{args.output}/mitochondrions")
            plasmid_path = f"{args.output}/mitochondrions/{filename}"
            writer(plasmid_path, mito)


def main():
    desc = "CHROmosomes PLASmids MITochondrions extractOR - \
Extracts and separates chromosome, plamids and mitochondrion from a fasta file"
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
