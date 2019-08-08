#! /usr/bin/env python3
# -*-coding:utf8-*-

import argparse
import os

import numpy as np
from Bio import SeqIO


def unit(nb_read):
    size = ["", "K", "M", "G", "T", "P", "E", "Z", "Y"]
    counter = 0
    while nb_read >= 1000:
        counter += 1
        nb_read /= 1000
    if counter > 0:
        return f"\t({round(nb_read,3)}{size[counter]})"
    else:
        return ''


def main():
    desc = "description here"
    parser = argparse.ArgumentParser(
        prog="coverage",
        description=desc
    )
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="Input originals files",
        type=str
    )
    parser.add_argument(
        "--reads_len",
        default=125,
        help="Lenght of reads",
        type=int
    )
    parser.add_argument(
        "--abundance",
        help="input the abundance file",
        type=str
    )
    parser_function = parser.add_mutually_exclusive_group()
    parser_function.add_argument(
        "-C",
        "--cov_max",
        help="Input the maximum coverage of contigs",
        type=int
    )
    parser_function.add_argument(
        "-c",
        "--cov_min",
        help="Input the minimum coverage of contigs",
        type=int
    )
    parser_function.add_argument(
        "-r",
        "--nb_reads",
        help="Input the number of reads",
        type=int
    )
    args = parser.parse_args()
    sequence_files = {}
    
    # loading seq ID and total lenght
    for seq in args.input:
        f = open(seq)
        with f:
            records = []
            lenght = 0
            fasta = SeqIO.parse(f, "fasta")
            for record in fasta:
                records.append(record.id)
                lenght += len(record.seq)
            sequence_files[os.path.basename(seq)] = {"lenght": lenght, "records": records}
    
            # setting abundance by seq_files
    if args.abundance:
        abd = {}
        with open(args.abundance) as abundancefile:
            for l in abundancefile.read().split("\n"):
                abd[l.split(" ")[0]] = float(l.split(" ")[-1])
        for seq in sequence_files.keys():
            total_abundance = 0
            for record in sequence_files[seq]["records"]:
                if record in abd.keys():
                    total_abundance += abd[record]
            sequence_files[seq]["abundance"] = total_abundance
    else:
        for seq in sequence_files.keys():
            sequence_files[seq]["abundance"] = 1/len(sequence_files.keys())
    
    # If number of reads (-r) in args
    if args.nb_reads:
        print(f"Reads : {args.nb_reads} {unit(args.nb_reads)}")
        for seq in sequence_files.keys():
            cov = (args.reads_len * args.nb_reads * sequence_files[seq]["abundance"])/sequence_files[seq]["lenght"]
            print(f"coverage {seq:<30} : {round(cov,2)} {unit(cov)}")
    
    # If minimum coverage (-c) in args 
    elif args.cov_min:
        min_abd = np.infty
        seq_name = ""
        for seq in sequence_files.keys():
            if sequence_files[seq]["abundance"] < min_abd:
                min_abd = sequence_files[seq]["abundance"]
                seq_name = seq
        nb_reads = round(((args.cov_min * sequence_files[seq_name]["lenght"])/(args.reads_len))*(1/sequence_files[seq_name]["abundance"]), 0)
        print(f"Reads : {nb_reads} {unit(nb_reads)}")
        for seq in sequence_files.keys():
            cov = (args.reads_len * nb_reads * sequence_files[seq]["abundance"])/sequence_files[seq]["lenght"]
            print(f"coverage {seq:<30} : {round(cov,2)} {unit(cov)}")

    # If maximum coverage (-C) in args
    elif args.cov_max:
        max_abd = -np.infty
        seq_name = ""
        for seq in sequence_files.keys():
            if sequence_files[seq]["abundance"] > max_abd:
                max_abd = sequence_files[seq]["abundance"]
                seq_name = seq
        nb_reads = round((args.cov_max * sequence_files[seq_name]["lenght"])/(args.reads_len * sequence_files[seq_name]["abundance"]),0)
        print(f"Reads : {nb_reads} ({unit(nb_reads)})")
        for seq in sequence_files.keys():
            cov = (args.reads_len * nb_reads * sequence_files[seq]["abundance"])/sequence_files[seq]["lenght"]
            print(f"coverage {seq:<30} : {round(cov,2)} {unit(cov)}")
    
    else:
        print("Any mode specified (-r, -c, -C)")
    print("Bye")

if __name__ == "__main__":
    main()
