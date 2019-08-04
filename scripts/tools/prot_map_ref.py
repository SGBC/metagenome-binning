#! /usr/bin/env python3.7
# -*-coding:utf8-*-

import subprocess
import os
import shutil
import argparse

REF_path = "samples/chromosomes/"


class SoftwareNotFoundError(Exception):
    """Exception to raise when a software is not in the path
    """

    def __init__(self, software):
        super().__init__(f"{software} not found in PATH")


def software_exists(software_name):
    """check if a command-line utility exists and is in the path
    """
    s = shutil.which(software_name)
    if s is not None:
        return s
    else:
        raise SoftwareNotFoundError(software_name)


def prokka(ref_path, ref_name, outdir):
    prokka = software_exists("prokka")
    args = [prokka, "--outdir", outdir]
    args += ["--force", "--prefix", ref_name, ref_path]
    subprocess.run(args)


def main():
    desc = "desc here"
    parser = argparse.ArgumentParser(
        prog="prot_map_ref",
        description=desc
        )
    parser.add_argument(
        "-i",
        "--input",
        metavar="samples/chromosomes/[A-Z]*.fna",
        type=str,
        required=True,
        help="Input contigs files in proteic fasta format",
        nargs="+"
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="prot_map_folder",
        type=str,
        help="Output_folder",
        default="samples/prot_map"
    )
    args = parser.parse_args()
    print("starting")
    os.makedirs(args.output, exist_ok=True)
    refs_path = args.input
    refs_name = [os.path.basename(i).split(".")[0] for i in refs_path]
    for i in range(0, len(refs_path)):
        print(f"run prokka on {refs_name[i]}")
        prokka(refs_path[i], refs_name[i], args.output)
    print("done")

if __name__ == "__main__":
    main()
