#! /usr/bin/env python3
# -*-coding:utf8-*-

import subprocess
import os
import shutil
import argparse
import os


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


def diamond_db(db_name, ref):
    diamond = software_exists("diamond")
    args = [diamond, "makedb", "--in", ref, "-d", db_name]
    subprocess.run(args)


def main():
    desc = "desc here"
    parser = argparse.ArgumentParser(
        prog="diamond_db",
        description=desc
        )
    parser.add_argument(
        "-i",
        "--input",
        metavar=".faa",
        type=str,
        required=True,
        help="Input contigs files in proteic fasta format",
        nargs="+"
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="diamond_db_folder",
        type=str,
        help="Output_folder",
        default="samples/diamond_db"
    )
    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    maps = args.input
    maps_name = [os.path.basename(i) for i in maps]
    dbs = [f"{args.output}/{i.split('.')[0]}" for i in maps_name]
    for i in range(0, len(maps_name)):
        diamond_db(dbs[i], maps[i])
    print("done")

if __name__ == "__main__":
    main()
