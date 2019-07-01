#! /usr/bin/env python3.7
# -*-coding:utf8-*-

import subprocess
import shutil
import os
import json
#import hashlib

from Bio import SeqIO

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


def prodigal(genome, output_prefix="bin", trans_table=11, output_dir=None):
    """wrapper function around prodigal
    """

    prodigal = software_exists("prodigal")
    if output_dir:
        out_prot = f"{output_dir}/{output_prefix}.faa"
        out_gff = f"{output_dir}/{output_prefix}.gff"
    else:
        out_prot = f"{output_prefix}.faa"
        out_gff = f"{output_prefix}.gff"
    args = {
        "input": genome,
        "trans_table": str(trans_table),
        "proteins": out_prot,
        "genes": out_gff
    }

    gff = subprocess.run([prodigal, "-i", args["input"], "-g",
                          args["trans_table"], "-m", "-f", "gff",
                          "-a", args["proteins"]],
                         capture_output=True)
    with open(args["genes"], "wb") as f:
        f.write(gff.stdout)

print("starting")
os.makedirs("metrics", exist_ok=True)
refs = os.listdir(REF_path)
refs_path = [f"{REF_path}{i}" for i in refs]
os.makedirs("metrics", exist_ok=True)
for i in range(0, len(refs)):
    print(f"run prodigal on {refs[i][:-4]}")
    prodigal(refs_path[i], refs[i][:-4], output_dir="metrics")
hashs = {}
for file in refs:
    print(f"Make hash lib for {file[:-4]}")
    f = open(f"metrics/{file[:-4]}.faa")
    fasta = SeqIO.parse(f, "fasta")
    for record in fasta:
        # hsh = hashlib.sha256(str(record.seq).encode()).hexdigest()
        hsh = record.seq.lower()
        if hsh in hashs.keys():
            hashs[hsh].append(file[:-4])
        else:
            hashs[hsh] = [file[:-4]]
print("save data")
to_write = json.dumps(hashs)
with open("metrics/hashs_ref", "w") as f:
    f.write(to_write)
