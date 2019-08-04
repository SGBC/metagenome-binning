#! /usr/bin/env python3.7
# -*-coding:utf8-*-

import json
import os

from Bio import SeqIO


print("Run chromo_compiler")
ref = {}

REF_PATH = "samples/chromosomes"
files = os.listdir(REF_PATH)

records = []
for file in files:
    fastafile = open(f"{REF_PATH}/{file}", 'r')
    fasta = SeqIO.parse(fastafile, "fasta")
    for record in fasta:
        records.append(record)
        if file[:-4] in ref:
            ref[file[:-4]].append(record.id)
        else:
            ref[file[:-4]] = [record.id]
SeqIO.write(records, f"{REF_PATH}/all_chromo.fna", "fasta")
file = open(f"{REF_PATH}/index_chromo", "w")
with file:
    file.write(json.dumps(ref))
print("done")
