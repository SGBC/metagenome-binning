#! /usr/bin/env python3.7
# -*-coding:utf8-*-

import os, json

from Bio import SeqIO

# nb_reads = 8000000
nb_reads = 19289384
cov_expected = 250
len_reads = 125


def cov(len_reads, nb_reads, prop, genlen):
    c = ((len_reads*nb_reads*prop)/genlen)
    return c


def reads(len_reads, cov, prop, genlen):
    r = ((cov*genlen)/(len_reads*prop))
    return r

bins_path = "samples/pur_set"
query = os.listdir(bins_path)
query_path = [f"{bins_path}/{i}" for i in query]
ref_index = "samples/chromosomes/index_chromo"
ref_index = json.loads(open(ref_index).read())
prop = {}
with open("scripts/abundance_file") as f:
    lines = f.readlines()
    for l in lines:
        lsplit = l.split()
        if lsplit:
            prop[lsplit[0]] = float(lsplit[1])
orga = {}

for k in ref_index.keys():
    orga[k] = 0
    prop[k] = 0
for q in query_path:
    file = open(q)
    with file:
        bin_record = []
        fasta = SeqIO.parse(file, "fasta")
        for record in fasta:
            sid = record.id
            lenght = len(str(record.seq))
            for k in ref_index.keys():
                if sid in ref_index[k]:
                    orga[k] += lenght
                    prop[k] += prop[sid]
for k in orga.keys():
    c = cov(len_reads, nb_reads, prop[k], orga[k])
    print(f"coverage {k:<30} :\t{round(c,2)}")

print("Number of read to obtain a coverage of 5:")

for k in orga.keys():
    r = reads(len_reads, cov_expected, prop[k], orga[k])
    print(f"nb reads total for cov={cov_expected} of {k:<30} :\t{round(r,2)}")
