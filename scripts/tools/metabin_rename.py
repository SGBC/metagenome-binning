#! /usr/bin/env python3
# -*- coding: utf-8

import os


path = "results/metabat/fasta_bins/"
files = os.listdir(path)
nw_names = []
for file in files:
    nw_names.append(file[4:])
for i, j in zip(files, nw_names):
    if 'bin.' in i:
        if os.path.isfile(f"{path}{j}"):
            os.remove(f"{path}{j}")
        os.rename(f"{path}{i}", f"{path}{j}")
