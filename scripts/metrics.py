#! /usr/bin/env python3.7
# -*-coding:utf8-*-

import subprocess
import shutil
import os
import json
# import hashlib

from plotly.offline import plot
from plotly.graph_objs import Heatmap
import plotly.graph_objs as go
from Bio import SeqIO

REF_path = "samples/chromosomes/"
META_path = "samples/metabat/fasta_bins/"
CONC_path = "samples/concoct/fasta_bins/"


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


def m_andi(ref_path, bins_path):
    andi = software_exists("andi")
    refs = os.listdir(ref_path)
    refs = [f"{ref_path}{i}" for i in refs]
    querry = os.listdir(bins_path)
    querry = [f"{bins_path}{i}" for i in querry]
    args = [andi, "-j"]
    args += refs + querry
    args += ["2> /dev/null"]
    andi_out = subprocess.run(args, capture_output=True)
    andi_out = andi_out.stdout.decode('ascii')
    andi_txt = [i.split(" ") for i in andi_out.split("\n")][1:-1]
    samples = [i[0] for i in andi_txt]
    result = [i[1:] for i in andi_txt]
    result = [[i for i in a if i != ""] for a in result]
    ref_res = samples[:len(refs)]
    bins_res = samples[len(refs):]
    matrix = []
    for r in ref_res:
        line = []
        for b in bins_res:
            value = result[samples.index(r)][samples.index(b)]
            if value == "nan":
                line.append(1)
            else:
                line.append(float(value))
        matrix.append(line)
    trace = Heatmap(z=matrix, x=[f"bin_{i}" for i in bins_res], y=ref_res, colorscale=[[0,"#33BF00"], [1,"#FFFFFF"]])
    data = [trace]
    plot(data)
    return(samples, result)


# code from wiz project !
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


def exploit_prod_data(refs, refs_path, hash_ref):
    for i in range(0, len(refs)):
        print(f"run prodigal on {refs[i]}")
        prodigal(refs_path[i], refs[i], output_dir="metrics")
    bins = {}
    for file in refs:
        print(f"Make hash lib for {file}")
        f = open(f"metrics/{file}.faa")
        fasta = SeqIO.parse(f, "fasta")
        for record in fasta:
            # hsh = hashlib.sha256(str(record.seq).encode()).hexdigest()
            hsh = record.seq.lower()
            if not file[:-4] in bins.keys():
                refs_count = {
                    "Bacillus_subtilis": 0,
                    "Cryptococcus_neoformans": 0,
                    "Enterococcus_faecalis": 0,
                    "Escherichia_coli": 0,
                    "Lactobacillus_fermentum": 0,
                    "Listeria_monocytogenes": 0,
                    "Pseudomonas_aeruginosa": 0,
                    "Saccharomyces_cerevisiae": 0,
                    "Salmonella_enterica": 0,
                    "Staphylococcus_aureus": 0,
                    "unidentified": 0
                }
                bins[file] = refs_count
            if hsh in hash_ref.keys():
                bins[file][hash_ref[hsh]] += 1
            else:
                bins[file]["unidentified"] += 1
    for k in bins.keys():
        total = 0
        for l in bins[k].keys():
            total += bins[k][l]
        for l in bins[k].keys():
            bins[k][l] = (bins[k][l]/total, bins[k][l])
    return bins


def draw_prod_data(bins):
    ref = [
        "Bacillus_subtilis",
        "Cryptococcus_neoformans",
        "Enterococcus_faecalis",
        "Escherichia_coli",
        "Lactobacillus_fermentum",
        "Listeria_monocytogenes",
        "Pseudomonas_aeruginosa",
        "Saccharomyces_cerevisiae",
        "Salmonella_enterica",
        "Staphylococcus_aureus",
        "unidentified"
        ]
    trace = []
    bins_names = list(bins.keys())
    for r in ref:
        data = []
        for s in bins_names:
            data.append(bins[s][r][0])
        trace.append(go.Bar(x=bins_names, y=data, name=r))
    layout = go.Layout(barmode="group")
    fig = go.Figure(data=trace, layout=layout)
    plot(fig)

hashs_ref = json.loads(open("metrics/hashs_ref").read())
os.makedirs("metrics", exist_ok=True)
print("Run metrics on Metabat bins")
refs = os.listdir(META_path)
refs_path = [f"{META_path}{i}" for i in refs]
for i in range(0, len(refs)):
    prodigal(refs_path[i], f"META_{i}", output_dir="metrics")
data = exploit_prod_data(refs, refs_path, hashs_ref)
draw_prod_data(data)
metabat = m_andi(REF_path, META_path)
print("Run metrics on Concoct bins")
refs = os.listdir(CONC_path)
refs_path = [f"{CONC_path}{i}" for i in refs]
for i in range(0, len(refs)):
    prodigal(refs_path[i], f"CONC_{i}", output_dir="metrics")
data = exploit_prod_data(refs, refs_path, hashs_ref)
draw_prod_data(data)
concoct = m_andi(REF_path, CONC_path)
