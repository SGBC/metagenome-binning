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

REF_path = "samples/chromosomes"
META_path = "samples/metabat/fasta_bins"
CONC_path = "samples/concoct/fasta_bins"
TNF_HCLUST = "samples/4NF_hclust_bins"
PNF_HCLUST = "samples/5NF_hclust_bins"
KM_clust = "samples/kmeans_clust"


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


def blastx(db, querry):
    diamond = software_exists("diamond")
    args = [diamond, "blastx", "-d", db, "-q", querry, "-f", "6", "qseqid", "positive"]
    # args = [diamond, "blastx", "-d", db, "-q", querry, "-f", "6", "qseqid", "nident"]
    output = subprocess.run(args, capture_output=True)
    out = output.stdout
    out = out.decode()
    out = out.split("\n")
    for l in out:
        out[out.index(l)] = l.split("\t")
    total = 0
    for i in out:
        # print(out)
        if len(i) > 1:
            total+=int(i[1])
    # print(f"{round(total,2)}%")
    return total


def diamond_db(db_name, ref):
    diamond = software_exists("diamond")
    args = [diamond, "makedb", "--in", ref, "-d", db_name]
    subprocess.run(args)


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


def draw(data, setname):
    bins_name = list(data.keys())
    bins_name.sort()
    bins_sect = list(data[bins_name[0]].keys())

    values = []
    for sect in bins_sect:
        trace = go.Bar(
            x=bins_name,
            y=[data[x][sect]*100 for x in bins_name],
            name=sect
        )
        values.append(trace)
    layout = go.Layout(
        barmode='stack',
        title=f"{setname} bins"
    )
    fig = go.Figure(data=values, layout=layout)
    plot(fig)


def tests(test_path, ref_path, metaname):
    print(f"Run metrics on {metaname} bins")
    bins = os.listdir(test_path)
    bins_path = [f"{test_path}/{i}" for i in bins]
    refs_gen = os.listdir(ref_path)
    refs_gen_path = [f"{ref_path}/{i}" for i in refs_gen]
    bin_ratio = {}
    for file in os.listdir("samples/prot_map"):
            os.remove(f"samples/prot_map/{file}")
    for file in os.listdir("samples/diamond_db"):
            os.remove(f"samples/diamond_db/{file}")
    for b in bins_path:
        print(f"metrics on {bins[bins_path.index(b)]} into {metaname}.")
        prodigal(b, bins[bins_path.index(b)][:-3], output_dir="samples/prot_map")
        diamond_db(f"samples/diamond_db/{bins[bins_path.index(b)][:-3]}", f"samples/prot_map/{bins[bins_path.index(b)][:-3]}.faa")
        ref_ratio = {}
        for ref in refs_gen_path:
            ratio = blastx(db=f"samples/diamond_db/{bins[bins_path.index(b)][:-3]}", querry=ref)
            ref_ratio[refs_gen[refs_gen_path.index(ref)]] = ratio
        total = 0
        for r in ref_ratio.keys():
            total += ref_ratio[r]
        if total > 0:
            for r in ref_ratio.keys():
                ref_ratio[r] = ref_ratio[r]/total
        bin_ratio[bins[bins_path.index(b)]] = ref_ratio
    # print(bin_ratio)
    draw(bin_ratio, metaname)
os.makedirs("metrics", exist_ok=True)
path = [KM_clust, META_path, CONC_path, TNF_HCLUST]
metaname = ["Kmeans_clust", "metabat", "concoct", "4NF_hclust"]
# for p, m in zip(path, metaname):
#     tests(p, REF_path, m)
#     # input("step")
tests(KM_clust, REF_path, metaname[0])
