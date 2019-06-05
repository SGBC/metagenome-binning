#! /usr/bin/env python3.7
# -*-coding:utf8-*-

import subprocess
import shutil
import os

from plotly.offline import plot
from plotly.graph_objs import Heatmap


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

metabat = m_andi(REF_path, META_path)
concoct = m_andi(REF_path, CONC_path)
