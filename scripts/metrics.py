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

REF_FOLDER = "samples/chromosomes/"
REF_path = "samples/chromosomes/all_chromo.fna"
META_path = "samples/metabat/fasta_bins"
CONC_path = "samples/concoct/fasta_bins"
TNF_HCLUST = "samples/4NF_hclust_bins"
PNF_HCLUST = "samples/5NF_hclust_bins"
KM_clust = "samples/kmeans_clust"
PUR_SET = "samples/pur_set"


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


def blastx(db, query):
    diamond = software_exists("diamond")
    # args = [diamond, "blastx", "-d", db, "-q", query, "-f", "6", "qseqid", "positive"]
    args = [diamond, "blastx", "-d", db, "-q", query, "-e", "0.0000000001", "-f", "6", "qseqid", "positive"]
    output = subprocess.run(args, capture_output=True)
    out = output.stdout
    out = out.decode()
    out = out.split("\n")
    out = [l.split("\t") for l in out]
    # for l in out:
    #     out[out.index(l)] = l.split("\t")
    total = 0
    id_values = {}
    for i in out:
        if len(i) > 1:
            total += int(i[1])
            if i[0] not in id_values:
                id_values[i[0]] = int(i[1])
            else:
                id_values[i[0]] += int(i[1])
    species = {}
    c_index_file = open("samples/chromosomes/index_chromo")
    c_index = json.loads(c_index_file.read())
    for key in id_values.keys():
        for c_ind in c_index.keys():
            if key in c_index[c_ind]:
                if c_ind in species.keys():
                    species[c_ind] += (id_values[key]/total)
                else:
                    species[c_ind] = (id_values[key]/total)
    # print(f"{round(total,2)}%")
    return species


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
    bins_sect = []
    for b in bins_name:
        for sect in data[b].keys():
            if sect not in bins_sect:
                bins_sect.append(sect)
    values = []
    for sect in bins_sect:
        y_data = []
        for x in bins_name:
            if sect in data[x].keys():
                y_data.append(data[x][sect]*100)
            else:
                y_data.append(0)
        trace = go.Bar(
            x=bins_name,
            y=y_data,
            name=sect
        )
        values.append(trace)
    layout = go.Layout(
        barmode='stack',
        title=f"{setname} bins"
    )
    fig = go.Figure(data=values, layout=layout)
    plot(fig)


def m_andi(ref_path, bins_path, name):
    andi = software_exists("andi")
    refs = os.listdir(ref_path)
    refs = [f"{ref_path}/{i}" for i in refs]
    query = os.listdir(bins_path)
    query = [f"{bins_path}/{i}" for i in query]
    args = [andi, "-j"]
    args += refs + query
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
    maxi = 0
    for r in ref_res:
        line = []
        for b in bins_res:
            value = result[samples.index(r)][samples.index(b)]
            if value == "nan":
                line.append(-10)
            else:
                maxi = max(maxi, float(value))
                line.append(float(value))
        matrix.append(line)
    outvalue = maxi * 1.001
    for x in range(0,len(matrix)):
        for y in range(0,len(matrix[x])):
            if matrix[x][y] == -10:
                matrix[x][y] = 1
            else:
                matrix[x][y] = matrix[x][y]/outvalue
    trace = Heatmap(z=matrix, x=[f"bin_{i}" for i in bins_res], y=ref_res, colorscale=[[0, "#33BF00"], [maxi/outvalue, "#FFFFFF"], [1,"#8c00bf"]])
    data = [trace]
    layout = go.Layout(
        title=f"andi on {name} set."
    )
    fig = go.Figure(data, layout)
    plot(fig)
    # return(samples, result)


# def finch_rs(test_path, setname):
#     bins = os.listdir(test_path)
#     bins_path = [f"{test_path}/{i}" for i in bins]
#     finch = software_exists("finch")
#     f_folder = "samples/finch_temp"
#     os.makedirs(f_folder, exist_ok=True)
#     #refs = [f"samples/pur_sketch/{i}" for i in os.listdir("samples/pur_sketch")]
#     refs = ["samples/pur_sketch/all.sk"]
#     bins_ratio = {}
#     for b in bins_path:
#         contigs = {}
#         total_len = 0
#         bin_name = bins[bins_path.index(b)]
#         for file in os.listdir(f_folder):
#             os.remove(f"{f_folder}/{file}")
#         fasta = SeqIO.parse(b, "fasta")
#         contig_name = []
#         for record in fasta:
#             # print(record.id)
#             SeqIO.write(record, f"{f_folder}/{record.id}", "fasta")
#             contigs[record.id] = {"len" : len(str(record.seq))}
#             total_len += contigs[record.id]["len"]
#             contig_name.append(record.id)
#         seqs = [f"{f_folder}/{i}" for i in os.listdir(f_folder)]
#         args = [finch, "sketch"] + seqs + ["-N","-o", f"{f_folder}/{bin_name[:-3]}.sk"]
#         subprocess.run(args)
#         args = [finch, "dist"] + refs + [f"{f_folder}/{bin_name[:-3]}.sk"]
#         for contig in contig_name:
#             os.remove(f"{f_folder}/{contig}")
#         f_out = subprocess.run(args, capture_output=True)
#         print(f_out.stdout, f_out.stderr)
#         stdout = f_out.stdout.decode()
#         data = json.loads(stdout)
#         for d in data:
#             # print(d["query"])
#             # print(contigs[d["query"]])
#             if "ref" in contigs[d["query"]].keys():
#                 if d["jaccard"] < contigs[d["query"]]["jaccard"]:
#                     contigs[d["query"]]["ref"] = d["reference"]
#                     contigs[d["query"]]["jaccard"] = d["jaccard"]
#             else:
#                 contigs[d["query"]]["ref"] = d["reference"]
#                 contigs[d["query"]]["jaccard"] = d["jaccard"]
#         species = {}
#         for contig in contigs.keys():
#             # print(contig, contigs[contig])
#             if "ref" in contigs[contig].keys():
#                 if contigs[contig]["ref"] in species.keys():
#                     species[contigs[contig]["ref"]] += contigs[contig]["len"]/total_len
#                 else:
#                     species[contigs[contig]["ref"]] = contigs[contig]["len"]/total_len
#             else:
#                 if "unknown" in species.keys():
#                     species["unknown"] += contigs[contig]["len"]/total_len
#                 else:
#                     species["unknown"] = contigs[contig]["len"]/total_len
#         bins_ratio[bin_name[:-3]] = species
#     return bins_ratio


def tests(test_path, ref_path, setname):
    print(f"Run metrics on {setname} bins")
    bins = os.listdir(test_path)
    bins_path = [f"{test_path}/{i}" for i in bins]
    bin_ratio = {}
    for file in os.listdir("samples/prot_map"):
            os.remove(f"samples/prot_map/{file}")
    for file in os.listdir("samples/diamond_db"):
            os.remove(f"samples/diamond_db/{file}")
    for b in bins_path:
        bin_name = bins[bins_path.index(b)]
        print(f"metrics on {bin_name} into {setname}.")
        prodigal(b, bin_name[:-3], output_dir="samples/prot_map")
        d_db_name = f"samples/diamond_db/{bins[bins_path.index(b)][:-3]}"
        d_db_src = f"samples/prot_map/{bins[bins_path.index(b)][:-3]}.faa"
        diamond_db(d_db_name, d_db_src)
        bin_ratio[bin_name] = blastx(db=f"samples/diamond_db/{bin_name[: -3]}", query=REF_path)
    draw(bin_ratio, setname)
    # bin_ratio = finch_rs(test_path, setname)
    # draw(bin_ratio, f"{setname} with finch-rs")
    m_andi(PUR_SET, test_path, setname)

os.makedirs("metrics", exist_ok=True)
path = [KM_clust, META_path, CONC_path, TNF_HCLUST, PUR_SET]
metaname = ["Kmeans_clust", "metabat", "concoct", "4NF_hclust", "Originals chromosomes"]
for p, m in zip(path, metaname):
    tests(p, REF_path, m)
# tests(path[4], REF_path, metaname[4])
