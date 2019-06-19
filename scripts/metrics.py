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
    # args = [diamond, "blastx", "-d", db, "-q", query,
    #         "-f", "6", "qseqid", "positive"]
    args = [diamond, "blastx", "-d", db, "-q", query,
            "-e", "0.0000000001", "-f", "6", "qseqid", "positive"]
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


def draw(data, setname, force_range=False, type='stack'):
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
    if force_range:
        layout = go.Layout(
            barmode=type,
            title=f"{setname} bins",
            yaxis=dict(range=[0, 100])
        )
    else:
        layout = go.Layout(
            barmode='stack',
            title=f"{setname} bins"
        )
    fig = go.Figure(data=values, layout=layout)
    # plot(fig)
    return plot(fig, include_plotlyjs=True, output_type='div')


def m_andi2(ref_path, bins_path, name):
    andi = software_exists("andi")
    refs = os.listdir(ref_path)
    refs = [f"{ref_path}/{i}" for i in refs]
    query = os.listdir(bins_path)
    query_path = [f"{bins_path}/{i}" for i in query]
    ref_index = "samples/chromosomes/index_chromo"
    ref_index = json.loads(open(ref_index).read())
    query_index = {}
    for q in query_path:
        file = open(q)
        with file:
            bin_record = []
            fasta = SeqIO.parse(file, "fasta")
            for record in fasta:
                sid = record.id
                lenght = len(str(record.seq))
                bin_record.append((sid, lenght))
            query_index[query[query_path.index(q)]] = bin_record
    args = [andi] + refs + query_path + ["2> /dev/null"]
    andi_out = subprocess.run(args, capture_output=True)
    andi_out = andi_out.stdout.decode('ascii')
    andi_txt = [i.split(" ") for i in andi_out.split("\n")][1:-1]
    samples = [i[0] for i in andi_txt]
    result = [i[1:] for i in andi_txt]
    result = [[i for i in a if i != ""] for a in result]
    bins_ratio = {}
    bins_compo = {}
    bins_contigs = {}
    for b in query_index.keys():
        bin_lenght = 0
        bin_contig, bin_compo, bin_ratio = {}, {}, {}
        for key in list(ref_index.keys()) + ["Unknow", "Conflicts"]:
            bin_compo[key] = 0
            bin_contig[key] = 0
        for q in query_index[b]:
            qid, lenght = q
            bin_lenght += lenght
            qid_dist = {}
            for r in ref_index.keys():
                values = []
                for rid in ref_index[r]:
                    value = result[samples.index(qid)][samples.index(rid)]
                    if value == 'nan':
                        values.append(100)
                    else:
                        values.append(float(value))
                qid_dist[r] = value  # sum(values)/len(values)
            min_value = min(qid_dist.values())
            key = []
            if min_value != 100:
                for k in qid_dist.keys():
                    if qid_dist[k] == min_value:
                        key.append(k)
                if len(key) == 1:
                    key = key[0]
                else:
                    key = "Conflicts"
            else:
                key = "Unknow"
            bin_compo[key] += lenght
            bin_contig[key] += 1
        for k in bin_compo.keys():
            bin_ratio[k] = bin_compo[k]/bin_lenght
        bins_ratio[b] = bin_ratio
        bins_compo[b] = bin_compo
        bins_contigs[b] = bin_contig
    for key in list(ref_index.keys()) + ["Unknow", "Conflicts"]:
        bin_ratio[key] = 0
    bins_ratio['total'] = bin_ratio
    total_contig_length = 0
    for b in bins_compo.keys():
        for k in bins_compo[b].keys():
            total_contig_length += bins_compo[b][k]
            bins_ratio["total"][k] += bins_compo[b][k]
    bins_compo["total"] = dict(bins_ratio["total"])
    for key in bins_ratio["total"].keys():
        bins_ratio["total"][key] = bins_ratio["total"][key]/total_contig_length
    graph1 = draw(bins_compo, f"Andi composition of each organisme on {name}")
    graph2 = draw(bins_ratio, f"Andi ratio of each organismes on {name}", True)
    graph3 = draw(precision_recall(bins_ratio, bins_compo),
                  f"Andi precision/recall of {name}", True, 'group')
    with open("graphs.html", "a") as html:
        html.write("<div></br>")
        html.writelines(graph1)
        html.write("</br>")
        html.writelines(graph2)
        html.write("</br>")
        html.writelines(graph3)
        html.write("</br></hr></div>")


def precision_recall(ratio, compos):
    G_tp = 0  # Global_true_positive number
    G_fp = 0  # Gloval_false_positive number
    precisions = []
    recalls = []
    to_write = []
    bins_pre_rec = {}
    for b in ratio.keys():
        if b != "total":
            bins_pre_rec[b] = {}
            B_tp = 0
            B_fp = 0
            max_ratio = max(ratio[b].values())
            key = ""
            for k in ratio[b].keys():
                if ratio[b][k] == max_ratio:
                    key = k
                    break
            for k in compos[b].keys():
                if k == key:
                    B_tp += compos[b][k]
                else:
                    B_fp += compos[b][k]
            # total = compos["total"][key]
            bins_pre_rec[b]["precision"] = ((B_tp)/(B_tp+B_fp))
            bins_pre_rec[b]["recall"] = (B_tp/compos["total"][key])
            to_write.append(f"Bin {b}\tprecision: {round((bins_pre_rec[b]['precision']*100), 2)} %")
            print(to_write[-1])
            to_write.append(f"Bin {b}\trecall: {round(bins_pre_rec[b]['recall']*100, 2)} %")
            print(to_write[-1])
            precisions.append(bins_pre_rec[b]["precision"])
            recalls.append(bins_pre_rec[b]["recall"])
            G_tp += B_tp
            G_fp += B_fp
    # to_write.append(f"Methode precision: {round(
    #   (((G_tp)/(G_tp+G_fp))*100),2)} %")
    bins_pre_rec["Mean"] = {}
    bins_pre_rec["Mean"]["precision"] = sum(precisions)/len(precisions)
    to_write.append(f"Mean precision: {round(bins_pre_rec['Mean']['precision'], 2)} %")
    print(to_write[-1])
    bins_pre_rec["Mean"]["recall"] = sum(recalls)/len(recalls)
    to_write.append(f"Mean recall: {round(bins_pre_rec['Mean']['recall'], 2)} %")
    print(to_write[-1])
    with open("logs.log", 'a') as logs:
        for i in to_write:
            logs.write(i+"\n")
    return bins_pre_rec


def tests(test_path, ref_path, setname):
    print(f"Run metrics on {setname} bins")
    with open("logs.log", 'a') as logs:
        logs.write(f"Run metrics on {setname} bins\n")
    # bins = os.listdir(test_path)
    # bins_path = [f"{test_path}/{i}" for i in bins]
    # bin_ratio = {}
    # for file in os.listdir("samples/prot_map"):
    #         os.remove(f"samples/prot_map/{file}")
    # for file in os.listdir("samples/diamond_db"):
    #         os.remove(f"samples/diamond_db/{file}")
    # for b in bins_path:
    #     bin_name = bins[bins_path.index(b)]
    #     print(f"metrics on {bin_name} into {setname}.")
    #     prodigal(b, bin_name[:-3], output_dir="samples/prot_map")
    #     d_db_name = f"samples/diamond_db/{bins[bins_path.index(b)][:-3]}"
    #     d_db_src = f"samples/prot_map/{bins[bins_path.index(b)][:-3]}.faa"
    #     diamond_db(d_db_name, d_db_src)
    #     bin_ratio[bin_name] = blastx(
    #       db=f"samples/diamond_db/{bin_name[: -3]}",
    #       query=REF_path)
    # draw(bin_ratio, setname)
    # bin_ratio = finch_rs(test_path, setname)
    # draw(bin_ratio, f"{setname} with finch-rs")
    m_andi2(PUR_SET, test_path, setname)


REF_FOLDER = "samples/chromosomes/"
REF_path = "samples/chromosomes/all_chromo.fna"
META_path = "samples/metabat/fasta_bins"
CONC_path = "samples/concoct/fasta_bins"
TNF_HCLUST = "samples/4NF_hclust_bins"
PNF_HCLUST = "samples/5NF_hclust_bins"
KM_clust = "samples/kmeans_clust"
PUR_SET = "samples/pur_set"


os.makedirs("metrics", exist_ok=True)
path = [KM_clust, META_path, CONC_path, TNF_HCLUST, PUR_SET]
metaname = ["Kmeans_clust", "metabat", "concoct",
            "4NF_hclust", "Originals chromosomes"]
for p, m in zip(path, metaname):
    tests(p, REF_path, m)
# tests(path[4], REF_path, metaname[4])
