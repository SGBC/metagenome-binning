#! /usr/bin/env python3.7
# -*-coding:utf8-*-

import subprocess
import shutil
import os
import json
import argparse
from hashlib import sha256
import pysam
from itertools import product
from time import sleep
import traceback

from sklearn.metrics import homogeneity_score, completeness_score, fowlkes_mallows_score

from plotly.offline import plot
from plotly.graph_objs import Heatmap, Layout, Figure
import plotly.graph_objs as go
from Bio import SeqIO

from binning_viewer import load, recentring, pal
import datetime


def compile_data(bins_paths, metrics_index):
    print("Compile data")
    bins = {}
    error = 0
    seq_number = 0
    for bin_path in bins_paths:
        contigs = []
        with open(bin_path) as f:
            fasta = SeqIO.parse(f, "fasta")
            for record in fasta:
                seq_number += 1
                if record.id in metrics_index.keys():
                    print(f"\r{record.id}", end=" "*10)
                    if sha256(str(record.seq).encode()).hexdigest() == metrics_index[record.id][2]:
                        contigs.append((record.id, metrics_index[record.id][0], len(str(record.seq)), metrics_index[record.id][1]))
                    else:
                        contigs.append((record.id, "Error hash", len(str(record.seq)), "E"))
                        error += 1
                else:
                    contigs.append((record.id, "No data", len(str(record.seq)), "E"))
                    error += 1
        bins[os.path.basename(bin_path)] = contigs
    print()
    bins_name = list(bins.keys())
    print(f"{seq_number} sequences loaded.")
    print(f"{error} sequences create error.")
    bins_name.sort()
    bins_sect = []
    for b in bins_name:
        for contig in bins[b]:
            sect = contig[1]
            if "." in sect:
                sect = sect.split(".")[0]
            if sect not in bins_sect:
                bins_sect.append(sect)
    bins_compo = {}
    bins_ratio = {}
    trust = {"certain": 0, "doubt": 0, "uncertain": 0}
    b_len = 0
    for b in bins_name:
        bins_compo[b] = {}
        for sect in bins_sect:
                bins_compo[b][sect] = 0
        for contig in bins[b]:
            b_len += contig[2]
            # print(contig[3])
            if contig[3] in ["ADL", "ALDO", "AODL", "ADO"]:
                trust['certain'] += contig[2]
            elif contig[3] in ["AUD", "ADU", "ALDOC", "AODLC", "ADOC"]:
                trust['doubt'] += contig[2]
            else:
                trust['uncertain'] += contig[2]
            bins_compo[b][contig[1]] += contig[2]
        bins_ratio[b] = {}
        for c in bins_compo[b].keys():
            bins_ratio[b][c] = bins_compo[b][c]/b_len
    # print("bincompo\n", bins_compo)
    # print("binratio\n", bins_ratio)
    print(f"trust:\n\tCertain: {round((trust['certain']/b_len)*100,2)}%\n\tDoubt: {round((trust['doubt']/b_len)*100,2)}%\n\tUncertain {round((trust['uncertain']/b_len)*100,2)}%")
    return (bins_compo, bins_ratio, seq_number, error, trust)


def precision_recall(compos):
    to_write = []
    global_compo = {}
    precision = {}
    recall = {}
    for b in compos.keys():
        for orga in compos[b].keys():
            if orga in global_compo.keys():
                global_compo[orga] += compos[b][orga]
            else:
                global_compo[orga] = compos[b][orga]
    for b in compos.keys():
        dominant_score = max(compos[b].values())
        dominant_name = ""
        for k in compos[b].keys():
            if compos[b][k] == dominant_score:
                dominant_name = k
        # print(compos[b])
        precision[b] = dominant_score/sum(compos[b].values())
        recall[b] = dominant_score/global_compo[dominant_name]
        to_write.append(f"Bin {b}\tprecision: {round((precision[b]*100), 2)} %")
        to_write.append(f"Bin {b}\trecall: {round(recall[b]*100, 2)} %")
    precision["Mean"] = sum(precision.values())/len(precision.keys())
    recall["Mean"] = sum(recall.values())/len(recall.keys())
    to_write.append(f"Mean precision: {round(precision['Mean']*100, 2)} %")
    to_write.append(f"Mean recall: {round(recall['Mean']*100, 2)} %")
    # for l in to_write:
        # print(l)
    with open("logs.log", 'a') as logs:
        for i in to_write:
            logs.write(i+"\n")
    return precision, recall


def homogeneity_compleness(metrics_index, bins_paths):
    bins_dict = {}
    for bin_file in bins_paths:
        with open(bin_file) as f:
            fasta = SeqIO.parse(f, "fasta")
            for record in fasta:
                bins_dict[record.id] = bins_paths.index(bin_file)
    orga = []
    bins_clusters = []
    metrics_clusters = []
    for key in bins_dict:
        bins_clusters.append(bins_dict[key])
        if not metrics_index[key][0] in orga:
            orga.append(metrics_index[key][0])
        metrics_clusters.append(orga.index(metrics_index[key][0]))
    homo = homogeneity_score(metrics_clusters, bins_clusters)
    compl = completeness_score(metrics_clusters, bins_clusters)
    fm_score = fowlkes_mallows_score(metrics_clusters, bins_clusters)
    print(f"Homogeneity score : {round(homo*100,2)}%\nCompletness score : {round(compl*100, 2)}%\nFowlkes_mallows score: {round(fm_score*100, 2)}%")
    return homo, compl, fm_score


def draw_stat(bins_paths, seq_nb, nb_err, trust, homo, compl, fm_score, setname, draw):
    trace = go.Table(
        header=dict(values=['Statistics of', setname], fill=dict(color=["#666666", "#666666"]), font=dict(color="#ffffff", size=12)),
        cells=dict(values=[
            [
                "Contigs loaded :",
                "Sequences that cause errors",
                "Contigs identified perfectly",
                "Contigs identified only one time",
                "Contigs undentified",
                "Homogeneity score",
                "Completness score",
                "Fowlkes_mallows score",
                "number of clusters"
            ]+["Paths" for i in range(len(bins_paths))],
            [
                seq_nb, nb_err,
                f"{round((trust['certain']/(trust['certain']+trust['doubt']+trust['uncertain']))*100,1)}% ({trust['certain']})",
                f"{round((trust['doubt']/(trust['certain']+trust['doubt']+trust['uncertain']))*100,1)}% ({trust['doubt']})",
                f"{round((trust['uncertain']/(trust['certain']+trust['doubt']+trust['uncertain']))*100,1)}% ({trust['uncertain']})",
                f"{homo*100}%",
                f"{compl*100}%",
                f"{fm_score*100}%",
                len(bins_paths)]+bins_paths
            ], fill=dict(color=[
                    ["#ffffcc", "#ffffcc", "#66cc00", "#ff6600", "#ff0000", "#ffffcc", "#ffffcc", "#ffffcc", "#cccc99"]+["#cccc99" for i in range(len(bins_paths))],
                    ["#ffffcc", "#ffffcc", "#66cc00", "#ff6600", "#ff0000", "#ffffcc", "#ffffcc", "#ffffcc", "#cccc99"]+["#cccc99" for i in range(len(bins_paths))]
                ]
            )
        )
    )
    fig = [trace]
    if draw:
        plot(fig)
        return None
    else:
        return plot(fig, include_plotlyjs=True, output_type='div')


def tests(bins_paths, metrics_index, setname, draw=False, bmap=False, nopal=False):
    compo, ratio, seq_nb, nb_err, trust = compile_data(bins_paths, metrics_index)
    homo, compl, fm_score = homogeneity_compleness(metrics_index, bins_paths)
    precision, recall = precision_recall(compo)
    results = []
        # print("Generating html")
    print("Generating graphs")
    results.append(draw_stat(bins_paths, seq_nb, nb_err, trust, homo, compl, fm_score, setname, draw))
    results.append(draw_pre_rec(precision, recall, setname, draw))
    # results.append(draw_ratio(ratio, setname, draw))
    results.append(draw_compo(compo, setname, draw))
    if bmap:
        results.append(draw_bins_map(bins_paths, draw, setname, nopal))
    return(results)


def draw_bins_map(bins_paths, draw, setname, nopal=False):
    print("Mapping bins")
    nucl_list = ["".join(i) for i in product("ATCG", repeat=4)]
    matrix, contig_name, gc = load(bins_paths, 4, nucl_list, nopal, True)
    bam = pysam.AlignmentFile("samples/mapping/reads.bam", "rb")
    cov = []
    print("Looking for coverage:")
    for i in range(len(contig_name)):
        r = 0
        if "separation" not in contig_name[i]:
            r = int(bam.count(contig_name[i]))
            if not r:
                r = 0
            print(f"\r{round((i/len(contig_name))*100,2):>6}% {contig_name[i]:<10}: {r}", end=" "*10)
        else:
            r = 0
        cov.append(r)
    print("\r100% SUCCES", " "*15)
    for l, c, g in zip(matrix, cov, gc):
        l.append(c)
        # l.append(g)
    matrix = recentring(matrix)
    now = datetime.datetime.now()
    print("Building bins map")
    if nopal:
        temp_nucl_list = []
        for i in nucl_list:
            if i not in temp_nucl_list and pal(i) not in temp_nucl_list:
                temp_nucl_list.append(i)
        nucl_list = list(temp_nucl_list)
    data = [Heatmap(z=matrix, y=contig_name, x=nucl_list+["COV"], colorscale='Viridis')]
    layout = Layout(title=f"{setname} - {now.day}/{now.month}/{now.year} | {now.hour}h{now.minute}")
    fig = Figure(data, layout)
    if draw:
        plot(fig)
        return None
    else:
        return plot(fig, include_plotlyjs=True, output_type='div')


def draw_compo(compo, setname, draw):
    bins_name = list(compo.keys())
    bins_name.sort()
    bins_sect = []
    for b in bins_name:
        for sect in compo[b].keys():
            if sect not in bins_sect:
                bins_sect.append(sect)
    values = []
    for sect in bins_sect:
        y_data = []
        for x in bins_name:
            if sect in compo[x].keys():
                y_data.append(compo[x][sect]*100)
            else:
                y_data.append(0)
        trace = go.Bar(
            x=bins_name,
            y=y_data,
            name=sect
        )
        values.append(trace)
        layout = go.Layout(
            barmode="stack",
            title=f"Composition (b) of each originals organismes in {setname} bins"
        )
    fig = go.Figure(data=values, layout=layout)
    if draw:
        plot(fig)
        return None
    else:
        return plot(fig, include_plotlyjs=True, output_type='div')


def draw_pre_rec(precision, recall, setname, draw):
    bins_name = list(precision.keys())
    bins_name.sort()
    values = []
    values.append(go.Bar(
        x=bins_name,
        y=[precision[k]*100 for k in bins_name],
        name='Sensitivity',
        # mode="markers"
    ))
    values.append(go.Bar(
        x=bins_name,
        y=[recall[k]*100 for k in bins_name],
        name='Specificity',
        # mode="markers"
    ))
    layout = go.Layout(
        barmode='group',
        yaxis=dict(range=[0, 100]),
        title=f"Sensitivity/specificity on {setname} bins"
    )
    fig = go.Figure(data=values, layout=layout)
    if draw:
        plot(fig)
        return None
    else:
        return plot(fig, include_plotlyjs=True, output_type='div')


def main():
    desc = "desc here"
    parser = argparse.ArgumentParser(
        prog="metrics",
        description=desc
        )
    parser.add_argument(
        "-i",
        "--input",
        metavar=".fna",
        type=str,
        required=True,
        help="Input contigs files in nucleic fasta format",
        nargs="+"
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="output_folder",
        type=str,
        help="path to output folder"
    )
    parser.add_argument(
        "-n",
        "--setname",
        metavar='"metabat set"',
        type=str,
        default="{Insert name here} set",
        required=True
    )
    parser.add_argument(
        "-m",
        "--metrics_index",
        metavar='"seq_index.json"',
        type=str,
        required=True
    )
    parser.add_argument(
        "--map",
        action="store_true",
        default=False,
        help="Create a heatmap of bins"
    )
    parser.add_argument(
        "--nopal",
        action="store_true",
        default=False,
        help="merge palyndromes"
    )
    
    args = parser.parse_args()

    metrics_index = json.loads(open(args.metrics_index).read())

    # print(ref_index)
    to_write = []
    draw = True
    if args.output:
        draw = False
    to_write = tests(args.input, metrics_index, args.setname, draw, args.map, args.nopal)

    if args.output:
        os.makedirs(args.output, exist_ok=True)
        with open(f"{args.output}/{args.setname}-report.html", "w") as save:
            save.write("<html><body>")
            for graph in to_write:
                save.writelines(graph)
                save.write("</br></hr>")
            save.write("</body></html>")


if __name__ == "__main__":
    try:
        main()
    except Exception as E:
        for i in range(3):
            print("\a", end='\r')
            sleep(0.33)
        print(traceback.format_exc())
    except KeyboardInterrupt:
        print("\nYou have kill me :'(  MUURRRRDDDDEEEERRRRR !!!!!\a")
    else:
        print("Done")
