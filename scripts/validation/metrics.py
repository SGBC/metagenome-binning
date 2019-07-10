#! /usr/bin/env python3.7
# -*-coding:utf8-*-

import subprocess
import shutil
import os
import json
import argparse

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
    args = [diamond, "blastx", "-d", db, "-q", query, '-f', '6', 'qseqid', 'nident']
    output = subprocess.run(args, capture_output=True)
    print(output.stderr)
    out = output.stdout
    out = out.decode()
    out = out.split("\n")
    out = [l.split("\t") for l in out]
    seq_id = {}
    for i in out:
        if len(i) > 1:
            if i[0] in seq_id.keys():
                seq_id[i[0]] += int(i[1])
            else:
                seq_id[i[0]] = int(i[1])
    return seq_id


def diamond_db(db_name, ref):
    diamond = software_exists("diamond")
    args = [diamond, "makedb", "--in", ref, "-d", db_name]
    subprocess.run(args)


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


def andi(ref_paths, bins_paths):
    andi = software_exists("andi")
    if isinstance(ref_paths, list):
        ref_paths = ref_paths
    else:
        ref_paths = [ref_paths]
    args = [andi] + ref_paths + bins_paths + ["2> /dev/null"]
    andi_out = subprocess.run(args, capture_output=True)
    print("End ANDI calculs, parsing")
    if andi_out.stderr:
        print("ANDI ERROR")
    # print(andi_out.stderr)
    andi_out = andi_out.stdout.decode('ascii')
    andi_txt = [i.split(" ") for i in andi_out.split("\n")][1:-1]
    samples = [i[0] for i in andi_txt]
    result = [i[1:] for i in andi_txt]
    result = [[i for i in a if i != ""] for a in result]
    return (samples, result)


def make_index(files_paths):
    query_index = {}
    for f in files_paths:
        with open(f) as file:
            fasta = SeqIO.parse(file, "fasta")
            bin_contigs = {}
            for record in fasta:
                bin_contigs[record.id] = {}
                bin_contigs[record.id]['lenght'] = len(str(record.seq))
            query_index[os.path.basename(f)] = bin_contigs
    return query_index


def andi_data(ref_paths, bins_paths, ref_index):
    print("Running ANDI")
    # ref_index = json.loads(open(ref_index).read())
    # print(ref_paths)
    # print(bins_paths)
    # print(ref_index)
    # input()
    query_index = make_index(bins_paths)
    andi_index, andi_matrix = andi(ref_paths, bins_paths)
    for b in query_index.keys():
        for contig in query_index[b].keys():
            contig_dist = {}
            for gen_name in ref_index.keys():
                values = []
                for ref_contig in ref_index[gen_name]:
                    value = andi_matrix[andi_index.index(contig)][andi_index.index(ref_contig)]
                    if value == "nan":
                        values.append(10000)
                    else:
                        values.append(float(value))
                contig_dist[gen_name] = min(values)
            min_value = min(contig_dist.values())
            keys = []
            if min_value != 10000:
                for key in contig_dist.keys():
                    if contig_dist[key] == min_value:
                        if ".fna" in key:
                            keys.append(key[:-4])
                        else:
                            keys.append(key)
                if len(keys) == 1:
                    query_index[b][contig]["andi"] = keys[0]
                else:
                    query_index[b][contig]["andi"] = keys
            else:
                query_index[b][contig]["andi"] = "Unknow"
    print("End ANDI")
    return query_index


def diamond_data(bins_paths, db_paths, query_index):
    print("Running Diamond")
    for b in query_index.keys():
        for contig in query_index[b].keys():
            query_index[b][contig]["diamond"] = {}
    for bins in bins_paths:
        for db in db_paths:
            genes_score = blastx(db, bins)
            for b in query_index.keys():
                for contig in query_index[b].keys():
                    db_name = os.path.basename(db)
                    # print(db_name.split(".")[0])
                    if contig in genes_score.keys():
                        if "." in db_name:
                            query_index[b][contig]["diamond"][db_name.split(".")[0]] = genes_score[contig]
                        else:
                            query_index[b][contig]["diamond"][db_name] = genes_score[contig]
    for b in query_index.keys():
        for contig in query_index[b].keys():
            # print(query_index[b][contig])
            maxi = 0
            for ref in query_index[b][contig]["diamond"].keys():
                maxi = max(maxi, query_index[b][contig]["diamond"][ref])
            if maxi == 0:
                query_index[b][contig]["diamond"] = "Unknow"
            else:
                refs = []
                for ref in query_index[b][contig]["diamond"].keys():
                    if query_index[b][contig]["diamond"][ref] == maxi:
                        refs.append(ref)
                if len(refs) > 1:
                    query_index[b][contig]["diamond"] = refs
                else:
                    query_index[b][contig]["diamond"] = refs[0]
    # print(query_index[b][contig])
    return query_index


def compile_data(query_index):
    print("Compile data")
    bins_name = list(query_index.keys())
    bins_name.sort()
    bins_sect = []
    for b in bins_name:
        for contig in query_index[b].keys():
            if not isinstance(query_index[b][contig]["andi"], list):
                sect = query_index[b][contig]["andi"]
            if "." in sect:
                sect = sect.split(".")[0]
            if sect not in bins_sect:
                bins_sect.append(sect)
            if not isinstance(query_index[b][contig]["diamond"], list):
                sect = query_index[b][contig]["diamond"]
            if "." in sect:
                sect = sect.split(".")[0]
            if sect not in bins_sect:
                bins_sect.append(sect)
    bins_compo = {}
    bins_ratio = {}
    for b in bins_name:
        b_len = 0
        bins_compo[b] = {}
        bins_compo[b]["andi"] = {}
        bins_compo[b]["diamond"] = {}
        bins_compo[b]["consensus"] = {}
        for sect in bins_sect:
            for i in bins_compo[b].keys():
                bins_compo[b][i][sect] = 0
        for contig in query_index[b].keys():
            # print(b, contig)
            # print(query_index[b])
            b_len += query_index[b][contig]["lenght"]
            bins_compo[b]["andi"][query_index[b][contig]['andi']] += query_index[b][contig]['lenght']
            bins_compo[b]["diamond"][query_index[b][contig]['diamond']] += query_index[b][contig]['lenght']
            if isinstance(query_index[b][contig]['andi'], list) and isinstance(query_index[b][contig]['diamond'], list):
                consensus = []
                for i in query_index[b][contig]['andi']:
                    for j in query_index[b][contig]['diamond']:
                        if i == j:
                            consensus.append(i)
                if len(consensus) == 1:
                    bins_compo[b]["consensus"][consensus[0]] += query_index[b][contig]['lenght']
                else:
                    bins_compo[b]["consensus"]["Conflicts"] += query_index[b][contig]['lenght']
            elif isinstance(query_index[b][contig]['andi'], list) and query_index[b][contig]['diamond'] in query_index[b][contig]['andi']:
                bins_compo[b]["consensus"][query_index[b][contig]['diamond']] += query_index[b][contig]['lenght']
            elif isinstance(query_index[b][contig]['diamond'], list) and query_index[b][contig]['andi'] in query_index[b][contig]['diamond']:
                bins_compo[b]["consensus"][query_index[b][contig]['andi']] += query_index[b][contig]['lenght']
            elif query_index[b][contig]['andi'] == query_index[b][contig]['diamond']:
                bins_compo[b]["consensus"][query_index[b][contig]['andi']] += query_index[b][contig]['lenght']
            elif query_index[b][contig]['andi'] == "Unknow" and query_index[b][contig]['diamond'] != "Unknow":
                bins_compo[b]["consensus"][query_index[b][contig]['diamond']] += query_index[b][contig]['lenght']
            elif query_index[b][contig]['diamond'] == "Unknow" and query_index[b][contig]["andi"] != "Unknow":
                bins_compo[b]["consensus"][query_index[b][contig]['andi']] += query_index[b][contig]['lenght']
            else:
                bins_compo[b]["consensus"]["Conflicts"] += query_index[b][contig]['lenght']
        bins_ratio[b] = {}
        for m in bins_compo[b].keys():
            bins_ratio[b][m] = {}
            for c in bins_compo[b][m].keys():
                bins_ratio[b][m][c] = bins_compo[b][m][c]/b_len
    print("bincompo\n", bins_compo)
    print("binratio\n", bins_ratio)
    return (bins_compo, bins_ratio)


def precision_recall(compos):
    to_write = []
    global_compo = {}
    precision = {}
    recall = {}
    for b in compos.keys():
        for orga in compos[b]["consensus"].keys():
            if orga in global_compo.keys():
                global_compo[orga] += compos[b]["consensus"][orga]
            else:
                global_compo[orga] = compos[b]["consensus"][orga]
    for b in compos.keys():
        dominant_score = max(compos[b]["consensus"].values())
        dominant_name = ""
        for k in compos[b]["consensus"].keys():
            if compos[b]["consensus"][k] == dominant_score:
                dominant_name = k
        print(compos[b])
        precision[b] = dominant_score/sum(compos[b]["consensus"].values())
        recall[b] = dominant_score/global_compo[dominant_name]
        to_write.append(f"Bin {b}\tprecision: {round((precision[b]*100), 2)} %")
        to_write.append(f"Bin {b}\trecall: {round(recall[b]*100, 2)} %")
    precision["Mean"] = sum(precision.values())/len(precision.keys())
    recall["Mean"] = sum(recall.values())/len(recall.keys())
    to_write.append(f"Mean precision: {round(precision['Mean']*100, 2)} %")
    to_write.append(f"Mean recall: {round(recall['Mean']*100, 2)} %")
    for l in to_write:
        print(l)
    with open("logs.log", 'a') as logs:
        for i in to_write:
            logs.write(i+"\n")
    return precision, recall


def tests(bins_paths, ref_paths, db_paths, ref_index, setname, draw=False):
    # bins_paths = [f"{bins_paths}/{i}" for i in os.listdir(bins_paths)]
    #ref_paths = [f"{ref_paths}/{i}" for i in os.listdir(ref_paths)]
    query_index = andi_data(ref_paths, bins_paths, ref_index)
    query_index = diamond_data(bins_paths, db_paths, query_index)
    compo, ratio = compile_data(query_index)
    precision, recall = precision_recall(compo)
    results = [draw_pre_rec(precision, recall, setname, draw)]
    results.append(draw_ratio(ratio, setname, draw))
    results.append(draw_compo(compo, setname, draw))
    return(results)


def draw_compo(compo, setname, draw):
    bins_name = list(compo.keys())
    bins_name.sort()
    bins_sect = []
    for b in bins_name:
        for sect in compo[b]["consensus"].keys():
            if sect not in bins_sect:
                bins_sect.append(sect)
    values = []
    for sect in bins_sect:
        y_data = []
        for x in bins_name:
            if sect in compo[x]['consensus'].keys():
                y_data.append(compo[x]["consensus"][sect]*100)
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


def draw_ratio(ratio, setname, draw):
    bins_name = list(ratio.keys())
    bins_name.sort()
    bins_sect = []
    for b in bins_name:
        for sect in ratio[b]["consensus"].keys():
            if sect not in bins_sect:
                bins_sect.append(sect)
    values = []
    for sect in bins_sect:
        y_data = []
        for x in bins_name:
            # print(ratio[x])
            if sect in ratio[x]["consensus"].keys():
                y_data.append(ratio[x]["consensus"][sect]*100)
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
            title=f"Composition (%) of each originals organismes in {setname} bins",
            yaxis=dict(range=[0, 100])
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
        name='Sensitivity'
    ))
    values.append(go.Bar(
        x=bins_name,
        y=[recall[k]*100 for k in bins_name],
        name='Recall'
    ))
    layout = go.Layout(
        barmode='group',
        yaxis=dict(range=[0, 100]),
        title=f"Sensitivity/recall on {setname} bins"
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
        "-r",
        "--ref",
        default="samples/chromosomes/all_chromo.fna",
        type=str
    )
    parser.add_argument(
        '--db',
        metavar="diamond_db_folder/*",
        type=str,
        help="Diamond db folder with db of ref proteom",
        default=[f"samples/diamond_db/{i}" for i in os.listdir("samples/diamond_db/")],
        nargs="+"
    )
    parser.add_argument(
        "--index",
        metavar="ref sequence name index",
        type=str,
        default="samples/chromosomes/index_chromo"
    )
    parser.add_argument(
        "-n",
        "--setname",
        metavar='"metabat set"',
        type=str,
        default="{Insert name here} set"
    )
    args = parser.parse_args()


    # path = [PUR_SET, META_path, CONC_path, TNF_HCLUST, KM_clust]
    # metaname = ["Originals chromosomes (test)", "metabat", "concoct",
    #             "Aglomerative_clustering", "Kmeans_clustering"]

    ref_index = args.index
    ref_index = json.loads(open(ref_index).read())

    # print(ref_index)
    to_write = []
    draw = True
    if args.output:
        draw = False

    to_write = tests(args.input, args.ref, args.db, ref_index, args.setname, draw)
    # else:
    #     for p, m in zip(path, metaname):
    #         to_write += tests(p, REF_path, DB_Path, ref_index, m)

    if args.output:
        os.makedirs(args.output, exist_ok=True)
        with open(f"{args.output}/results.html", "w") as save:
            save.write("<html><body>")
            for graph in to_write:
                save.writelines(graph)
                save.write("</br></hr>")
            save.write("</body></html>")


# REF_FOLDER = "samples/chromosomes/"
# REF_path = "samples/chromosomes/all_chromo.fna"
# META_path = "samples/metabat/fasta_bins"
# CONC_path = "samples/concoct/fasta_bins"
# TNF_HCLUST = "samples/4NF_hclust_bins"
# PNF_HCLUST = "samples/5NF_hclust_bins"
# KM_clust = "samples/kmeans_clust"
# PUR_SET = "samples/pur_set"
# SMALL_SET = "samples/test"
# DB_Path = [f"samples/diamond_db/{i}" for i in os.listdir("samples/diamond_db/")]

if __name__ == "__main__":
    main()

