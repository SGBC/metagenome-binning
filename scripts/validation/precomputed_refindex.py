#! /usr/bin/env python3.7
# -*-coding:utf8-*-

import subprocess
import shutil
import os
import json
import argparse
from hashlib import sha256
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
    # print(output.stderr)
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
                    if contig in genes_score.keys():
                        if "." in db_name:
                            query_index[b][contig]["diamond"][db_name.split(".")[0]] = genes_score[contig]
                        else:
                            query_index[b][contig]["diamond"][db_name] = genes_score[contig]
    for b in query_index.keys():
        for contig in query_index[b].keys():
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
    return query_index


def compile_data(query_index):
    print("Compile data")
    bins_name = list(query_index.keys())
    bins_name.sort()
    bins_sect = []
    final_index = {}
    for b in bins_name:
        for contig in query_index[b].keys():
            sect = ""
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
    for b in bins_name:
        for contig in query_index[b].keys():
            cid, specie, code = contig,"",""
            if isinstance(query_index[b][contig]['andi'], list) and isinstance(query_index[b][contig]['diamond'], list):
                consensus = []
                for i in query_index[b][contig]['andi']:
                    for j in query_index[b][contig]['diamond']:
                        if i == j:
                            consensus.append(i)
                if len(consensus) == 1:
                    specie = consensus[0]
                    code = "ADL"
                    final_index[cid] = (specie, code)
                else:
                    specie = "Conflict"
                    code = "ADLC"
            elif isinstance(query_index[b][contig]['andi'], list):
                if query_index[b][contig]['diamond'] in query_index[b][contig]['andi']:
                    specie = query_index[b][contig]['diamond']
                    if specie != "Unknow":
                        code = "ALDO"
                    else:
                        code = "ALDUC"
                else:
                    specie = query_index[b][contig]['diamond']
                    # specie = "Conflict"
                    code = "ALDOC"
            elif isinstance(query_index[b][contig]['diamond'], list):
                if query_index[b][contig]['andi'] in query_index[b][contig]['diamond']:
                    specie = query_index[b][contig]['andi']
                    if specie != "Unknow":
                        code = "AODL"
                    else:
                        code = "AUDLC"
                else:
                    specie = query_index[b][contig]['andi']
                    # specie = "Conflict"
                    code = "AODLC"
            elif query_index[b][contig]['andi'] == query_index[b][contig]['diamond']:
                specie = query_index[b][contig]["andi"]
                if specie != "Unknow":
                    code = "ADO"
                else:
                    code = "U"
            elif query_index[b][contig]['andi'] == "Unknow":
                specie = query_index[b][contig]['diamond']
                code = "AUD"
            elif query_index[b][contig]['diamond'] == "Unknow":
                specie = query_index[b][contig]['andi']
                code = "ADU"
            else:
                bins_compo[b]["consensus"]["Conflicts"] += query_index[b][contig]['lenght']
            final_index[cid] = [specie, code]
    return final_index


def hash_calc(final_index, ref_path):
    for file in ref_path:
        with open(file) as fasta_file:
            fasta = SeqIO.parse(fasta_file, "fasta")
            for record in fasta:
                final_index[record.id].append(sha256(str(record.seq).encode()).hexdigest())
    return final_index


def tests(bins_paths, ref_paths, db_paths, ref_index):
    query_index = andi_data(ref_paths, bins_paths, ref_index)
    query_index = diamond_data(bins_paths, db_paths, query_index)
    final_index = compile_data(query_index)
    final_index = hash_calc(final_index, bins_paths)
    for k in final_index.keys():
        print(k, "\t", final_index[k][1], "\t", final_index[k][0])
    return(final_index)


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
        help="path to output folder",
        default="results"
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

    ref_index = args.index
    ref_index = json.loads(open(ref_index).read())

    to_write = ""

    to_write = json.dumps(tests(args.input, args.ref, args.db, ref_index))

    os.makedirs(args.output, exist_ok=True)
    with open(f"{args.output}/seq_index.json", "w") as save:
        save.write(to_write)

if __name__ == "__main__":
    main()
