#! /usr/bin/env python3
# -*- coding: utf-8

from multiprocessing import Pool
import json, sys

from Bio import SeqIO
import numpy
from scipy.sparse import coo_matrix
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from plotly.offline import plot
from plotly.figure_factory import create_distplot as distplot
from plotly.figure_factory import create_dendrogram as dendrogram
from plotly.graph_objs import Figure
from plotly.graph_objs import Layout
from plotly.graph_objs import Scatter
from plotly.graph_objs import Heatmap


def nucl_frequence(sequence, kmer, nucl_list):
    nf_dict = {}
    length = len(sequence)-(kmer-1)
    buffer = str(sequence[:(kmer-1)])
    for nucl in sequence[(kmer-1):]:
        buffer += str(nucl).capitalize()
        if buffer in nf_dict.keys():
            nf_dict[buffer] += 1
        else:
            nf_dict[buffer] = 1
        buffer = str(buffer[1:])
    nf_list = []
    for nucl in nucl_list:
        if nucl in nf_dict.keys():
            nf_list.append(nf_dict[nucl]/length)
        else:
            nf_list.append(0)
        # matrice = coo_matrix(nf_list)
    return nf_list


def load(files, nucl_list, kmer, output):
    print("Loading contigs")
    contigs = {}
    c_index = {}
    count = 0
    c_2 = 0
    for file in files:
        f = open(file, 'r')
        with f:
            fasta = SeqIO.parse(f, 'fasta')
            for record in fasta:
                print("\rImported : ", c_2, end="")
                contigs[record.id] = nucl_frequence(record.seq, kmer, nucl_list)
                #c_index[record.id] = f"m-{count}.json"
                c_2 += 1
                # if len(contigs) == 1000:
                #     matrice = json.dumps(contigs)
                #     with open(f"{output}/matrice/m-{count}.json", "w") as mf:
                #         mf.write(matrice)
                #    count += 1
                #    contigs = {}
            # if len(contigs):
            #     matrice = json.dumps(contigs)
            #     with open(f"{output}/matrice/m-{count}.json", "w") as mf:
            #         mf.write(matrice)
            #     count += 1
            #     contigs = {}
    print("\nLoading succesfully")
    return contigs


def dendrogram_tetra(c_index, cpu, output):
    """
    Draws a heatmap representing the proximity of the contigs with respect
    to their tetranucleic composition and sorts them by a hierarchy
    clustering. This graph is generated only if at least 2 contigs are
    submitted. The generated graph is returned in HTML or if it could not
    be generated an information message is returned.
    """

    # logger.debug("drawing dendrogram_tetra graph")

    # if more than one contig
    if len(c_index) > 1:
        distances = distance_calculation(c_index, output, cpu)
        id_cells = [contig.uid for contig in contigs]
        cells_value = []

        for id_row in id_cells:
            col_val = []

            for id_col in id_cells:
                if id_col == id_row:
                    col_val.append(0)

                else:
                    if (id_col, id_row) in distances.keys():
                        col_val.append(distances[(id_col, id_row)])

                    else:
                        col_val.append(distances[(id_row, id_col)])

            cells_value.append(col_val)
        cells_value = numpy.array(cells_value)

        # init fig and create upper dendro
        fig = dendrogram(cells_value, labels=id_cells, orientation='bottom')

        # create the colorscale legend
        for i in range(len(fig['data'])):
            fig['data'][i]['yaxis'] = 'y2'

        # create side dendrogram
        fig_side = dendrogram(cells_value, orientation="right")
        for i in range(len(fig_side['data'])):
            fig_side["data"][i]["xaxis"] = "x2"

        # add side dendrogram to fig
        fig.add_traces(fig_side["data"])

        # create heatmap
        dendro_leaves = fig_side['layout']["yaxis"]["ticktext"]
        dendro_leaves = list(map(int, dendro_leaves))
        data_dist = pdist(cells_value)
        heat_data = squareform(data_dist)
        heat_data = heat_data[dendro_leaves, :]
        heat_data = heat_data[:, dendro_leaves]
        heatmap = [
            Heatmap(
                x=dendro_leaves,
                y=dendro_leaves,
                z=heat_data,
                colorscale=[[0, "#080D0F"], [1, "#7BCBF5"]]
                # colorscale=[[0, "#00F532"], [1, "#A80017"]]
            )
        ]
        heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
        heatmap[0]['y'] = fig_side['layout']['yaxis']['tickvals']

        # add heatmap data to fig
        fig.add_traces(heatmap)

        # edit layout
        fig['layout'].update({
            # 'width': 800,
            # 'height': 800,
            'showlegend': False,
            "hovermode": "closest"
            })

        # edit x axis
        fig['layout']['xaxis'].update({
            "domain": [.15, 1],  # heatmap proportion in axe x
            "mirror": False,
            "showgrid": False,
            "showline": False,
            "zeroline": False,
            "ticks": ""
        })

        # edit x axis2
        fig['layout'].update({'xaxis2': {
            'domain': [0, .15],  # side dendrogram proportion in axe x
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': ""
        }})

        # edit y axis
        fig['layout']["yaxis"].update({
            "domain": [0, .85],  # heatmap proportion in axe y
            "mirror": False,
            "showgrid": False,
            "showline": False,
            "zeroline": False,
            "showticklabels": False,
            "ticks": ""
        })

        # edit y axis2
        fig["layout"].update({'yaxis2': {
            "domain": [.825, .975],  # upper dendrogram proportion in axe x
            "mirror": False,
            "showgrid": False,
            "showline": False,
            "zeroline": False,
            "showticklabels": False,
            "ticks": ""
        }})

        # edit title of the graph
        fig['layout'].update(
            title="proximity between two sequences based on their \
                composition tetra nucleic"
        )

        # developement line to draw directly the graph during the run
        plot(fig)

        # logger.debug("drawing succesful")

        # return the graph in html format with javascript
        # return plot(fig, include_plotlyjs=True, output_type='div')


def manhattan_distance(seq1, seq2):
    """
    A function that calculates the Manhattan distance between two
    dictionaries representative of the proportion of tetranucleotides
    in each sequence.
    https://fr.wikipedia.org/wiki/Distance_de_Manhattan
    """
    total = 0
    val_seq1, val_seq2 = 0, 0
    for i, j in zip(seq1, seq2):
        total += abs(j-i)
    return total


def distance_calculation(c_index, output, nb_process=6):
    """
    Function that prepares the tasks for the pool calculates
    and launches this one
    """
    mult = 500000
    results = []
    task = task_list(c_index, nb_process*mult)
    with open(f"{output}/matrice/manhattan", "w") as man:
        man.write("\n")
    print("starting calcul manhattan")
    for i in range(0, (len(c_index)**2-len(c_index)), nb_process*mult):
        tasks = next(task)
        # print(tasks[0])
        results = []
        with Pool(processes=nb_process) as pool:
            results = pool.map(compute_dist, tasks)
        # print(f"\n{len(results)} {sys.getsizeof(results)}")
        print("\rWritting data", 5*".", 15*" ", end="")
        with open(f"{output}/matrice/manhattan", "a") as man:
            man.write(json.dumps(results))
        results = []
    print("end calcul manhattan")
    exit()
    dict_results = {}
    for result in results:
        contig_id, value = result
        dict_results[contig_id] = value
    return dict_results


def task_list(c_index, nb):
    contigs = list(c_index.keys())
    i = 0
    total = len(contigs)**2 - len(contigs)
    tasks = []
    for contig_i in contigs[:-1]:
        for contig_j in contigs[contigs.index(contig_i)+1:]:
            nf_i, nf_j = {}, {}
            task = (
                    (
                        contig_i,
                        contig_j
                    ),
                    (
                        c_index[contig_i],
                        c_index[contig_j]                        
                    )
                )
            tasks.append(task)
            i += 1
            if len(tasks) == nb:
                print(f"\r{round((i/total)*100,2):>8}% | {i:>8}/{total}", end="")
                out = list(tasks)
                tasks = []
                yield out


def compute_dist(task):
    """
    Function used in the parallelization of Manhattan distance
    calculation tasks between two sequences.
    """
    id_task, values = task
    result = manhattan_distance(values[0], values[1])
    return (id_task, result)
