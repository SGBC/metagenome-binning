#! /usr/bin/env python3
# -*- coding: utf-8

# import
import os
import argparse

import toolsbox


from sklearn.cluster import AgglomerativeClustering as AC
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform, pdist
from sklearn.decomposition import PCA, NMF
from multiprocessing import Pool
from numpy import arange
import pysam

from sklearn.metrics import silhouette_score
from sklearn.metrics import calinski_harabasz_score as ch_score
from sklearn.metrics import davies_bouldin_score as db_score

from plotly.offline import plot
import plotly.graph_objs as go


"""
<!> Warning: <!>
This code is an experimental prototype.
It probably contains several bugs and is optimized to melt fragile computers
with many useless lines.
Used this one at your own risk.
Maybe in the future this code will be redone and made simpler.
Kisses, the wicked programmer <3
"""


# loading files
def clustering(nf_matrix, c_tables, matrix, dist):
    cluster = AC(affinity="cityblock", compute_full_tree=True, linkage="average", n_clusters=None, distance_threshold=dist).fit_predict(nf_matrix)
    var, nb = toolsbox.variance(cluster, c_tables, matrix)
    print(nb, var)
    return (dist, var, nb, cluster)


# main function
def main():
    desc = "desc here"
    parser = argparse.ArgumentParser(
        prog="hclust",
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
        "-k",
        "--kmer",
        default=4,
        type=int
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="output_folder",
        type=str,
        default="results/hclust",
        help="path to output folder"
    )
    parser.add_argument(
        '--pca',
        metavar="0.95",
        type=float,
        help="% of variance conserved in set"
    )
    parser.add_argument(
        "--nmf",
        metavar="25",
        type=int
    )
    parser.add_argument(
        "--bam",
        metavar='file.bam',
        type=str
    )
    args = parser.parse_args()

    kmer = args.kmer
    bam = None
    if args.bam:
        bam = pysam.AlignmentFile(args.bam, "rb")
    nf_matrix, c_tables, nucl_list = toolsbox.load(args.input, kmer, bam)
    print("Pairwise distances calculation")
    output = f"{args.output}"
    os.makedirs(f"{output}", exist_ok=True)
    if args.pca:
        print("PCA activated")
        pca = PCA(n_components=args.pca)
        nf_matrix = pca.fit_transform(nf_matrix)
    if args.nmf:
        print("NMF activated")
        nmf = NMF(n_components=args.nmf)
        nf_matrix = nmf.fit_transform(nf_matrix)
    matrix = squareform(pdist(nf_matrix, "cityblock"))
    vars_clust = []
    def_cluster = None
    print("Clustering")
    seuils = [0.025] #[10, 1, 0.1, 0.01]
    bound_up, bound_down = 52, 1
    print(f"Distance threshold : [{round(bound_down,3)};{round(bound_up,3)};{seuils[0]}]")
    clusters = {}
    to_write = []
    for s in range(len(seuils)):
        with Pool(processes=(os.cpu_count())-1) as pool:
            queue = [(nf_matrix, c_tables, matrix, i) for i in arange(bound_down, bound_up, seuils[s]) if i not in clusters.keys()]
            var_clust = pool.starmap(clustering, queue)
        var_clust.sort()
        for clust in var_clust:
            clusters[clust[0]] = list(clust[3])
            cluster = list(clust[3])
            c = []
            for i in cluster:
                if i not in c:
                    c.append(i)
            var = toolsbox.ball_hall(cluster, matrix)
            dunn = "none"
            sil = "none"
            ch = "none"
            db = "none"
            if len(c) > 1:
                dunn = round(toolsbox.dunn(cluster, c_tables, matrix), 3)
                sil = round(silhouette_score(matrix, cluster, metric='cityblock'), 3)
                ch = round(ch_score(matrix, cluster), 3)
                db = round(db_score(matrix, cluster), 3)
            print(f"{len(c):>3}: Ball-Hall = {round(var,3):<8} Dunn = {dunn:<5} Silhouette coef = {sil:<5} CH_score = {ch:<5} DB_score = {db:<5}")
            to_write.append(f"seuil {clust[0]:<5} cluster {len(c):>5}: Ball-Hall = {round(var,3):<8} Dunn = {dunn:<5} Silhouette coef = {sil:<5} CH_score = {ch:<5} DB_score = {db:<5}")

        variances = [i[1] for i in var_clust]
        mean_variance = sum(variances)/len(variances)
        var_clust = [i[:3] for i in var_clust]
        vars_clust += var_clust
        x = [i[0] for i in var_clust]
        y = [i[1] for i in var_clust]
        if toolsbox.plateau(x, y):
            dist = toolsbox.plateau(x, y)
            dists = [i[0] for i in var_clust]
            print(f"\a<!> {dist} is the distance selected, she produce {var_clust[dists.index(dist)][2]} clusters <!>")
            print("<!> Reason : low slope <!>")
            print("Last clustering")
            def_cluster = clusters[dist]
            break
        elif s == len(seuils)-1:
            var = [i[1] for i in var_clust]
            dist = var_clust[var.index(min(var))][0]
            print(f"\a<!> {dist} is distances selected, she produce {var_clust[var.index(min(var))][2]} clusters <!>")
            print("<!> Reason : maximum depth <!>")
            def_cluster = clusters[dist]
            break
        else:
            var = [i[1] for i in var_clust]
            bound_up = var_clust[var.index(min(var))][0]
            bound_down = bound_up-seuils[s]
            print(f"Distance threshold : [{round(bound_down,3)};{round(bound_up,3)};{seuils[s+1]}], Clusters: [{var_clust[-1][2]};{var_clust[0][2]}]")
    with open("results/hclust_rawdata", "w") as f:
        for i in to_write:
            f.write(str(i+"\n"))

    vars_clust = [(i, j, k) for (j, i, k) in vars_clust]
    toolsbox.save(args.input, def_cluster, c_tables, output)
    vars_clust = [(i, j, k) for (j, i, k) in vars_clust]
    # vars_clust.sort()
    # trace = go.Scatter(
    #     y=[i[1] for i in vars_clust],
    #     x=[i[0] for i in vars_clust],
    #     mode="lines",
    #     name="variance by number of cluster"
    # )

    # layout = go.Layout(
    #     title=f"Mean variance by distance treshold",
    #     xaxis=dict(title="distance treshold"),
    #     yaxis=dict(title="Variance")
    # )
    # fig = dict(data=[trace], layout=layout)
    # plot(fig,  filename='graph0.html')
    # variance_prime = [(vars_clust[i+1][1] - vars_clust[i][1]) for i in range(len(vars_clust)-1)]
    # trace1 = go.Scatter(
    #     y=[i[2] for i in vars_clust],
    #     x=[i[0] for i in vars_clust],
    #     mode="lines",
    #     name="distance threshold by number of cluster"
    # )
    # layout = go.Layout(
    #     title=f"number of clusters by distance treshold",
    #     xaxis=dict(title="distance treshold"),
    #     yaxis=dict(title="number of cluster")
    # )
    # fig = dict(data=[trace1], layout=layout)
    # plot(fig, filename='graph1.html')
    # trace2 = go.Scatter(
    #     y=[(i[0]*i[1]*i[2]) for i in vars_clust],
    #     x=[i[0] for i in vars_clust],
    #     mode="lines",
    #     name="distance threshold by number of cluster"
    # )
    # layout = go.Layout(
    #     title=f"number of clusters by distance treshold",
    #     xaxis=dict(title="distance treshold"),
    #     yaxis=dict(title="number of cluster*variance*distancethreshold")
    # )
    # fig = dict(data=[trace2], layout=layout)
    # plot(fig, filename='graph2.html')


if __name__ == "__main__":
    main()
