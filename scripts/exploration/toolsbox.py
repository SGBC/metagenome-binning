#! /usr/bin/env python3
# -*- coding: utf-8


from itertools import product, combinations

from Bio import SeqIO
import numpy as np
from scipy.spatial.distance import pdist


def load(files, kmer, bam):
    vars_list = ["".join(i) for i in product("ATCG", repeat=kmer)]
    matrix, contigs_index = [], []
    count = 0
    for file in files:
        f = open(file)
        with f:
            fasta = SeqIO.parse(f, "fasta")
            for record in fasta:
                contig_freq = nucleotides_frequences(record.seq, kmer, vars_list)
                if bam:
                    contig_freq.append(bam.count(record.id))
                contigs_index.append(record.id)
                matrix.append(contig_freq)
                count += 1
                print(f"\rImported contigs : {count}", end="")
    print("\nLoading succesfully")
    if bam:
        vars_list.append("Cov")
    np_matrix = np.array(recentring(matrix))
    return np_matrix, contigs_index, vars_list


def min_max(value, mini, maxi):
    return (value-mini)/(maxi-mini)


def recentring(matrix):
    maxi_matrix = list(matrix[0])
    mini_matrix = list(matrix[0])
    for liste in matrix:
        for i in range(len(liste)):
            mini_matrix[i] = min(liste[i], mini_matrix[i])
            maxi_matrix[i] = max(liste[i], maxi_matrix[i])
    recentred_matrix = []
    for l in matrix:
        r_l = [min_max(c, mini, maxi) for c, mini, maxi in zip(l, mini_matrix, maxi_matrix)]
        recentred_matrix.append(r_l)
    return recentred_matrix


def nucleotides_frequences(seq, kmer, vars_list):
    seq = str(seq).upper()
    var_dict = {}
    length = len(seq)-(kmer-1)
    buffer = str(seq[:(kmer-1)])
    for n in seq[(kmer-1):]:
        buffer += str(n)
        if buffer in var_dict.keys():
            var_dict[buffer] += 1
        else:
            var_dict[buffer] = 1
        buffer = str(buffer[1:])
    matrix = []
    for i in vars_list:
        if i in var_dict.keys():
            matrix.append(var_dict[i]/length)
        else:
            matrix.append(0)
    return matrix


def save(files, clust, contigs, output):
    print("sorting contigs")
    bins = {}
    for file in files:
        f = open(file)
        with f:
            fasta = SeqIO.parse(f, "fasta")
            for record in fasta:
                bnumber = clust[contigs.index(record.id)]
                if bnumber in bins.keys():
                    bins[bnumber].append(record)
                else:
                    bins[bnumber] = [record]
    print("Writing data")
    for k in bins.keys():
        SeqIO.write(bins[k], f"{output}/{k}.fa", "fasta")


# def variance(cluster, c_tables, matrix):
#     bins = {}
#     for clust, c_id in zip(cluster, c_tables):
#         if clust in bins.keys():
#             bins[clust].append(c_id)
#         else:
#             bins[clust] = [c_id]
#     variances = []
#     bins_keys = list(bins.keys())
#     bins_keys.sort()
#     for b in bins_keys:
#         count = 0
#         total = 0
#         for i in bins[b]:
#             for j in bins[b][bins[b].index(i)+1:]:
#                 total += matrix[c_tables.index(i)][c_tables.index(j)]
#                 count += 1
#         if count > 0:
#             mean = total/count
#             square_count = 0
#             for i in bins[b]:
#                 for j in bins[b][bins[b].index(i)+1:]:
#                     square_count += (matrix[c_tables.index(i)][c_tables.index(j)]-mean)**2
#             variance = square_count/count
#         else:
#             variance = 1
#         variances.append(variance)
#     if len(bins_keys) > 0:
#         return sum(variances)/len(variances), len(bins_keys)
#     else:
#         return 10**10, 0

def ball_hall(cluster, matrix):
    clusters = {}
    for i in range(len(cluster)):
        if cluster[i] in clusters.keys():
            clusters[cluster[i]].append(matrix[i])
        else:
            clusters[cluster[i]] = [matrix[i]]
    variances = []
    for clust in clusters.keys():
        variances.append(vector_variance(clusters[clust]))
    return sum(variances)/len(variances)


def vector_variance(x):
    x_mean = [0 for i in range(len(x[0]))]
    for i in x:
        x_mean += i
    x_mean /= len(x)
        # for j in range(len(x_mean)):
        #     x_mean[j] += x[i][j]
    # for i in range(len(x_mean)):
    #     x_mean[i] = x_mean[i]/len(x)
    total = 0
    for i in x:
        total += (pdist([i, x_mean], "cityblock"))**2
    # print(total/len(x))
    return float(total/len(x))


# def variance(x):
#     x_mean = sum(x)/len(x)
#     total = 0
#     for i in x:
#         total += (x_mean-i)**2
#     return total/len(x)


def local_min(vars_list):
    min_max = []
    for i in range(0, len(vars_list)):
        if i == 0:
            min_max.append(vars_list[i])
        elif i == len(vars_list)-1:
            min_max.append(vars_list[i])
        elif vars_list[i][0] < min(vars_list[i-1][0], vars_list[i+1][0]) or vars_list[i][0] > max(vars_list[i-1][0], vars_list[i+1][0]):
            min_max.append(vars_list[i])
    print(min_max)
    return min_max[-1]


def plateau(x, y):
    # vars_list = dist, var, nb
    plateaux = []
    if abs(y[-1] - x[1])/len(x) <= 2:
        return x[0]
    for i in range(len(x)-1):
        if abs(y[i+1] - y[i])/2 <= 2:
            for j in range(i+1, len(x)):
                if (abs(y[j]-y[i][1])/(j-i)) > 2:
                    plateaux += [(i, j)]
                    break
    if len(plateaux) == 1:
        return x[plateaux[0][0]]
    elif len(plateaux) != 0:
        scores = [i[1]-i[0] for i in plateaux]
        return x[plateaux[scores.index(max(scores))][0]]
    else:
        return None


def dunn(cluster, c_tables, matrix):
    max_inside = 0
    min_outside = 10**10
    contigs = {}
    clust_name = []
    for c, contig in zip(cluster, c_tables):
        contigs[contig] = c
        if c not in clust_name:
            clust_name.append(c)
    if len(clust_name) > 1:
        for i, j in combinations(contigs.keys(),2):
            if contigs[i] == contigs[j]:
                max_inside = max(max_inside, matrix[c_tables.index(i)][c_tables.index(j)])
            else:
                min_outside = min(min_outside, matrix[c_tables.index(i)][c_tables.index(j)])
        return min_outside/max_inside
    else:
        return 0
