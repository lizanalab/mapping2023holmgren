import re
import glob
import os
import pandas as pd
import seaborn as sns

sns.set()
import matplotlib.pyplot as plt

from scipy.interpolate import griddata
from sklearn.manifold import MDS
from sklearn.preprocessing import normalize

import numpy as np


def get_gammas(chromosome, alphas):
    gammas = {}

    gamma_ = re.compile(r".+gamma(\d+)_alpha")

    for alpha in alphas:
        files = glob.glob(f"matlab/output/A1_chr{chromosome}_gamma*_alpha{alpha}_partitions.csv")
        gammas_ = [int(gamma_.match(file).group(1)) for file in files]
        gammas[alpha] = sorted(set(gammas_))

    return gammas


def log(p):
    if p < 1e-16:
        return 0.0
    return np.log2(float(p))


def plogp(p):
    return -p * log(p)


def perplexity(a):
    sum_a = sum(a)
    h = sum(plogp(p / sum_a) for p in a)
    return 2 ** h


def effective_size(a):
    sum_a = sum(a)
    return sum_a / perplexity(a)


def load_modularity_scores(chromosome, gammas, alpha=108):
    mod = pd.read_csv(f"matlab/output/A1_chr{chromosome}_alpha{alpha}_modularity.csv", sep=" ", header=None).T
    mod.set_axis(gammas, axis=1, inplace=True)

    return mod


def load_partitions(chromosome, gamma, alpha=108, max_clusters=None):
    partitions = pd.read_csv(f"output/A1_chr{chromosome}_gamma{gamma}_alpha{alpha}_partitions.txt", skiprows=1, sep=" ",
                             dtype=np.int32)

    n_clusters = partitions.ClusterId.max()

    if max_clusters is not None and n_clusters > max_clusters:
        partitions.drop(partitions[partitions.ClusterId > max_clusters].index, inplace=True)
        n_clusters = partitions.ClusterId.max()
        print(f"Including only the best {n_clusters} clusters with {partitions.shape[0]} partitions.")

    return partitions, n_clusters


def partition_validation(filename, threshold=0.45, seed=123):
    from collections import defaultdict
    import os

    name = os.path.splitext(os.path.basename(filename))[0]
    out_name = f"./output/{name}.txt"
    command = f"./partition-validation/partition-validation -s {seed} -t {threshold} {filename} {out_name}"
    os.system(command)

    clusters = defaultdict(list)

    with open(out_name) as f:
        next(f), next(f)
        for line in f:
            line = line.strip()
            cluster, partition = map(int, line.split())
            clusters[cluster - 1].append(partition - 1)

    return dict(clusters)


def run_significance_clustering(agg_file, result_file):
    binary = "./significance-clustering/target/release/significance-clustering"
    cmd = f"{binary} {agg_file} {result_file}"
    os.system(cmd)


def write_partition(result_file, out_file):
    with open(result_file) as f:
        lines = f.readlines()
        lines = lines[1:]  # skip header

    flow = 1.0 / float(len(lines))
    rank = 1

    with open(out_file, "w") as f:
        f.write(f"#{result_file}\n")
        f.write("#path flow name node_id\n")
        for line in lines:
            line = line.strip()
            path, node_id, _frac = line.split(" ")
            node_id = int(node_id)

            insignificant = path.endswith(";")
            sep = "" if insignificant else ":"

            path += sep + str(rank)
            rank += 1

            if ";" in path:
                path += ";"

            f.write(f"{path} {flow} \"{node_id}\" {node_id}\n")


def clamp(value, min_value=0, max_value=1):
    return max(min_value, min(max_value, value))


def plot_solution_landscape(solution_landscape, n_clusters, title=None, out_name=None):
    x = solution_landscape.xcoord.values
    y = solution_landscape.ycoord.values
    modularity = solution_landscape.modularity.values

    def plot_contour(x, y, quality, resolution=1000, contour_method='cubic'):
        resolution = f"{resolution}j"
        X, Y = np.mgrid[min(x):max(x):complex(resolution), min(y):max(y):complex(resolution)]
        points = list(zip(x, y))
        Z = griddata(points, quality, (X, Y), method=contour_method)
        return X, Y, Z

    X, Y, Z = plot_contour(x, y, modularity)

    palette = sns.light_palette("navy", reverse=False, as_cmap=True)
    sns.set(style="whitegrid", font_scale=2)
    f, ax = plt.subplots(figsize=(10, 10))
    sns.despine(f, left=True, bottom=True)
    p0 = ax.contourf(X, Y, Z, cmap=palette, alpha=0.5)

    max_cluster_size = solution_landscape.clustersize.max()

    best = solution_landscape[solution_landscape.index <= n_clusters]

    p1 = sns.scatterplot(x="xcoord", y="ycoord",
                         hue="modularity",
                         size="clustersize",
                         edgecolor=palette(100),
                         sizes=(10, clamp(max_cluster_size, min_value=10, max_value=1000)),
                         alpha=1,
                         palette=palette,
                         legend="full",
                         data=best.iloc[::-1])

    plt.axis('equal')
    if title is not None:
        plt.title(title, pad=30)
    plt.xlabel(None)
    plt.ylabel(None)
    plt.xticks([])
    plt.yticks([])

    f.colorbar(p0, shrink=0.5, anchor=(0.5, 0.4), label="Q")

    handles, labels = ax.get_legend_handles_labels()

    for handle in handles[n_clusters + 1:]:
        handle.set_color(palette(np.max(modularity)))

    labels[n_clusters + 1] = 'Cluster size'

    plt.legend([handles[n_clusters + 1], handles[n_clusters + 2], handles[-1]],
               [labels[n_clusters + 1], labels[n_clusters + 2], labels[-1]],
               bbox_to_anchor=(1.065, 1),
               loc='upper left',
               frameon=False,
               ncol=1,
               labelspacing=1.0)

    if out_name is not None:
        plt.savefig(out_name, bbox_inches="tight")

    return ax
    #
    # kwargs = {
    #     'horizontalalignment': 'left',
    #     'size': 'small',
    #     #'color': palette(0),
    #     'color': 'black',
    #     'weight':'semibold'
    # }
    #
    # for cluster_id in range(n_clusters):
    #     cluster_size, partitionid, modularity, x, y = solution_landscape.iloc[cluster_id]
    #
    #     x_size = cluster_size / max_cluster_size
    #
    #     size_scale = 0.01
    #     x_offset = 0.03
    #     y1_offset = -0.03
    #     y2_offset = 0.007
    #
    #     p1.text(x + size_scale * x_size + x_offset,
    #             y + y1_offset,
    #             np.around(modularity, decimals=3),
    #             **kwargs)
    #
    #     p1.text(x + size_scale * x_size + x_offset,
    #             y + y2_offset,
    #             f"{int(partitionid)} ({int(cluster_size)})",
    #             **kwargs)
    #
    # if out_name is not None:
    #     plt.savefig(out_name, bbox_inches="tight")
