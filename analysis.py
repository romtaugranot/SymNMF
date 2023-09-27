import math
import sys
import pandas as pd
import numpy as np
import random
import symnmfmodule
from kmeans import k_means_algorithm
from sklearn.metrics import silhouette_score

np.random.seed(0)


def get_cluster_labels(assignment, data_points):
    num_datapoints = len(data_points)

    # Initialize cluster_labels as an array of zeros
    cluster_labels = np.zeros(num_datapoints, dtype=int)

    # Find the cluster of each datapoint
    for i, cluster_points in enumerate(assignment):
        for datapoint in cluster_points:
            index = data_points.index(datapoint)
            cluster_labels[index] = i

    return cluster_labels


def main():
    # Read arguments
    k = int(sys.argv[1])
    file_path = sys.argv[2]

    # Read data points
    data_points = np.genfromtxt(file_path, dtype=float, encoding=None, delimiter=",")
    data_points = data_points.tolist()

    # Compute H
    W = symnmfmodule.norm(data_points)
    n = len(W)
    m = sum(map(sum, W)) / sum(map(len, W))
    max_val = 2 * math.sqrt(m) / math.sqrt(k)
    H = np.random.uniform(low=0, high=max_val, size=(n, k))
    H = H.tolist()
    H = symnmfmodule.symnmf(W, H, n, k)
    H = np.array(H)

    # Compute labels of data points for both algorithms
    symnmf_labels = np.argmax(H, axis=1)
    centroids, kmeans_assignment = k_means_algorithm(k, 300, file_path)
    kmeans_labels = get_cluster_labels(kmeans_assignment, data_points)

    # Compute Silhouette coefficient for both algorithms
    symnmf_score = "{:.4f}".format(silhouette_score(data_points, symnmf_labels))
    kmeans_score = "{:.4f}".format(silhouette_score(data_points, kmeans_labels))

    # Print result
    print("nmf: {}".format(symnmf_score))
    print("kmeans: {}".format(kmeans_score))

if __name__ == "__main__":
    args = sys.argv

    if len(args) != 3:
        print("An Error Has Occurred\n")
        exit()

    try:
        k = int(args[1])
        file_path = args[2]
    except ValueError:
        print("An Error Has Occurred\n")
        exit()

main()
