import sys
import math
import random
import numpy as np
import pandas as pd
import symnmfmodule
from sklearn.metrics import silhouette_score

np.random.seed(0)

# *** Code for kmeans.py *** #

eps = 0.001


# *** Math functions *** #

# pre-condition: vector1's length == vector2's length
def euclidean_distance(vector1, vector2):
    dis = 0
    for i in range(len(vector1)):
        dis += (vector2[i] - vector1[i]) ** 2
    return math.sqrt(dis)


# pre-condition: vector1's length == vector2's length
def sum_vectors(vector1, vector2):
    vector = [0 for _ in range(len(vector1))]
    for i in range(len(vector1)):
        vector[i] = float("{:.4f}".format(vector1[i] + vector2[i]))
    return vector


def divide_by_scalar(vector, scalar):
    new_vector = [0 for _ in range(len(vector))]
    for i in range(len(vector)):
        new_vector[i] = float("{:.4f}".format(vector[i] / scalar))
    return new_vector


def read_file():
    iter = 200  # default value.
    i = 3  # represents the index of the file in the args array, assuming iter was provided.
    data_points = []
    flag_iter = True
    if len(sys.argv) == 4:  # iter was provided
        try:
            iter = int(sys.argv[2])
        except ValueError:
            flag_iter = False
    else:  # iter was not provided
        i -= 1
    flag_file = True
    try:
        with open(sys.argv[i]) as file:
            for line in file:
                vector = line.rstrip().split(',')
                vector = [eval(i) for i in vector]
                data_points += [vector]
    except OSError:
        flag_file = False
    return flag_file, flag_iter, data_points, iter


def initialize():
    # default values
    K = 0
    flag_K = True
    try:
        K = int(sys.argv[1])
    except ValueError:
        flag_K = False
    flag_file, flag_iter, data_points, iter = read_file()
    return flag_K, flag_iter, flag_file, K, iter, data_points


# pre-condition: old_centroids' length == new_centroids' length
def compute_flag_delta(old_centroids, new_centroids):  # returns True if and only if every delta is less than eps.
    flag_delta = True
    for i in range(len(old_centroids)):
        delta = euclidean_distance(old_centroids[i], new_centroids[i])
        if delta >= eps:
            flag_delta = False
            break
    return flag_delta


def find_min_index(data_point, centroids):
    min_dis = math.inf
    min_index = -1
    for k in range(len(centroids)):
        dis = euclidean_distance(data_point, centroids[k])
        if dis < min_dis:
            min_dis = dis
            min_index = k
    return min_index


# returns an array (of length K) of arrays containing datapoints assigning each one to a centroid.
# i.e. assignment[i][j] means that cluster indexed i was the closest cluster to the data point at assignment[i][j].
def assign_datapoints_to_clusters(data_points, centroids):
    assignment = [[] for _ in range(len(centroids))]
    for i in range(len(data_points)):
        min_index = find_min_index(data_points[i], centroids)
        assignment[min_index] += [data_points[i]]
    return assignment


# returns an array (of length K) of updated centroids.
def update_centroids(assignment, K, d):
    updated = []
    for i in range(K):
        new_centroid = [0 for _ in range(d)]
        for j in range(len(assignment[i])):
            new_centroid = sum_vectors(new_centroid, assignment[i][j])
        new_centroid = divide_by_scalar(new_centroid, len(assignment[i]))
        updated += [new_centroid]
    return updated


def check_arguments(flag_K, flag_iter, flag_file, K, N, iter):
    flag_K = flag_K and 1 < K < N
    if flag_K:
        # print("K:", K)
        pass
    else:
        print("Invalid number of clusters!")

    flag_iter = flag_iter and 1 < iter < 1000
    if flag_iter:
        # print("iter:", iter)
        pass
    else:
        print("Invalid maximum iteration!")

    if not flag_file:
        print("NA")

    return flag_K and flag_iter and flag_file


def print_results(clusters):
    for i in range(len(clusters)):
        print(str(clusters[i]).replace(" ", "")[1:-1])


def k_means_algorithm(K, max_iter, file_path):
    flag_K = 1 < K
    flag_iter = 1 < max_iter < 1000
    flag_file = True  # Initialize flag_file here

    data_points = []
    try:
        with open(file_path) as file:
            for line in file:
                vector = line.rstrip().split(',')
                vector = [eval(i) for i in vector]
                data_points += [vector]
    except OSError:
        flag_file = False

    N = len(data_points)
    d = len(data_points[0])

    arguments_ok = flag_K and flag_iter and flag_file

    if arguments_ok:
        centroids = [data_points[i] for i in range(K)]
        iteration_number = 0
        flag_delta = compute_flag_delta(centroids, [[0 for _ in range(d)] for _ in range(K)])

        while not (iteration_number == max_iter or flag_delta):
            iteration_number += 1
            assignment = assign_datapoints_to_clusters(data_points, centroids)
            new_centroids = update_centroids(assignment, K, d)
            flag_delta = compute_flag_delta(centroids, new_centroids)
            centroids = new_centroids

        return centroids, assignment


# *** Code for analysis.py *** #

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
