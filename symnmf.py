import sys
import numpy as np
import math
import csv
import symnmfmodule


def print_matrix(matrix):
    str1 = ''
    for row in matrix:
        str1 = ','.join(map(lambda x: "%.4f" % x, row))
        print(str1)


# note in 2.9.2 in project_NMF.pdf we are told we don't need to validate args.
def algorithm():
    np.random.seed(0)
    k = int(sys.argv[1])
    goal = sys.argv[2]
    file_path = sys.argv[3]
    data_points = np.genfromtxt(file_path, dtype=float, encoding=None, delimiter=",")  # A two-dimensional numpy array
    data_points = data_points.tolist() 
    

    result_mat = None
    if goal == 'sym':
        result_mat = symnmfmodule.sym(data_points)
    elif goal == 'ddg':
        result_mat = symnmfmodule.ddg(data_points)
        n = len(result_mat)
        result_mat = [[0 if i != j else result_mat[i] for j in range(n)] for i in range(n)]
    elif goal == 'norm':
        result_mat = symnmfmodule.norm(data_points)
    else: # goal == 'symmnf'
        W = symnmfmodule.norm(data_points)
        n = len(W)
        m = sum(map(sum, W)) / sum(map(len, W)) # average of entries of W.
        r = 2 * math.sqrt(m) / math.sqrt(k)
        H = np.random.uniform(low=0, high=r, size=(n,k))
        H = H.tolist()
        result_mat = symnmfmodule.symnmf(W, H, n, k)
    print_matrix(result_mat)

algorithm()


