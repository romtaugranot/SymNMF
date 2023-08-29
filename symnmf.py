import sys
import numpy as np
import csv
import symnmfmodule


# note in 2.9.2 in project_NMF.pdf we are told we don't need to validate args.
def algorithm():
    np.random.seed(0)
    k = int(sys.argv[1])
    goal = sys.argv[2]
    file_path = sys.argv[3]

    if goal == 'sym':
        A = symnmfmodule.sym(file_path)
    elif goal == 'ddg':
        D = symnmfmodule.ddg(file_path)
    elif goal == 'norm':
        W = symnmfmodule.norm(file_path)
    else: # goal == 'symmnf'
        H = symnmfmodule.symnmf(file_path)


algorithm()
