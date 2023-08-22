import sys
import numpy as np
import symnmfmodule


# note in 2.9.2 in project_NMF.pdf we are told we don't need to validate args.
def algorithm():
    np.random.seed(0)
    k = int(sys.argv[1])
    goal = sys.argv[2]
    file_name = sys.argv[3]
    X = ...

    if goal == 'sym':
        A = symnmfmodule.sym(X)
    elif goal == 'ddg':
        D = symnmfmodule.ddg(X)
    elif goal == 'norm':
        W = symnmfmodule.norm(X)
    else: # goal == 'symmnf'
        H = symnmfmodule.symnmf(X)


algorithm()
