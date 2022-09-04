import numpy as np
import pandas as pd

def get_q(data, n, m):
    q = np.zeros((m, n, n))

    # init links and tables
    for mi in range(m):
        t = 0
        for ni in range(n):
            q[mi, ni, 0:n] = data[t:t+n, mi]
            t = t + n
    return q

def eva_ubqp_sep(q, mi, ind):
    n = q.shape[1]
    fit = 0
    for i in range(n):
        if ind[i] == 0:
            for j in range(0, i+1):
                if ind[j] == 0:
                    fit = fit + q[mi, i, j]
    return fit

def eva_ubqp(data, ind):
    m = data.shape[1]
    n = int(np.sqrt(data.shape[0]))
    q = get_q(data, n, m)
    r1 = eva_ubqp_sep(q, 0, ind)
    r2 = eva_ubqp_sep(q, 1, ind)
    print(r1, r2)


data = pd.read_csv('mubqp_50.csv', header=None, sep=',').values
ind = [1,0,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,1,0,0,1,1,0,0,1,1,1,1,0,0,1,0,1,1,0,1,0,0,0,1,1]
eva_ubqp(data, ind)