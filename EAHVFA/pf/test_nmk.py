import numpy as np
import pandas as pd

# links = np.zeros((2, 18, 2+1))
# print(links.shape)



def get_links_tables(data, n, m, k):
    links = np.zeros((m, n, k+1))
    tables = np.zeros((m, n, 2**(k+1)))

    # init links and tables
    for mi in range(m):
        t = 0
        q = n * (k + 1)
        for ni in range(n):
            links[mi, ni, 0:(k+1)] = data[t:t+k+1, mi]
            tables[mi, ni, 0:2 ** (k + 1)] = data[q:q + 2 ** (k + 1), mi]
            t = t + k + 1
            q = q + 2 ** (k + 1)
    return links, tables

def eva_nmk_sep(links, tables, mi, ind):
    acc = 0.0
    n = links.shape[1]
    for i in range(n):
        acc = acc + tables[mi, i, sigma(links, mi, ind, i)]
    return acc / float(n)


def sigma(links, mi, ind, ni):
    n = 1
    acc = 0
    k = links.shape[2] - 1
    for j in range(k+1):
        if ind[int(links[mi, ni, j])] == 0:
            acc = acc | n
        n = n << 1
    return int(acc)

def eva_nmk(links, tables, ind):
    r1 = eva_nmk_sep(links, tables, 0, ind)
    r2 = eva_nmk_sep(links, tables, 1, ind)

    print(r1, r2)



data = pd.read_csv('rmnk_64.csv', header=None, sep=',').values
links, tables = get_links_tables(data, 64, 2, 2)
ind = [0,1,0,0,0,0,1,0,1,1,1,1,1,0,1,1,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,1,1,0,0,1,1,1,0,1,1,1,0,1,1,0]
eva_nmk(links, tables, ind)
