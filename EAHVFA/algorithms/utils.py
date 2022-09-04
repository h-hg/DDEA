import numpy as np
import random
import copy
from benchmark.problems import get_fitness

# ------------------------------- Reference --------------------------------
# Liu S, Wang H, Yao W. A surrogate-assisted evolutionary algorithm with hypervolume
# triggered fidelity adjustment for noisy multiobjective integer programming. Applied Soft Computing, 2022.
# ------------------------------- Copyright --------------------------------
# Copyright (c) 2022 HandingWangXD Group. Permission is granted to copy and
# use this code for research, noncommercial purposes, provided this
# copyright notice is retained and the origin of the code is cited. The
# code is provided "as is" and without any warranties, express or implied.
# ------------------------------- Developer --------------------------------
# This code is written by Shulei Liu. Email: shuleiliu@126.com


def sample(pro, size, dim, fidelity):
    if pro == 'uf1' or pro == 'uf2' or pro == 'uf6':
        x1_low = [0]
        x1_up = [1]
        xor_low = np.zeros(dim-1) - 1
        xor_up = np.zeros(dim-1) + 1
        x_low = x1_low + list(xor_low)
        x_up = x1_up + list(xor_up)
    elif pro == 'uf4':
        x1_low = [0]
        x1_up = [1]
        xor_low = np.zeros(dim - 1) - 2
        xor_up = np.zeros(dim - 1) + 2
        x_low = x1_low + list(xor_low)
        x_up = x1_up + list(xor_up)
    elif pro == 'uf8' or pro == 'uf9':
        x1_low = [0, 0]
        x1_up = [1, 1]
        xor_low = np.zeros(dim - 2) - 2
        xor_up = np.zeros(dim - 2) + 2
        x_low = x1_low + list(xor_low)
        x_up = x1_up + list(xor_up)
    elif pro == 'zdt4':
        x1_low = [0]
        x1_up = [1]
        xor_low = np.zeros(dim - 1) - 5
        xor_up = np.zeros(dim - 1) + 5
        x_low = x1_low + list(xor_low)
        x_up = x1_up + list(xor_up)
    else:
        x_low = list(np.zeros(dim))
        x_up = list(np.zeros(dim)+1)


    x = []
    for i in range(size):
        ind = []
        for i in range(dim):
            ind.append(random.sample(range(int(x_low[i]*100), int(x_up[i]*100), 2), 1)[0] / 100)
        x.append(ind)

    # 真实评价
    f1 = []
    f2 = []
    f3 = []

    for xi in x:
        if pro[0:3] == 'zdt' or pro[0:3] == 'wfg' or (pro[0:3] in ['uf1', 'uf2', 'uf4', 'uf6']):
            y1, y2 = get_fitness(xi, pro, fidelity)
            f1.append(y1)
            f2.append(y2)
        else:
            y1, y2, y3 = get_fitness(xi, pro, fidelity)
            f1.append(y1)
            f2.append(y2)
            f3.append(y3)


    return x, f1, f2, f3

# generate the first neighborhood according to ind
def neighborhood_one(ind, x_1, x_other, n_size):
    x_1 = list(x_1)
    x_other = list(x_other)
    ns = []
    m = random.sample(range(0, len(ind)), n_size)
    for mi in m:
        s = copy.deepcopy(ind)
        value = s[mi]
        while value == s[mi]:
            if mi == 0:
                value = random.sample(x_1, 1)[0]
            else:
                value = random.sample(x_other, 1)[0]
        s[mi] = value
        ns.append(copy.deepcopy(s))
    return ns

# generate the second neighborhood according to ind
def neighborhood_two(ind, x_1, x_other, n_size):
    x_1 = list(x_1)
    x_other = list(x_other)
    ns = []
    m1 = random.sample(range(0, len(ind)), n_size)
    m2 = random.sample(range(0, len(ind)), n_size)
    for mi in range(n_size):
        s = copy.deepcopy(ind)
        while m1[mi] == m2[mi]:
            m1[mi] = random.sample(range(0, len(s)), 1)[0]
        value1 = s[m1[mi]]
        while value1 == s[m1[mi]]:
            if m1[mi] == 0:
                value1 = random.sample(x_1, 1)[0]
            else:
                value1 = random.sample(x_other, 1)[0]
        s[m1[mi]] = value1
        value2 = s[m2[mi]]
        while value2 == s[m2[mi]]:
            if m2[mi] == 0:
                value2 = random.sample(x_1, 1)[0]
            else:
                value2 = random.sample(x_other, 1)[0]
        s[m2[mi]] = value2
        ns.append(copy.deepcopy(s))
    return ns


# randomlt select a solution from a neighborhood ns
def random_select_from_one_set(ns):
    if len(ns) == 0:
        return []
    else:
        idx = np.random.randint(0, len(ns))
        s_t = copy.deepcopy(ns[idx])
    return s_t

# The method to calculated IGD
# p the set of reference points
# Q the set of optimal solutions
def IGD(P,Q):
    size = P.shape[0]
    sum_dis = 0
    for i in range(size):
        min_dis = 0
        n = np.mat(P[i,])
        for j in range(Q.shape[0]):
            m = np.mat(Q[j,])
            dis = np.sqrt(np.sum(np.square(m-n)))
            if j == 0:
                min_dis = dis
            elif dis < min_dis:
                min_dis = dis
        sum_dis = sum_dis + min_dis
    return sum_dis / size

# The method to calculated IGD plus
# r the set of reference points
# s the set of optimal solutions
def IGD_Plus(r,s):
    size = r.shape[0]
    sum_dis = 0
    for i in range(size):
        min_dis = np.inf
        for j in range(s.shape[0]):
            plus_dis = 0
            for k in range(s.shape[1]):
                if s[j, k] - r[i, k] > 0:
                    plus_dis = plus_dis + ((s[j, k] - r[i, k]) **2)
                else:
                    plus_dis = plus_dis + 0
            plus_dis = np.sqrt(plus_dis)
            if plus_dis < min_dis:
                min_dis = plus_dis
        sum_dis = sum_dis + min_dis
    return sum_dis / size
