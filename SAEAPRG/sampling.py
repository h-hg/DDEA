#    Authors: Shulei Liu, Handing Wang, Wei Peng, Wen Yao
#    Xidian University, China
#    Defense Innovation Institute, Chinese Academy of Military Science, China
#    EMAIL: shuleiliu@126.com, hdwang@xidian.edu.cn
#    WEBSITE: https://sites.google.com/site/handingwanghomepage
#    DATE:  February 2022
# ------------------------------------------------------------------------
# This code is part of the program that produces the results in the following paper:
#
# Shulei Liu, Handing Wang, Wei Peng, Wen Yao, A Surrogate-Assisted Evolutionary Feature Selection Algorithm with Parallel Random Grouping for High-Dimensional Classification, IEEE Transactions on Evolutionary Computation, 2022.
#
# You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
# ------------------------------------------------------------------------

import pandas as pd
import numpy as np
import random
import basic_operators as bo

# Sampling evaluation data to train surrogate models
# name, the name of datasets
# 采集数据的大小
def sample(name):
    # 1、读取数据集
    # print_data = []
    file_name  = 'data/' + name + '.csv'
    data = pd.read_csv(file_name, sep=",", header=None)
    data = data.values
    if data.shape[1]-1 < 5000:
        size = 2000
    elif data.shape[1]-1 >= 5000 and data.shape[1]-1 < 10000:
        size = 4000
    else:
        size = 5000
    pop = cbsm(size, data.shape[1]-1)
    # pop = bi.random_sample(size, data.shape[1]-1)
    # pop = bi.boundary_sample(size, data.shape[1] - 1)
    for k in range(len(pop)):
        error_rate = bo.get_error_rate(data, pop[k], pro=name)
        pop[k].append(error_rate)
        # if k % 200 == 0:
        #     print(k)
    return np.array(pop)


# constraint based sampling method
# size, the total number of samples
# cols, the dimension
def cbsm(size, cols):
    pop = []
    c = int(cols*0.15)
    p = c / cols
    n1 = int(size*0.05)
    n2 = size - n1
    random_number = np.abs(np.random.normal(0, 0.5, size)) + 0.5
    for i in range(n1):
        ind = []
        for k in range(cols):
            r = np.random.random()
            if r <= p:
                ind.append(1)
            else:
                ind.append(0)
        pop.append(ind)
    for j in range(n2):
        ind = []
        for m in range(cols):
            r = np.random.random()
            if r <= 0.5:
                ind.append(0)
            else:
                ind.append(1)
        sum = np.sum(ind)
        if sum > c:
            ones = np.where(np.array(ind) == 1)[0]
            rn = random_number[j]
            if rn >= 1:
                rn = np.random.random()
            u = sum - c + c * rn
            number = np.random.randint(sum - c, int(u + 1))
            deleted_idx = random.sample(ones.tolist(), number)
            for did in deleted_idx:
                ind[did] = 0
        pop.append(ind)
    return pop
