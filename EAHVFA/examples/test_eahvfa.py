import sys
sys.path.append('..')
import multiprocessing
from algorithms.EAHVFA import eahvfa

# ------------------------------- Reference --------------------------------
# Liu S, Wang H, Yao W. A surrogate-assisted evolutionary algorithm with hypervolume
# triggered fidelity adjustment for noisy multiobjective integer programming[J]. Applied Soft Computing, 2022.
# ------------------------------- Copyright --------------------------------
# Copyright (c) 2022 HandingWangXD Group. Permission is granted to copy and
# use this code for research, noncommercial purposes, provided this
# copyright notice is retained and the origin of the code is cited. The
# code is provided "as is" and without any warranties, express or implied.
# ------------------------------- Developer --------------------------------
# This code is written by Shulei Liu. Email: shuleiliu@126.com

# Parameter setting
# pro, test problem, its evaluation function is defined in benchmark/problems.py
# dim, the dimension of decision variables

# pros = ['zdt1', 'zdt2', 'zdt3', 'zdt4', 'zdt6', 'dtlz1', 'dtlz2', 'dtlz3', 'dtlz4', 'dtlz5', 'dtlz6', 'dtlz7', 'wfg1', 'wfg2', 'wfg5']
# parallel
pros = ['zdt1', 'zdt4', 'dtlz1', 'dtlz2']
dims = [10, 30, 50]
for pro in pros:
    for dim in dims:
        pool = multiprocessing.Pool(processes=10)
        for i in range(10):
            pool.apply_async(eahvfa, (pro, dim))
        pool.close()
        pool.join()
    print(pro)


# not parallel
# eahvfa('zdt1', 10)