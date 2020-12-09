# Input:
# dimension            -No. of Decision Variables
# datanum              -The size of offline data
# lower_bound          -Lower bound of test function
# upper_bound          -Upper bound of test function
#
# Output:
# Execution Time
# optimum           -The final optimal solution.
#
####    Authors:    Pengfei Huang, Handing Wang, Wenping Ma
####    Xidian University, China
####    EMAIL:      pfeihuang@foxmail.com, hdwang@xidian.edu.cn
####    DATE:       December 2019
# ------------------------------------------------------------------------
# This code is part of the program that produces the results in the following paper:
#
# Pengfei Huang,Handing Wang,Wenping Ma,Stochastic Ranking for Offline Data-Driven Evolutionary Optimization Using Radial Basis Function Networks with Multiple Kernels,2019 IEEE Symposium Series on Computational Intelligence.
#
# You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
# ------------------------------------------------------------------------
# coding:utf-8
import numpy as np
import time
from RBFN import RBFN
from Latin import latin
import Test_Functions as fun
from GA import GA

if __name__ == '__main__':
    starttime = time.perf_counter()
    dimension = 10
    fun = fun.ellipsoid
    lower_bound = -5.12
    upper_bound = 5.12
    datanum = 11 * dimension

    x = latin(datanum, dimension, lower_bound, upper_bound)
    y = fun(x)
    model = [0] * 4
    model[0] = RBFN(input_shape=dimension, hidden_shape=int(np.sqrt(datanum)), kernel='gaussian')
    model[1] = RBFN(input_shape=dimension, hidden_shape=int(np.sqrt(datanum)), kernel='reflect')
    model[2] = RBFN(input_shape=dimension, hidden_shape=int(np.sqrt(datanum)), kernel='mul')
    model[3] = RBFN(input_shape=dimension, hidden_shape=int(np.sqrt(datanum)), kernel='inmul')
    for i in range(4):
        model[i].fit(x, y)

    max_iter = 100
    ga = GA(pop_size=100, dimension=dimension, lower_bound=lower_bound, upper_bound=upper_bound)
    ga.init_Population()
    for i in range(max_iter):
        ga.crossover(ga.pc)
        ga.mutation(ga.pm)
        ga.pop = np.unique(ga.pop, axis=0)
        for j in range(len(model)):
            temp = model[j].predict(ga.pop)
            if j == 0:
                fit_value = temp
            else:
                fit_value = np.column_stack((fit_value, temp))
        fit_value = fit_value.reshape((len(ga.pop), len(model)))
        ga.selection(fit_value)
    optimum = ga.pop[-1]
    endtime = time.perf_counter()
    print('Optimal solution :', optimum)
    print('Execution Time :', endtime - starttime)
