# Input:
# dimension            -No. of Decision Variables
# datanum              -The size of offline data
# lower_bound          -Lower bound of test function
# upper_bound          -Upper bound of test function
#
# Output:
# Execution Time
# optimum           -The final optimal solution.
#### Authors: Pengfei Huang, Handing Wang, Yaochu Jin
#### Xidian University, China and University of Surrey, United Kingdom.
#### EMAIL: pfeihuang @ foxmail.com, hdwang @ xidian.edu.cn
#### WEBSITE: https://sites.google.com/site/handingwanghomepage
#### DATE: December 2020
# ------------------------------------------------------------------------
# This code is part of the program that produces the results in the following paper:
#
# Pengfei Huang,Handing Wang,Yaochu Jin,Offine Data-Driven Evolutionary Optimization Based on Tri-Training, Swarm and Evolutionary Computation, Accepted.
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

def updatemodel(data):
    global numx0, numy0, numx1, numy1, numx2, numy2
    pre0 = model[0].predict(data)
    pre1 = model[1].predict(data)
    pre2 = model[2].predict(data)

    error = abs(pre1-pre2)
    seq = np.ravel(np.where(error == np.min(error)))[0]
    xtemp = np.row_stack((numx0, data[seq]))
    ytemp = np.append(numy0, (pre1[seq]+pre2[seq])/2)
    model[0].fit(xtemp, ytemp)
    # print('model0update')
    error = abs(pre0-pre2)
    seq = np.ravel(np.where(error == np.min(error)))[0]
    xtemp = np.row_stack((numx1, data[seq]))
    ytemp = np.append(numy1, (pre0[seq]+pre2[seq])/2)
    model[1].fit(xtemp, ytemp)
    # print('model1update')
    error = abs(pre0-pre1)
    seq = np.ravel(np.where(error == np.min(error)))[0]
    xtemp = np.row_stack((numx2, data[seq]))
    ytemp = np.append(numy2, (pre0[seq]+pre1[seq])/2)
    model[2].fit(xtemp, ytemp)
    # print('model2update')

def resetmodel(x,y):
    global numx0, numx1, numx2, numy0, numy1, numy2
    shuffledata = np.column_stack((y, x))
    np.random.shuffle(shuffledata)
    newx = shuffledata[:, 1:]
    newy = shuffledata[:, :1]
    numx0 = newx[:traindata, ]
    numy0 = newy[:traindata, ]
    numx1 = newx[traindata:2 * traindata, ]
    numy1 = newy[traindata:2 * traindata, ]
    numx2 = newx[datanum - traindata:, ]
    numy2 = newy[datanum - traindata:, ]

    model[0].fit(numx0, numy0)
    model[1].fit(numx1, numy1)
    model[2].fit(numx2, numy2)

if __name__ == '__main__':
    starttime = time.perf_counter()
    dimension = 10
    fun = fun.ellipsoid
    lower_bound = -5.12
    upper_bound = 5.12
    datanum = 11 * dimension

    x = latin(datanum, dimension, lower_bound, upper_bound)
    y = fun(x)
    model = [0] * 3
    traindata = int(datanum / 3)
    model[0] = RBFN(input_shape=dimension, hidden_shape=int(np.sqrt(traindata)), kernel='gaussian')
    model[1] = RBFN(input_shape=dimension, hidden_shape=int(np.sqrt(traindata)), kernel='gaussian')
    model[2] = RBFN(input_shape=dimension, hidden_shape=int(np.sqrt(traindata)), kernel='gaussian')
    resetmodel(x,y)

    max_iter = 100
    ga = GA(pop_size=100, dimension=dimension, lower_bound=lower_bound, upper_bound=upper_bound)
    ga.init_Population()
    for i in range(max_iter):
        updatemodel(ga.pop)
        ga.crossover(ga.pc)  
        ga.mutation(ga.pm) 
        ga.pop = np.unique(ga.pop, axis=0)
        for j in range(0, 3):
            temp = model[j].predict(ga.pop)
            if j == 0:
                fit_value = temp
            else:
                fit_value = fit_value + temp
        fit_value = fit_value.reshape((len(ga.pop), 1))
        ga.selection(fit_value)  
        resetmodel(x,y)

    optimum = ga.first[-1]
    endtime = time.perf_counter()
    print('Optimal solution :', optimum)
    print('Execution Time :', endtime - starttime)
