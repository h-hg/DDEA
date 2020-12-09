#**************************************************************************************************
#Last Edited: 12/10/2020
#**************************************************************************************************
# Implementation of Classifier-Assisted Level-Assisted Learning Swarm Optimizer (CA-LLSO)
#
# Feng-Feng Wei, Wei-Neng Chen, Qiang Yang, Jeremiah Deng, Xiao-Nan Luo, Hu Jin and Jun Zhang, 
# "A Classifier-Assisted Level-Based Learning Swarm Optimizer for Expensive Optimization", 
# IEEE Transactions on Evolutionary Computation, accepted, 2020.
#
# This paper and this code should be referenced whenever they are used to
# generate results for the user's own research.
#
# This matlab code was written by Feng-Feng Wei
# School of Computer Science and Engineering, South China University of Technology
#
# Please refer with all questions, comments, bug reports, etc. to fengfeng_scut@163.com
#
#**************************************************************************************************
#
# REQUIREMENTS:
# The following libraries should be intalled before running the code:
#     python v2.7.15 or newer versions
#     numpy v1.15.4 or new versions
#     scipy v1.1.0 or newer versions
#     sklearn v0.20.3
#     pyDOE v0.3.8
#
#**************************************************************************************************

import numpy as np
import time
import math
import random
from sklearn.ensemble import GradientBoostingClassifier
from pyDOE import *
from numpy import linalg as LA
#import matplotlib.pyplot as plt
np.set_printoptions(threshold=np.inf)

fi = 5 # the test problem 0 1 2 3 4 5
D = 200 # test dimensionality

if fi == 0:
    bound = 5.12 # the boundary of the test problem
elif fi == 1:
    bound = 2.048
elif fi == 2:
    bound = 32.768
elif fi == 3:
    bound = 600
elif fi == 4 or fi == 5:
    bound = 5

if D>=50:
    NP = 200 # the population size
    count_l1 = 150
elif D<50:
    NP = 100
    count_l1 = 75
max_gen = 250 # the maxinal evolution generation
max_times = 6 # the total independent run times
NL = 4 # the number of layers
phi = 0.4
pop_x = None
pop_v = np.zeros((NP, D))
pop_y = None
db_x = None
db_v = np.zeros((NP,D))
db_y = None
rank = np.zeros((NP, ))
ranked_x = np.zeros((NP, D))
ranked_v = np.zeros((NP, D))
ranked_y = np.zeros((NP, ))

offspring_rank = np.zeros((NP, ))
generated_level_1 = np.zeros((NP, D))
generated_level_1_v = np.zeros((NP, D))
eval_per_gen = int(1000/max_gen)

O_ori = None
MD = np.empty((NP*10, D))

def Ellipsoid(p):#5.12 f123 D203050
    d = np.arange(0, D) + 1
    y = np.square(p)
    return np.dot(y, d).transpose()

def Rosenbrock(p):#2.048 f456
    y = np.sum(100 * np.square(p[:,1:] - np.square(p[:,:-1])) + np.square(p[:,:-1] - 1), axis = 1)
    return y

def Ackley(p):#32.768 f789  f12
    tmp1 = np.sqrt(np.sum(np.square(p), axis = 1) / D)
    tmp2 = np.sum(np.cos(2 * np.pi * p), axis = 1)
    y = -20 * np.exp(-0.2 * tmp1) - np.exp(tmp2 / D) + 20 + np. exp(1)
    return y

def Griewank(p):#600 f101112
    tmp1 = np.sum(np.square(p), axis = 1) / 4000
    tmp2 = np.cos(p / np.sqrt(np.arange(0, D) + 1))
    y = 1 + tmp1 - np.prod(tmp2, axis = 1)
    return y

def ShiftedRotatedRastrigin(p):#5 f13 D30 
    O_or = O_ori[0:D]
    O = O_or.reshape(1,D)
    z = np.dot((p - O), MD30)
    y = np.sum((np.square(z) - 10 * np.cos(2 * np.pi * p) + 10), axis = 1) - 330
    return y

def LAckley(p):#32.768 f12
    tmp1 = np.sqrt(np.sum(np.square(p), axis = 1) / D)
    tmp2 = np.sum(np.cos(2 * np.pi * p), axis = 1)
    y = -20 * np.exp(-0.2 * tmp1) - np.exp(tmp2 / D) + 20 + np. exp(1)
    return y

def LRastrigin(p):#5 f34 
    y = np.sum((np.square(p) - 10 * np.cos(2 * np.pi * p) + 10), axis = 1)
    return y

def LSphere(p):#f56
    y = np.sum(np.square(p), axis=1)
    return y

def LWeierstrass(p):#f78
    ps, D = p.shape
    kmax = 20
    a = 0.5
    b = 3
    y = []
    for i in range(p.shape[0]):
        acc = 0
        for j in range(D):
            for k in range(kmax):
                acc += a**k + math.cos(2*math.pi*b**k*(p[i, j] + 0.5))
                acc -= a**k * math.cos(2*math.pi*b**k*0.5)
        y.append(acc)
    return np.array(y)

def LGriewank(p):#600 f910
    tmp1 = np.sum(np.square(p), axis = 1) / 4000
    tmp2 = np.cos(p / np.sqrt(np.arange(0, D) + 1))
    y = 1 + tmp1 - np.prod(tmp2, axis = 1)
    return y

def RHCF(p): #[-5, 5] f_bias = 10
    final_y = 0
    subfunctions = [LAckley, LAckley, LRastrigin, LRastrigin, LSphere, LSphere, LWeierstrass, LWeierstrass, LGriewank, LGriewank]
    bias = 10
    sigma = np.array([0.1, 2, 1.5, 1.5, 1, 1, 1.5, 1.5, 2, 2])
    lambd = np.array([0.1 * 5 / 32, 5 / 32, 2 * 1, 1, 2 * 5 / 100, 5 / 100, 2 * 10, 10, 2 * 5 / 60, 5 / 60])
    lambd = np.repeat(lambd.reshape(10,1), D, axis = 1)
    O = O_ori[:, 0:D]
    bias = [0,100,200,300,400,500,600,700,800,900]
    O[9,:] = 0
    weight = np.empty((p.shape[0],10))
    for i in range(10):
        oo = np.repeat(O[i,:].reshape(1, D), p.shape[0], axis = 0)
        weight[:,i] = np.exp(-np.sum(np.square(p - oo), axis = 1) / (2*D*sigma[i]*sigma[i]))
    maxweight = np.amax(weight, axis = 1)
    for i in range(p.shape[0]):
        for j in range(10):
            if weight[i,j] != maxweight[i]:
                weight[i,j] = weight[i,j] * (1 - np.power(maxweight[i], 10))
    sumweight = np.sum(weight, axis = 1)
    weight = weight / np.repeat(sumweight.reshape(-1, 1), 10, axis=1)

    tmp_y = np.ones((1,D)) * 5
    for i in range (10):
        oo = np.repeat(O[i,:].reshape(1,D), p.shape[0], axis = 0)
        cur_subfunc = subfunctions[i]
        f = cur_subfunc(np.matmul((p - oo) / (np.repeat(lambd[i,:].reshape(1, D), p.shape[0], axis = 0)), MD[i * D:(i+1) * D,:]))
        fmax = cur_subfunc(np.matmul(tmp_y / lambd[i,:].reshape(1, D), MD[i * D:(i+1) * D,:]))
        f1 = 2000 * f / fmax
        final_y = final_y + weight[:,i] * (f1 + bias[i])
    return final_y

def rot_matrix(D, c):
    A = np.random.normal(0, 1, (D, D))
    P, _ = LA.qr(A)
    A = np.random.normal(0, 1, (D, D))
    Q, _ = LA.qr(A)
    u = np.random.random((1,D))
    D = np.power(c, (u-np.min(u)) / (np.max(u) - np.min(u)))
    D = np.squeeze(D)
    D = np.diag(D)
    M = np.dot(np.dot(P,D),Q)
    return M

def update(ori_x, ori_v, x_1, x_2):
    r1 = np.random.random((1,D))
    r2 = np.random.random((1,D))
    r3 = np.random.random((1,D))
    new_v = r1 * ori_v + r2 * (x_1 - ori_x) + phi * r3 * (x_2 - ori_x)
    new_x = ori_x + new_v
    for i in range(D):
        if new_x[0,i] < -bound:
            new_x[0,i] = -bound
        if new_x[0,i] > bound:
            new_x[0,i] = bound
    return new_x, new_v

functions1 = [Ellipsoid, Rosenbrock, Ackley, Griewank, ShiftedRotatedRastrigin, RHCF]

if __name__ == "__main__":
    selected_function = functions1[fi]
    if fi == 4:
        if D <= 100:
            O_ori = np.loadtxt('rastrigin_func_data.txt')
        else:
            O_ori = -5 + 10 * np.random.random((1,D))
        if D == 10:
            MD30 = np.loadtxt('rastrigin_M_D10.txt')
        elif D == 30:
            MD30 = np.loadtxt('rastrigin_M_D30.txt')
        elif D == 50:
            MD30 = np.loadtxt('rastrigin_M_D50.txt')
        else:
            c = 2
            MD30 = rot_matrix(D, c)
    if fi == 5:
        if D <= 100:
            O_ori = np.loadtxt('hybrid_func2_data.txt')
        else:
            O_ori = -5 + 10 * np.random.random((10,D))
        if D == 10:
            MD = np.loadtxt('hybrid_func2_M_D10.txt')
        elif D == 30:
            MD = np.loadtxt('hybrid_func2_M_D30.txt')
        elif D == 50:
            MD = np.loadtxt('hybrid_func2_M_D50.txt')
        else:
            c = [2, 3, 2, 3, 2, 3, 20, 30, 200, 300]
            #MD = rot_matrix(D, c)
            MD = np.empty((10 * D, D))
            for i in range(10):
                MD[i * D:(i+1) * D,:] = rot_matrix(D, c[i])

    start_time = time.time()
    best_ever = np.zeros((max_times, max_gen))
    con_count = np.zeros((max_times, max_gen))
    for times in range(max_times):
        #initialization
        db_x = lhs(D, samples=NP, criterion = 'center') * bound * 2 - bound
        db_v = np.zeros((NP,D))
        db_y = selected_function(db_x)
        #iteration
        for itera in range(max_gen): 
            ORIGIN = np.concatenate((db_x,db_v,db_y.reshape(db_y.shape[0],1)),axis=1)
            ORIGIN = ORIGIN[ORIGIN[:,-1].argsort()]
            pop_x = ORIGIN[0:NP,0:D]
            pop_v = ORIGIN[0:NP,D:D*2]
            pop_y = ORIGIN[0:NP,-1]
            level1_count = 0
            #sort_y = np.argsort(pop_y)
            #print(sort_y)
            for i in range(NL-1):
                ranked_x[int(NP/NL)*i:int(NP/NL)*(i+1),:] = pop_x[int(NP/NL)*i:int(NP/NL)*(i+1),:]
                ranked_v[int(NP/NL)*i:int(NP/NL)*(i+1),:] = pop_v[int(NP/NL)*i:int(NP/NL)*(i+1),:]
                ranked_y[int(NP/NL)*i:int(NP/NL)*(i+1)] = pop_y[int(NP/NL)*i:int(NP/NL)*(i+1), ]
                rank[int(NP/NL)*i:int(NP/NL)*(i+1)] = i
            ranked_x[int(NP/NL)*(NL-1):,:] = pop_x[int(NP/NL)*(NL-1):,:]
            ranked_v[int(NP/NL)*(NL-1):,:] = pop_v[int(NP/NL)*(NL-1):,:]
            ranked_y[int(NP/NL)*(NL-1):,] = pop_y[int(NP/NL)*(NL-1):,]
            rank[int(NP/NL)*(NL-1):] = NL - 1
            clf = GradientBoostingClassifier().fit(ranked_x, rank)
            #update
            offspring_x = np.zeros((NP, D))
            offspring_v = np.zeros((NP, D))
            offspring_rank = np.zeros((NP, ))
            for i in range(NP-(NL-1)*int(NP/NL)):
                rl1 = np.random.randint(0, NL-1)
                rl2 = np.random.randint(0, NL-1)
                while rl1 == rl2:
                    rl2 = np.random.randint(0, NL-1)
                if rl1 > rl2:
                    tmp = rl1
                    rl1 = rl2
                    rl2 = tmp
                k1 = np.random.randint(0, int(NP/NL))
                k2 = np.random.randint(0, int(NP/NL))
                lev_1 = rl1 * int(NP/NL) + k1
                lev_2 = rl2 * int(NP/NL) + k2
                offspring_x[-(i+1),:], offspring_v[-(i+1),:] = update(ranked_x[-(i+1),:], ranked_v[-(i+1),:], ranked_x[lev_1,:], ranked_x[lev_2,:])
                offspring_rank[-(i+1),] = clf.predict(offspring_x[-(i+1),:].reshape(1, D))
                if level1_count < count_l1 and offspring_rank[-(i+1),] == 0:
                    generated_level_1[level1_count,:] = offspring_x[-(i+1),:]
                    generated_level_1_v[level1_count,:] = offspring_v[-(i+1),:]
                    level1_count = level1_count + 1

            #update the third to the last-1 level
            for i in range(2, NL - 1):
                for j in range(int(NP/NL)):
                    rl1 = np.random.randint(0, i)
                    rl2 = np.random.randint(0, i)
                    while rl1 == rl2:
                        rl2 = np.random.randint(0, i)
                    if rl1 > rl2:
                        tmp = rl1
                        rl1 = rl2
                        rl2 = tmp
                    k1 = np.random.randint(0, int(NP/NL))
                    k2 = np.random.randint(0, int(NP/NL))
                    cur = i * int(NP/NL) + j
                    lev_1 = rl1 * int(NP/NL) + k1
                    lev_2 = rl2 * int(NP/NL) + k2

                    offspring_x[cur,:], offspring_v[cur,:] = update(ranked_x[cur,:], ranked_v[cur,:], ranked_x[lev_1,:], ranked_x[lev_2,:])
                    offspring_rank[cur,] = clf.predict(offspring_x[cur,:].reshape(1, D))
                    if level1_count < count_l1 and offspring_rank[cur,] == 0:
                        generated_level_1[level1_count,:] = offspring_x[cur,:]
                        generated_level_1_v[level1_count,:] = offspring_v[cur,:]
                        level1_count = level1_count + 1

            #update the second level
            for i in range(int(NP/NL)):
                r1 = np.random.randint(0, int(NP/NL))
                r2 = np.random.randint(0, int(NP/NL))
                while r1 == r2:
                    r2 = np.random.randint(0, int(NP/NL))
                if r1 > r2:
                    tmp = r1
                    r1 = r2
                    r2 = tmp
                j = i + int(NP/NL)
                offspring_x[j,:], offspring_v[j,:] = update(ranked_x[j,:], ranked_v[j,:], ranked_x[r1,:], ranked_x[r2,:])
                offspring_rank[j,] = clf.predict(offspring_x[j,:].reshape(1, D))
                if level1_count < count_l1 and offspring_rank[j,] == 0:
                    generated_level_1[level1_count,:] = offspring_x[j,:]
                    generated_level_1_v[level1_count,:] = offspring_v[j,:]
                    level1_count = level1_count + 1
            #print(offspring_x.shape,offspring_v.shape,offspring_rank.shape)
            #copy the first level
            offspring_x[0:int(NP/NL),:] = ranked_x[0:int(NP/NL),:]
            offspring_v[0:int(NP/NL),:] = ranked_v[0:int(NP/NL),:]
            offspring_rank[0:int(NP/NL),] = 0
            # L1-exploitation
            while level1_count < count_l1:
                con_count[times,itera] = con_count[times,itera]+1
                #print(offspring_x.shape,offspring_v.shape,offspring_rank.shape)
                OFF = np.concatenate((offspring_x,offspring_v,offspring_rank.reshape(offspring_rank.shape[0],1)),axis=1)
                OFF = OFF[OFF[:,-1].argsort()]
                offspring_x = OFF[:,0:D]
                offspring_v = OFF[:,D:D*2]
                offspring_rank = OFF[:,-1]
                nl0 = np.sum(offspring_rank==0)
                nl1 = np.sum(offspring_rank==1)
                nl2 = np.sum(offspring_rank==2)
                nl3 = np.sum(offspring_rank==3)
                nl = [nl0, nl0+nl1, nl0+nl1+nl2, nl0+nl1+nl2+nl3]
                #print(level1_count, nl)
                #update the last level
                if nl3 > 0 and level1_count < count_l1:
                    for i in range(nl3):
                        k1 = np.random.randint(0, nl0+nl1+nl2)
                        k2 = np.random.randint(0, nl0+nl1+nl2)
                        while k1 == k2:
                            k2 = np.random.randint(0, nl0+nl1+nl2)
                        if (offspring_rank[k1] == offspring_rank[k2]) and nl[int(offspring_rank[k1])] == (NP-nl3):
                            lev_1 = k1
                            lev_2 = k2
                        else:
                            while offspring_rank[k1] == offspring_rank[k2]:
                                k2 = np.random.randint(0, nl0+nl1+nl2)
                            if offspring_rank[k1] > offspring_rank[k2]:
                                lev_1 = k2
                                lev_2 = k1
                            else:
                                lev_1 = k1
                                lev_2 = k2
                        offspring_x[-(i+1),:], offspring_v[-(i+1),:] = update(offspring_x[-(i+1),:], offspring_v[-(i+1),:], offspring_x[lev_1,:], offspring_x[lev_2,:])
                        offspring_rank[-(i+1),] = clf.predict(offspring_x[-(i+1),:].reshape(1, D))
                        if level1_count < count_l1 and offspring_rank[-(i+1),] == 0:
                            generated_level_1[level1_count,:] = offspring_x[-(i+1),:]
                            generated_level_1_v[level1_count,:] = offspring_v[-(i+1),:]
                            level1_count = level1_count + 1
                #update the third to the last-1 level
                for i in range(2, 3):
                    if nl2 > 0 and level1_count < count_l1:
                        for j in range(nl2):
                            cur = nl0+nl1
                            cur = cur + j
                            if nl1 == 0:
                                k1 = np.random.randint(0, nl0)
                                k2 = np.random.randint(0, nl0)
                                while r1 == r2:
                                    r2 = np.random.randint(0, nl0)
                                lev_1 = k1
                                lev_2 = k2
                            else:
                                k1 = np.random.randint(0, nl0)
                                k2 = np.random.randint(0, nl1)
                                lev_1 = k1
                                lev_2 = k2+nl0

                            offspring_x[cur,:], offspring_v[cur,:] = update(offspring_x[cur,:], offspring_v[cur,:], offspring_x[lev_1,:], offspring_x[lev_2,:])
                            offspring_rank[cur,] = clf.predict(offspring_x[cur,:].reshape(1, D))
                            if level1_count < count_l1 and offspring_rank[cur,] == 0:
                                generated_level_1[level1_count,:] = offspring_x[cur,:]
                                generated_level_1_v[level1_count,:] = offspring_v[cur,:]
                                level1_count = level1_count + 1
                            
                #update the second level
                if nl1 > 0 and level1_count < count_l1:
                    for i in range(nl1):
                        r1 = np.random.randint(0, nl0)
                        r2 = np.random.randint(0, nl0)
                        while r1 == r2:
                            r2 = np.random.randint(0, nl0)
                        j = i + nl0
                        offspring_x[j,:], offspring_v[j,:] = update(offspring_x[j,:], offspring_v[j,:], offspring_x[r1,:], offspring_x[r2,:])
                        offspring_rank[j,] = clf.predict(offspring_x[j,:].reshape(1, D))
                        if level1_count < count_l1 and offspring_rank[j,] == 0:
                            generated_level_1[level1_count,:] = offspring_x[j,:]
                            generated_level_1_v[level1_count,:] = offspring_v[j,:]
                            level1_count = level1_count + 1
            dist = np.zeros((count_l1, int(NP/NL)))
            for i in range(count_l1):
                for j in range(int(NP/NL)):
                    dist[i,j] = np.sqrt(np.sum(np.power(generated_level_1[i,:]-ranked_x[j,:], 2)))

            max_dist = np.amax(dist, axis=1)
            #max_dist = np.sum(dist, axis=1) / int(NP/NL)
            mini_max_dist_ind = np.argsort(max_dist)
            #print(mini_max_dist_ind)
            x_tail = np.zeros((eval_per_gen, D))
            v_tail = np.zeros((eval_per_gen, D))
            x_tail[0:eval_per_gen-1,:] = generated_level_1[mini_max_dist_ind[0:eval_per_gen-1],:]
            v_tail[0:eval_per_gen-1,:] = generated_level_1_v[mini_max_dist_ind[0:eval_per_gen-1],:]
            rand_num1 = np.random.randint(eval_per_gen-1, count_l1)
            x_tail[eval_per_gen-1,:] = generated_level_1[mini_max_dist_ind[rand_num1],:]
            v_tail[eval_per_gen-1,:] = generated_level_1_v[mini_max_dist_ind[rand_num1],:]
            #print(selected_function(generated_level_1[mini_max_dist_ind,:]))
            #dynamic increase the database
            count = 0
            for ii in range(eval_per_gen):
                flag = 0
                for jj in range(db_x.shape[0]):
                    if all(x_tail[ii] == db_x[jj]):
                        flag = 1
                        count = count + 1
                        break
                if flag == 0:
                    db_x = np.concatenate((db_x, x_tail[ii,:].reshape(1,D)), axis=0)
                    db_v = np.concatenate((db_v, v_tail[ii,:].reshape(1,D)), axis=0)
                    db_y = np.concatenate((db_y, selected_function(x_tail[ii,:].reshape(1,D))), axis=0)
            #print(db_x.shape, db_v.shape, db_y.shape)
            best_ever[times,itera] = np.amin(db_y)
            print(itera, min(db_y), count)
            #print(con_count[times, itera])
        print(times, best_ever[times,-1])

    end_time = time.time()
    print('time cost: ', end_time - start_time, 's')
    best_ever = best_ever.transpose()
    best_average = np.mean(best_ever, axis=1)
    best_std = np.std(best_ever, axis=1)
    print(best_average[-1], best_std[-1])
    best_2_save = np.concatenate((best_ever, best_average.reshape(max_gen,1)),axis=1)
    #np.savetxt('new_L1explotation_f%d_D%d_eval%d_1L%d_ite%d.txt' %(fi+1, D, eval_per_gen, count_l1, max_times), best_2_save)
    #np.savetxt('convergence_speed_f%d_D%d_eval%d_1L%d_ite%d.csv' %(fi+1, D, eval_per_gen, count_l1, max_times), np.mean(con_count, axis=0).transpose())
