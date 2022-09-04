import sys
sys.path.append('..')
import pandas as pd
import numpy as np
import copy
from deap import creator, base, tools
from operator import attrgetter
from benchmark.problems import evaluate_true
import hvwfg
import csv
from sklearn.ensemble import RandomForestRegressor
from itertools import chain
import algorithms.utils as us

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

# evaluate a population using surrogate models
def evaluate(pop, p_namme, rf1, rf2, rf3):
    if len(pop) == 1:
        pop = [pop]
    if p_namme[0:3] == 'zdt' or p_namme[0:3] == 'wfg' or p_namme == 'uf1' or p_namme == 'uf2' or p_namme == 'uf4' or p_namme == 'uf6':
        f1 = rf1.predict(pop)
        f2 = rf2.predict(pop)
        for i in range(len(pop)):
            pop[i].fitness.values = [f1[i], f2[i]]
    else:
        f1 = rf1.predict(pop)
        f2 = rf2.predict(pop)
        f3 = rf3.predict(pop)
        for i in range(len(pop)):
            pop[i].fitness.values = [f1[i], f2[i], f3[i]]

# generate a individual using a random method
def my_randomint(x_1, x_other, size):
    ind = []
    x_1_index = np.random.randint(0, x_1.size)
    ind.append(np.round(x_1[x_1_index], 2))
    for i in range(1, size):
        x_other_index = np.random.randint(0, x_other.size)
        ind.append(np.round(x_other[x_other_index], 2))
    return ind

# point-based variation operator
def point_reproduction(parent, cidx, midx, x_1, x_other):
    offspring = []
    middle = int(len(parent)/2)
    for i in range(middle):
        # 两个父代个体
        p1 = copy.deepcopy(parent[i])
        p2 = copy.deepcopy(parent[middle+i])
        # 交叉
        if np.random.random() < cidx:
            idx = np.random.uniform(0, 1, len(p1))
            for t in range(len(idx)):
                if idx[t] < 0.3:
                    tmp_value = p1[t]
                    p1[t] = np.round(p2[t], 2)
                    p2[t] = np.round(tmp_value, 2)
        # 变异
        for j in range(len(p1)):
            if np.random.uniform(0, 1) < midx:
                if j == 0:
                    x_index = np.random.randint(0, x_1.size)
                    p1[j] = np.round(x_1[x_index], 2)
                else:
                    x_index = np.random.randint(0, x_other.size)
                    p1[j] = np.round(x_other[x_index], 2)
        for j in range(len(p2)):
            if np.random.uniform(0, 1) < midx:
                if j == 0:
                    x_index = np.random.randint(0, x_1.size)
                    p2[j] = np.round(x_1[x_index], 2)
                else:
                    x_index = np.random.randint(0, x_other.size)
                    p2[j] = np.round(x_other[x_index], 2)
        offspring.append(p1)
        offspring.append(p2)
    return offspring

# pop, set of nondominated solutions
# pro, optimization problem
# x_1, discrete values of the first decision variable
# x_other, discrete values of the other decision variables
# rf1, rf2, rf3 surrogate models
def semi_vns(pop, pro, x_1, x_other,  rf1, rf2, rf3):
    offspring = []
    for i in range(len(pop)):
        s = copy.deepcopy(pop[i])
        # the number of neighborhoods
        k = 1
        tmp_off = []
        while(k<=2):
            # generate the k^th neighborhood of s
            if k == 1:
                ns = us.neighborhood_one(s, x_1, x_other, 5)
            elif k == 2:
                ns = us.neighborhood_two(s, x_1, x_other, 10)
            tmp_off = tmp_off + ns
            # Shaking. Generate a point s' at random from the k^th neighborhood of s
            s_t = us.random_select_from_one_set(ns)
            nst1 = us.neighborhood_one(s_t, x_1, x_other, 5)
            nst2 = us.neighborhood_two(s_t, x_1, x_other, 10)
            tmp_off = tmp_off + nst1 + nst2
            k = k + 1
        evaluate(tmp_off, pro, rf1, rf2, rf3)
        del_dominated_solutions(tmp_off, pop[i])
        offspring = offspring + copy.deepcopy(tmp_off)
    return offspring


def del_dominated_solutions(pop, ind):
    i = len(pop) -1
    while i >= 0:
        if is_dominated(pop[i], ind) == False:
            del pop[i]
        i = i - 1
    return pop


def is_dominated(ind1, ind2):
    wvalues1 = ind1.fitness.wvalues
    wvalues2 = ind2.fitness.wvalues
    not_equal = False
    for self_wvalue, other_wvalue in zip(wvalues1, wvalues2):
        if self_wvalue > other_wvalue:
            return False
        elif self_wvalue < other_wvalue:
            not_equal = True
    return not_equal

# the method used to calculate HV according to nondominated solutions
def get_hv(ps, ref, obj_cout):
    # 计算HV
    if obj_cout == 2:
        pf = np.zeros((len(ps), 2))
        for p_i in range(len(ps)):
            pf[p_i, 0] = ps[p_i].fitness.values[0]
            pf[p_i, 1] = ps[p_i].fitness.values[1]
    elif obj_cout == 3:
        pf = np.zeros((len(ps), 3))
        for p_i in range(len(ps)):
            pf[p_i, 0] = ps[p_i].fitness.values[0]
            pf[p_i, 1] = ps[p_i].fitness.values[1]
            pf[p_i, 2] = ps[p_i].fitness.values[2]
    hv = hvwfg.wfg(pf, ref)
    return hv

def is_equal(p, q):
    flag = True
    if len(p) != len(q):
        flag = False
    else:
        i = 0
        while i < len(p):
            if p[i] != q[i]:
                flag = False
                break
            else:
                i = i + 1
    return flag

def get_equal_idx(pop, p):
    idx = -1
    for i in range(len(pop)):
        if is_equal(p, pop[i]):
           idx = i
           break
    return idx


def sel_by_non_cd(fronts, k):
    chosen = list(chain(*fronts[:-1]))
    k = k - len(chosen)
    if k > 0:
        sorted_front = sorted(fronts[-1], key=attrgetter("fitness.crowding_dist"), reverse=True)
        chosen.extend(sorted_front[:k])
    return chosen

# optimization process of EAHVFA
def eahvfa(pro, NDim):
    # define the optimization problem
    if pro[0:3] == 'zdt' or pro[0:3] == 'wfg' or (pro[0:3] in ['uf1', 'uf2', 'uf4', 'uf6']):
        creator.create('MultiObjMin', base.Fitness, weights=(-1.0, -1.0))
        obj_cout = 2
    else:
        creator.create('MultiObjMin', base.Fitness, weights=(-1.0, -1.0, -1.0))
        obj_cout = 3
    creator.create('Individual', list, fitness=creator.MultiObjMin)

    # 定义参数，toolbox存储基本的参数
    toolbox = base.Toolbox()
    # population size
    popSize = 100
    # sample size
    sample_size = 200
    toolbox.cidx = 0.9
    toolbox.midx =  1 / NDim
    # the interval of x is set to 0.02.
    x_1 = np.arange(0, 1.01, 0.02)
    if pro == 'zdt4':
        x_other = np.arange(-5, 5.01, 0.02)
    else:
        x_other = np.arange(0, 1.01, 0.02)
    x_1 = np.round(x_1, 2)
    x_other = np.round(x_other, 2)

    toolbox.register('attr_int', my_randomint, x_1, x_other, size=NDim)
    toolbox.register('individual', tools.initIterate, creator.Individual, toolbox.attr_int)
    toolbox.register('population', tools.initRepeat, list, toolbox.individual)

    pre = 0
    # the difference of HVs
    d_value = 0
    count = 0
    # 阈值
    threshold = 0
    # the lowest fidelity level
    fidelity = 1
    gen = 1
    # random sampling
    x, y1, y2, y3 = us.sample(pro, sample_size, NDim, fidelity)

    # train surrogate models with fidelity=1
    rf1 = RandomForestRegressor(n_estimators=1000)
    rf2 = RandomForestRegressor(n_estimators=1000)
    rf3 = []
    rf1.fit(x, y1)
    rf2.fit(x, y2)
    if pro[0:3] == 'dtl':
        rf3 = RandomForestRegressor(n_estimators=1000)
        rf3.fit(x, y3)

    # initialize population
    pop = toolbox.population(sample_size)
    for pi in range(len(pop)):
        for pj in range(len(pop[pi])):
            pop[pi][pj] = x[pi][pj]
            pop[pi][pj] = x[pi][pj]
        if obj_cout == 2:
            pop[pi].fitness.values = [y1[pi], y2[pi]]
        elif obj_cout == 3:
            pop[pi].fitness.values = [y1[pi], y2[pi], y3[pi]]

    # the reference point used to calculate HV
    if obj_cout == 2:
        ref = np.array([np.max(y1) * 5, np.max(y2) * 5])
    elif obj_cout == 3:
        ref = np.array([np.max(y1) * 5, np.max(y2) * 5, np.max(y3) * 5])
    # the number of Pareto solutions
    f_0 = 0


    while fidelity <= 10:

        if gen == 1 and fidelity == 1:
            # All the individuals use point-based variation operator to generate new solutions
            matting_pool = tools.selTournament(pop, k=popSize, tournsize=2)
            offspring = point_reproduction(matting_pool, toolbox.cidx, toolbox.midx, x_1, x_other)
            evaluate(offspring, pro, rf1, rf2, rf3)
        else:
            # the nondominated solutions using Semi-VNS to find promising solutions
            off1 = semi_vns(pop[0: f_0], pro,  x_1, x_other, rf1, rf2, rf3)
            if len(off1) == 0:
                off1 = point_reproduction(pop[0: f_0], toolbox.cidx, toolbox.midx, x_1, x_other)
            if f_0 < len(pop):
                matting_pool = tools.selTournament(pop[f_0: len(pop)], k=len(pop)-f_0, tournsize=2)
                off2 = point_reproduction(matting_pool, toolbox.cidx, toolbox.midx, x_1, x_other)
                if len(off2) != 0:
                    evaluate(off2, pro, rf1, rf2, rf3)
            else:
                off2 = []
            offspring = off1 + off2

        combinedPop = pop + offspring
        fronts = tools.emo.sortNondominated(combinedPop, k=popSize)
        ps = fronts[0]
        f_0 = len(ps)
        hv = get_hv(ps, ref, obj_cout)
        for front in fronts:
            tools.emo.assignCrowdingDist(front)

        d_value = hv - pre
        if d_value < 0:
            d_value = 0
        pre = hv
        if gen == 2:
            threshold = d_value / 10
            if threshold == 0:
                threshold = 1.0e-5

        if d_value < threshold or d_value <= 1.0e-5:
            count = count + 1
        else:
            count = 0

        if (count == 5 and fidelity < 10) or ( gen >= 50 and fidelity < 10):
            # update fidelity level
            fidelity = fidelity + 1
            gen = 1
            if fidelity <= 10:
                # Sampling
                po_size = int(fidelity * 20)
                random_size = sample_size - po_size
                s1 = sel_by_non_cd(fronts, k=po_size)
                if random_size > 0:
                    s2 = toolbox.population(random_size)
                else:
                    s2 = []
                samples = s1 + s2
                # evaluate sampled solutions using real fitness function
                # fidelity, updated fidelity
                evaluate_true(samples, pro, fidelity)
                # retrain surrogate models
                x_new = []
                y1_new = []
                y2_new = []
                y3_new = []

                for i in range(len(samples)):
                    x_new.append(list(samples[i]))
                    y1_new.append(samples[i].fitness.values[0])
                    y2_new.append(samples[i].fitness.values[1])
                    if obj_cout == 3:
                        y3_new.append(samples[i].fitness.values[2])

                rf1 = RandomForestRegressor(n_estimators=1000)
                rf2 = RandomForestRegressor(n_estimators=1000)
                if obj_cout == 3:
                    rf3 = RandomForestRegressor(n_estimators=1000)
                rf1.fit(x_new, y1_new)
                rf2.fit(x_new, y2_new)
                if obj_cout == 3:
                    rf3.fit(x_new, y3_new)

                combinedPop = copy.deepcopy(samples)

                fronts = tools.emo.sortNondominated(combinedPop, k=popSize)
                ps = fronts[0]
                for front in fronts:
                    tools.emo.assignCrowdingDist(front)

                # update reference points
                if obj_cout == 2:
                    ref = np.array([np.max(y1_new) * 5, np.max(y2_new) * 5])
                elif obj_cout == 3:
                    ref = np.array([np.max(y1_new) * 5, np.max(y2_new) * 5, np.max(y3_new) * 5])
                hv = get_hv(ps, ref, obj_cout)
                d_value = hv
                pre = hv
                count = 0
        gen = gen + 1


        if (count == 5 and fidelity == 10) or (gen >= 30 and fidelity == 10):
            fidelity = fidelity + 1

        pop = sel_by_non_cd(fronts, k=popSize)



    final_ps = tools.emo.sortNondominated(pop, k=len(pop), first_front_only=True)[0]
    evaluate_true(final_ps, pro, t=30)
    true_fronts = tools.emo.sortNondominated(final_ps, k=len(pop), first_front_only=True)
    true_ps = true_fronts[0]

    f1 = []
    f2 = []
    f3 = []
    if obj_cout == 2:
        p = np.zeros((len(true_ps), 2))
        for i in range(len(true_ps)):
            f1.append(true_ps[i].fitness.values[0])
            f2.append(true_ps[i].fitness.values[1])
            p[i, 0] = f1[i]
            p[i, 1] = f2[i]
    else:
        p = np.zeros((len(true_ps), 3))
        for i in range(len(true_ps)):
            f1.append(true_ps[i].fitness.values[0])
            f2.append(true_ps[i].fitness.values[1])
            f3.append(true_ps[i].fitness.values[2])
            p[i, 0] = f1[i]
            p[i, 1] = f2[i]
            p[i, 2] = f3[i]

    ref_pf = pd.read_csv('../pf/pf_'+pro+'.csv', sep=',', header=None).values

    # igd = meth.IGD(ref_pf, p)
    igd_plus = us.IGD_Plus(ref_pf, p)

    title = str('../results/eahvfa_igd_' + pro + '_' + str(NDim)+'.csv')
    # title2 = str('../results/eahvfa_ps_' + pro + '_' + str(NDim)+'_addi.csv')


    with open(title, 'a+') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerow([igd_plus])


