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

import numpy as np
import copy
import sys
sys.path.append('..')
from rbfn import RBFN as RBFN
import csv
import multiprocessing
from datetime import datetime
import dill
import basic_operators as bo
from deap import creator, base, tools

class SAEAPRG:

    def __init__(self, data, fitness, pro):
        self.data = np.array(data)
        self.fitness = np.array(fitness)
        self.pro = pro
        self.dim = self.data.shape[1]
        self.c = int(self.dim*0.15)

    # optimize each sub-problem
    def opt_subproblems(self, sub_group, sub_c, pop, q1=None, q2=None):
        # 定义一些基本参数
        popSize = 200
        # 交叉概率
        cidx = 0.9
        # 当前迭代次数
        gen = 0
        # 最大迭代次数
        max_gen = 15
        # 训练代理模型
        model = RBFN(len(sub_group), int(len(sub_group) / 2))
        selected_data, selected_value = bo.random_select(self.data, self.fitness, 25 * len(sub_group))
        # 分拣决策变量
        train_data = bo.sort_out_instance_data(selected_data, sub_group)
        # 训练代理模型
        model.fit(np.array(train_data), selected_value)
        # 初始迭代的初始种群，种群的大小为原问题的大小
        sub_pop = bo.sort_out_ind(copy.deepcopy(pop), sub_group)
        # 评价种群
        bo.evaluate(sub_pop, model, sub_c)
        # 存储模型
        if q1 != None and q2 != None:
            q1.put(dill.dumps(model))
            q2.put(sub_group)
        while gen < max_gen:
            offspring = bo.reproduction(sub_pop, cidx, 1 / len(sub_group))
            bo.evaluate(offspring, model, sub_c)
            combinedPop = sub_pop + offspring
            sub_pop = tools.selBest(combinedPop, k=popSize)
            gen = gen + 1

        for ind in sub_pop:
            if np.sum(ind) <= sub_c:
                optima_fitness = ind.fitness.values[0]
                break
        if q1 != None and q2 != None:
            return sub_pop, sub_group, optima_fitness
        else:
            return sub_pop, optima_fitness, model

    # dims存储了需要优化的维度
    # 优化dims中存储的相关维度
    # 其余维度使用local_best的相关维度的值填充
    # data, target表示在决策空间中采样的真实数据，用来训练代理模型
    def ga(self, pop):
        # 每个分组的大小
        group_size = 200
        # 对决策变量随机分成group_num组，每个分组之间是互斥的
        groups = bo.random_grouping(self.dim, group_size)
        # 每个子问题的约束
        sub_contains = []
        # 根据子问题的维度均分约束
        for group_i in range(len(groups)):
            tmp = int(len(groups[group_i]) * self.c / self.dim)
            if tmp > 1:
                sub_contains.append(int(len(groups[group_i]) * self.c / self.dim))
            else:
                sub_contains.append(1)
        # 处理所有子问题的约束之和小于原问题约束的情况
        if np.sum(sub_contains) < self.c:
            diff = self.c - np.sum(sub_contains)
            while diff > 0:
                id = np.random.randint(0, len(sub_contains))
                sub_contains[id] = sub_contains[id] + 1
                diff = diff - 1
        # 多进程优化各个子问题
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
        model_que = multiprocessing.Manager().Queue(len(groups))
        group_que = multiprocessing.Manager().Queue(len(groups))
        results = []
        for i in range(len(groups)):
            results.append(pool.apply_async(self.opt_subproblems, (groups[i], sub_contains[i], pop, model_que, group_que)))
        pool.close()
        pool.join()
        pops = []
        tmp_groups = []
        optima_fitnesses = []
        for res in results:
            sub_pop, group, optima_fitness = res.get()
            pops.append(copy.deepcopy(sub_pop))
            tmp_groups.append(copy.deepcopy(group))
            optima_fitnesses.append(optima_fitness)
        full_pop = bo.merge_ind(pops, tmp_groups, self.dim)
        r = np.random.randint(0, model_que.qsize())
        optima_idx = np.argmin(optima_fitnesses)
        for ri in range(model_que.qsize()):
            if ri == r:
                random_model = dill.loads(model_que.get())
                random_group = group_que.get()
                break
            else:
                model_que.get()
                group_que.get()
        return full_pop, random_model, random_group, tmp_groups[optima_idx], pops[optima_idx][0], optima_fitness


    def optimize(self):
        # 定义问题
        creator.create('FitnessMin', base.Fitness, weights=(-1.0,))
        # dc，违反约束的程度；g，当前的分组；
        creator.create('Individual', list, fitness=creator.FitnessMin, dc=None, g = None)
        # 定义参数，toolbox存储基本的参数
        toolbox = base.Toolbox()
        toolbox.register('attr_int', bo.int_ind, size=self.dim)
        toolbox.register('individual', tools.initIterate, creator.Individual, toolbox.attr_int)
        toolbox.register('population', tools.initRepeat, list, toolbox.individual)
        pop = toolbox.population(len(self.data))

        max_fit= 0
        max_cycles = 30

        for i in range(len(pop)):
            for j in range(len(pop[i])):
                pop[i][j] = self.data[i, j]
            pop[i].dc = np.sum(pop[i]) - self.c
            if np.sum(pop[i]) - self.c > 0:
                pop[i].fitness.values = [self.fitness[i] + (5 * pop[i].dc)]
            else:
                pop[i].fitness.values = [self.fitness[i]]

            if pop[i].fitness.values[0] > max_fit:
                 max_fit = pop[i].fitness.values[0]

        pop = tools.selBest(pop, 200)
        max_fit = 10
        models = []
        groups = []
        archive = []
        af_dict = {}

        for i in range(max_cycles):
            pop, random_model, random_group, optima_group, optima_ind, optima_fit = self.ga(copy.deepcopy(pop))
            # 用于构成集成模型
            models.append(random_model)
            groups.append(random_group)

            for ind in pop:
                if np.sum(ind) <= self.c:
                    archive.append(copy.deepcopy(ind))
                    break

            for gi in range(len(optima_group)):
                if optima_ind[gi] == 1:
                    g = optima_group[gi]
                    if g not in af_dict.keys():
                        af_dict[g] = (max_fit - optima_fit) * (1 + i / max_cycles)
                    else:
                        af_dict[g] = af_dict[g] + (max_fit - optima_fit) * (1 + i / max_cycles)

        all = pop + archive
        optima = bo.evaluate_by_ensembles(all, models, groups, self.c)

        p1_i = 0
        while p1_i < len(pop):
            if np.sum(pop[p1_i]) <= self.c:
                p1 = copy.deepcopy(pop[p1_i])
                break
            else:
                p1_i = p1_i + 1

        p2 = copy.deepcopy(all[optima])
        c1 = 0
        c2 = 0
        for pi in range(len(p1)):
            if p1[pi] == 1 and pi in af_dict.keys():
                c1 = c1 + af_dict[pi]
            if p2[pi] == 1 and pi in af_dict.keys():
                c2 = c2 + af_dict[pi]

        if c1 < c2:
            o3 = bo.get_fitness_by_classifier([p2], self.pro)
        else:
            o3 = bo.get_fitness_by_classifier([p1], self.pro)

        # 写入结果
        save_file3 = 'results/saeaprg.csv'
        with open(save_file3, 'a+') as csvFile3:
            writer = csv.writer(csvFile3)
            writer.writerow([self.pro, o3[0]])
