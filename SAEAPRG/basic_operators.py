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
import random
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import precision_score, accuracy_score
import copy
import sys
from deap.tools import cxUniform
import pandas as pd


# get the classification performance of a individual
# each binary coded individual is a selection strategy
# data, instance data
# ind, binary coded selection strategy
def get_error_rate(data, ind, pro=None):
    size = data.shape[0]
    increment = int(size / 10)
    # 打乱原来数据的顺序
    new_data = shuffle(data)
    x = new_data[:, :-1]
    y = new_data[:, -1]
    sub_x = get_training_data(x, ind)
    error_rate = 0
    accs = 0
    for c in range(10):
        # 2、划分训练集和测试集，10折交叉验证
        if c == 0:
            # 前9份是训练集，最后1份是测试集
            train_x = sub_x[0:size-increment, :]
            test_x = sub_x[size-increment:size, :]
            train_y = y[0:size-increment]
            test_y = y[size-increment:size]
        elif c == 9:
            train_x = sub_x[increment:size, :]
            test_x = sub_x[0:increment, :]
            train_y = y[increment:size]
            test_y = y[0:increment]
        else:
            train_x = np.vstack((sub_x[0:(c+1)*increment, :], sub_x[(c+2)*increment:size, :]))
            train_y = np.array(list(y[0:(c + 1) * increment]) + list(y[(c + 2) * increment:size]))
            test_x = sub_x[(c+1)*increment : (c+2)*increment, :]
            test_y = y[(c+1)*increment : (c+2)*increment]

        # train_x, test_x, train_y, test_y = ms.train_test_split(sub_x, y, test_size=0.3)
        knn = KNeighborsClassifier(5)
        knn.fit(train_x, train_y)
        pred_y = knn.predict(test_x)
        # indic = imbalance_indices(pred_y, test_y, pro)
        indic = precision_score(test_y, pred_y, average='macro')
        error_rate = error_rate + (1 - indic)
        accs = accs + accuracy_score(test_y, pred_y)
    return error_rate / 10


# Shuffle the order of instance data
def shuffle(data):
    size = data.shape[0]
    idx = random.sample(range(size), size)
    new_data = []
    for id in idx:
        new_data.append(list(data[id,:]))
    return np.array(new_data)


# 根据ind的取值获取相应的训练数据
# get a vertical subset of instance data according to a selection strategy
# data, instance data
# ind, a feature selection strategy(binary coded)
def get_training_data(data, ind):
    sub = []
    for i in range(len(ind)):
        if ind[i] == 1:
            sub.append(data[:,i].tolist())
    sub = np.array(sub)
    return sub.T


# randomly select a horizontal subset of instance data
# data, the set of features
# labels, the set of labels
# size, the size of selected subset
def random_select(data, labels, size):
    if size >= len(data):
        return np.array(data), np.array(labels)
    else:
        selected_data = []
        selected_value = []
        idx = random.sample(range(0, len(data)), size)
        for id in idx:
            selected_data.append(data[id])
            selected_value.append(labels[id])
        return np.array(selected_data), np.array(selected_value)


# get a subset of instance data according to a random group
# data, the set of features
# group, a random group
def sort_out_instance_data(data, group):
    ps = []
    for p in data:
        ind = []
        for gr in group:
            ind.append(p[gr])
        ps.append(ind)
    return np.array(ps)

# 根据维度挑选出
def sort_out_ind(pop, group):
    group = np.sort(group)
    # new_pop = []
    for p in pop:
        m = len(group) - 1
        n = len(p) - 1
        while n >= 0:
            if n != group[m]:
                del p[n]
            else:
                if m > 0:
                    m = m - 1
            n = n - 1
        # new_pop.append(p)
    return pop


# 使用代理模型评价个体，注意：子问题的维度升为原问题的维度时，使用0填充
# evaluate the population using a surrogate model
# pop, a population
# model, a trained surrogate model
# c, constraint boundary
def evaluate(pop, model, c):
    temp_pop = ind_to_list(pop)
    # 使用模型预测
    val = model.predict(temp_pop)
    # 表示不可行解的个数
    # count = 0
    for i in range(len(pop)):
        pop[i].dc = np.sum(pop[i]) - c
        if np.sum(pop[i]) - c > 0:
            pop[i].fitness.values = [val[i] + (5 * (np.sum(pop[i]) - c))]
        else:
            pop[i].fitness.values = [val[i]]

# 将元素为individual的list转换为array数组，便于代理模型预测
def ind_to_list(pop):
    c = np.zeros((len(pop), len(pop[0])))
    for i in range(len(pop)):
        for j in range(len(pop[i])):
            c[i,j] = pop[i][j]
    return c

# crossover and mutation
def reproduction(parent, cidx, midx):
    middle = int(len(parent) / 2)
    offspring = []
    for i in range(middle):
        # 两个父代个体
        p1 = copy.deepcopy(parent[i])
        p2 = copy.deepcopy(parent[middle+i])
        idx = np.random.uniform(0, 1)
        if idx < cidx:
            o1, o2 = cxUniform(p1, p2, 0.5)
            offspring.append(copy.deepcopy(o1))
            offspring.append(copy.deepcopy(o2))
        else:
            offspring.append(p1)
            offspring.append(p2)
    for j in range(len(offspring)):
        for k in range(len(offspring[j])):
            idx = np.random.uniform(0, 1)
            if idx < midx:
                if offspring[j][k] == 1:
                    offspring[j][k] = 0
                else:
                    offspring[j][k] = 1
    return offspring


# 将dim维的决策变量随机分为g组
# 每组维度之间没有耦合
# random grouping
# dim, the size of all the features
# length, the size of each group
def random_grouping(dim, length):
    groups = []
    # 决策变量的位置
    variables = np.arange(0, dim).tolist()
    while variables:
        if len(variables) > length:
            dims = random.sample(variables, length)
        else:
            dims = random.sample(variables, len(variables))
        groups.append(list(np.sort(dims)))
        dims.sort(reverse=True)
        for k in dims:
            variables.remove(k)
    return groups

# There is no coupling between the groups
# merge the subpopulation to the population
# pops, all the subpopulations
# groups, random groups
# dim, the size of all the features
def merge_ind(pops, groups, dim):
    # 合并后的种群
    pop = copy.deepcopy(pops[0])
    # 扩展到原问题的维度
    for i in range(len(pop)):
        count = dim - len(pop[i])
        for j in range(count):
            pop[i].append(0)
    # 逐个合并
    for g in range(len(groups)):
        group = groups[g]
        p = pops[g]
        for k in range(len(p)):
            for l in range(len(group)):
                pop[k][group[l]] = p[k][l]
    return pop


# 初始化个体
def int_ind(size):
    ind = []
    for i in range(size):
        r = np.random.random()
        if r <= 0.5:
            ind.append(0)
        else:
            ind.append(1)
    return ind

# 由于维度太高，全局模型不可用，使用局部模型集成的思想评价个体
# using ensemble model to re-evaluate saved individuals
def evaluate_by_ensembles(pop, models, groups, c):

    if len(models) != len(groups):
        print('Error Input!')
        return -1
    else:
        fitness = np.zeros(len(pop))
        for i in range(len(models)):
            model = models[i]
            group = groups[i]
            tmp = sort_out_instance_data(pop, group)
            val = model.predict(np.array(tmp))
            fitness = fitness + val
        for i in range(len(pop)):
            if np.sum(pop[i]) - c > 0:
                fitness[i] = sys.maxsize
        return np.argmin(fitness)


# 使用分类器评价每个个体
def get_fitness_by_classifier(pop, name):
    file_name = 'data/' + name + '.csv'
    data = pd.read_csv(file_name, sep=",", header=None)
    data = data.values
    fitness = []
    for d in pop:
        fitness.append(get_error_rate(data, d, name))
    return fitness
