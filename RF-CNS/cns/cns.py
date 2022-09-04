# Competitive neighborhood search

#### Authors: Lei Han, Handing Wang
#### Xidian University, China.
#### EMAIL: lhan_1@stu.xidian.edu.cn, hdwang@xidian.edu.cn
#### WEBSITE: https://sites.google.com/site/handingwanghomepage
#### DATE: January 2021
# ------------------------------------------------------------------------
# This code is part of the program that produces the results in the following paper:
#
# Lei Han, Handing Wang, A random forest assisted evolutionary algorithm using competitive neighborhood search for expensive constrained combinatorial optimization. Memetic Computing, accepted.
#
# You are free to use it for non-commercial purposes. However, we do not offer any forms of guarantee or warranty associated with the code. We would appreciate your acknowledgement.
# ------------------------------------------------------------------------
import random
import matplotlib.pyplot as plt
import copy
from sklearn.ensemble import RandomForestRegressor


class CNS:
    def __init__(self, max_iter, n_var, n_pop, eval_func, knapsack_size):
        self.max_iter = max_iter
        self.n_var = n_var
        self.n_pop = n_pop
        self.knapsack_size = knapsack_size
        self.eval_func = eval_func
        self.iter_val = [0 for i in range(max_iter)]
        self.iter_wgh = [0 for i in range(max_iter)]

        self.pop_train = []
        self.pop_train_val = []
        self.pop_train_wgh = []
        self.pop_train_vlt = []

        self.pop = []
        self.pop_val = []
        self.pop_wgh = []
        self.pop_vlt = []

        self.pop1 = []
        self.pop1_val = []
        self.pop1_wgh = []
        self.pop1_vlt = []

        self.pop2 = []
        self.pop2_val = []
        self.pop2_wgh = []
        self.pop2_vlt = []

        self.pop_offspring = []
        self.pop_offspring_val = []
        self.pop_offspring_wgh = []
        self.pop_offspring_vlt = []

        self.ind_best = []
        self.ind_best_val = 0
        self.ind_best_wgh = 0
        self.ind_best_vlt = 0

        self.ind_candidate = []
        self.ind_candidate_val = []
        self.ind_candidate_wgh = []
        self.ind_candidate_vlt = []

        self.rfr_wgh = RandomForestRegressor(n_estimators=30, random_state=0)
        self.rfr_val = RandomForestRegressor(n_estimators=30, random_state=0)

    def initialize(self):
        for i in range(self.n_pop):
            rand = random.random()
            if rand < 0.05:
                rand = 0.05
            elif rand > 0.95:
                rand = 0.95
            else:
                pass

            while True:
                temp = []
                for j in range(self.n_var):
                    r = random.random()
                    if r < rand:
                        temp_1 = 1
                    else:
                        temp_1 = 0
                    temp.append(temp_1)
                if temp not in self.pop_train:
                    break
            self.pop_train.append(temp)
        self.pop_train_val, self.pop_train_wgh, self.pop_train_vlt = self.eval_func(self.pop_train, self.knapsack_size)
        self.pop = copy.deepcopy(self.pop_train)
        self.pop_val = copy.deepcopy(self.pop_train_val)
        self.pop_wgh = copy.deepcopy(self.pop_train_wgh)
        self.pop_vlt = copy.deepcopy(self.pop_train_vlt)

        self.ind_best = self.pop_train[0]
        self.ind_best_val = self.pop_train_val[0]
        self.ind_best_wgh = self.pop_train_wgh[0]
        self.ind_best_vlt = self.pop_train_vlt[0]

        for i in range(self.n_pop - 1):
            if self.pop_train_vlt[i + 1] < self.ind_best_vlt:
                self.ind_best = self.pop_train[i + 1]
                self.ind_best_val = self.pop_train_val[i + 1]
                self.ind_best_wgh = self.pop_train_wgh[i + 1]
                self.ind_best_vlt = self.pop_train_vlt[i + 1]
            elif self.pop_train_vlt[i + 1] == self.ind_best_vlt:
                if self.pop_train_val[i + 1] > self.ind_best_val:
                    self.ind_best = self.pop_train[i + 1]
                    self.ind_best_val = self.pop_train_val[i + 1]
                    self.ind_best_wgh = self.pop_train_wgh[i + 1]
                    self.ind_best_vlt = self.pop_train_vlt[i + 1]

    def train_model(self):
        self.rfr_val.fit(self.pop_train, self.pop_train_val)
        self.rfr_wgh.fit(self.pop_train, self.pop_train_wgh)

    def model_predict(self, x):
        val = self.rfr_val.predict(x)
        wgh = self.rfr_wgh.predict(x)
        vlt = []
        for i in range(len(wgh)):
            if wgh[i] > self.knapsack_size:
                vlt.append(wgh[i] - self.knapsack_size)
            else:
                vlt.append(0)
        return val, wgh, vlt

    def neighborhood_search_1(self):
        self.pop1 = []
        self.pop1_val = []
        self.pop1_wgh = []
        self.pop1_vlt = []
        for i in range(self.n_pop):
            ind = copy.deepcopy(self.pop[i])
            index = random.randint(0, self.n_var - 1)
            if ind[index] == 0:
                ind[index] = 1
            else:
                ind[index] = 0
            self.pop1.append(ind)
        self.pop1_val, self.pop1_wgh, self.pop1_vlt = self.model_predict(self.pop1)

    def neighborhood_search_2(self):
        self.pop2 = []
        self.pop2_val = []
        self.pop2_wgh = []
        self.pop2_vlt = []
        for i in range(self.n_pop):
            ind = copy.deepcopy(self.pop[i])
            index = random.sample(range(self.n_var), 2)
            for j in range(2):
                if ind[index[j]] == 0:
                    ind[index[j]] = 1
                else:
                    ind[index[j]] = 0
            self.pop2.append(ind)
        self.pop2_val, self.pop2_wgh, self.pop2_vlt = self.model_predict(self.pop2)

    def neighborhood_competition(self):
        self.pop_offspring = []
        self.pop_offspring_val = []
        self.pop_offspring_wgh = []
        self.pop_offspring_vlt = []
        for i in range(self.n_pop):
            if self.pop1_vlt[i] < self.pop2_vlt[i]:
                self.pop_offspring.append(self.pop1[i])
                self.pop_offspring_val.append(self.pop1_val[i])
                self.pop_offspring_wgh.append(self.pop1_wgh[i])
                self.pop_offspring_vlt.append(self.pop1_vlt[i])
            elif self.pop1_vlt[i] == self.pop2_vlt[i]:
                if self.pop1_val[i] > self.pop2_val[i]:
                    self.pop_offspring.append(self.pop1[i])
                    self.pop_offspring_val.append(self.pop1_val[i])
                    self.pop_offspring_wgh.append(self.pop1_wgh[i])
                    self.pop_offspring_vlt.append(self.pop1_vlt[i])
                else:
                    self.pop_offspring.append(self.pop2[i])
                    self.pop_offspring_val.append(self.pop2_val[i])
                    self.pop_offspring_wgh.append(self.pop2_wgh[i])
                    self.pop_offspring_vlt.append(self.pop2_vlt[i])
            else:
                self.pop_offspring.append(self.pop2[i])
                self.pop_offspring_val.append(self.pop2_val[i])
                self.pop_offspring_wgh.append(self.pop2_wgh[i])
                self.pop_offspring_vlt.append(self.pop2_vlt[i])

    def selection(self):
        index = random.sample(range(self.n_pop), self.n_pop)
        for i in range(self.n_pop):
            if self.pop_vlt[i] > self.pop_offspring_vlt[index[i]]:
                self.pop[i] = self.pop_offspring[index[i]]
                self.pop_val[i] = self.pop_offspring_val[index[i]]
                self.pop_wgh[i] = self.pop_offspring_wgh[index[i]]
                self.pop_vlt[i] = self.pop_offspring_vlt[index[i]]
            elif self.pop_vlt[i] == self.pop_offspring_vlt[index[i]]:
                if self.pop_val[i] < self.pop_offspring_val[index[i]]:
                    self.pop[i] = self.pop_offspring[index[i]]
                    self.pop_val[i] = self.pop_offspring_val[index[i]]
                    self.pop_wgh[i] = self.pop_offspring_wgh[index[i]]
                    self.pop_vlt[i] = self.pop_offspring_vlt[index[i]]

    def update_model(self):
        pop_sort = self.feasibility_ranking()
        for i in range(len(pop_sort)):
            self.ind_candidate = pop_sort[i]
            if self.ind_candidate not in self.pop_train:
                break
        if self.ind_candidate in self.pop_train:
            for i in range(self.n_var):
                temp = copy.deepcopy(pop_sort[0])
                if temp[i] == 0:
                    temp[i] = 1
                else:
                    temp[i] = 0
                if temp not in self.pop_train:
                    self.ind_candidate = temp
                    break
        if self.ind_candidate in self.pop_train:
            for i in range(self.n_var - 1):
                temp_1 = copy.deepcopy(pop_sort[0])
                if temp_1[i] == 0:
                    temp_1[i] = 1
                else:
                    temp_1[i] = 0
                for j in range(self.n_var - i - 1):
                    temp_2 = copy.deepcopy(temp_1)
                    if temp_2[i + j + 1] == 0:
                        temp_2[i + j + 1] = 1
                    else:
                        temp_2[i + j + 1] = 0
                    if temp_2 not in self.pop_train:
                        self.ind_candidate = temp_2
                        break
                if self.ind_candidate not in self.pop_train:
                    break

        self.ind_candidate_val, self.ind_candidate_wgh, self.ind_candidate_vlt = self.eval_func([self.ind_candidate],
                                                                                                self.knapsack_size)

        self.pop_train.append(self.ind_candidate)
        self.pop_train_val.append(self.ind_candidate_val[0])
        self.pop_train_wgh.append(self.ind_candidate_wgh[0])
        self.pop_train_vlt.append(self.ind_candidate_vlt[0])
        self.rfr_val.fit(self.pop_train, self.pop_train_val)
        self.rfr_wgh.fit(self.pop_train, self.pop_train_wgh)

    def feasibility_ranking(self):
        pop_sort = copy.deepcopy(self.pop)
        pop_sort_val = copy.deepcopy(self.pop_val)
        pop_sort_vlt = copy.deepcopy(self.pop_vlt)
        for i in range(self.n_pop - 1):
            mark = 0
            for j in range(self.n_pop - i - 1):
                if pop_sort_vlt[j] == pop_sort_vlt[j + 1]:
                    if pop_sort_val[j] < pop_sort_val[j + 1]:
                        temp = pop_sort[j]
                        temp_val = pop_sort_val[j]
                        temp_vlt = pop_sort_vlt[j]

                        pop_sort[j] = pop_sort[j + 1]
                        pop_sort_val[j] = pop_sort_val[j + 1]
                        pop_sort_vlt[j] = pop_sort_vlt[j + 1]

                        pop_sort[j + 1] = temp
                        pop_sort_val[j + 1] = temp_val
                        pop_sort_vlt[j + 1] = temp_vlt

                        mark = 1
                elif pop_sort_vlt[j] > pop_sort_vlt[j + 1]:
                    temp = pop_sort[j]
                    temp_val = pop_sort_val[j]
                    temp_vlt = pop_sort_vlt[j]

                    pop_sort[j] = pop_sort[j + 1]
                    pop_sort_val[j] = pop_sort_val[j + 1]
                    pop_sort_vlt[j] = pop_sort_vlt[j + 1]

                    pop_sort[j + 1] = temp
                    pop_sort_val[j + 1] = temp_val
                    pop_sort_vlt[j + 1] = temp_vlt

                    mark = 1
            if mark == 0:
                break
        return pop_sort

    def update_best(self):
        if self.ind_candidate_vlt[0] < self.ind_best_vlt:
            self.ind_best = self.ind_candidate
            self.ind_best_val = self.ind_candidate_val[0]
            self.ind_best_wgh = self.ind_candidate_wgh[0]
            self.ind_best_vlt = self.ind_candidate_vlt[0]
        elif self.ind_candidate_vlt[0] == self.ind_best_vlt:
            if self.ind_candidate_val[0] > self.ind_best_val:
                self.ind_best = self.ind_candidate
                self.ind_best_val = self.ind_candidate_val[0]
                self.ind_best_wgh = self.ind_candidate_wgh[0]
                self.ind_best_vlt = self.ind_candidate_vlt[0]

    def show_result(self):
        plt.plot(range(self.max_iter), self.iter_wgh, 'blue', label='Weight')
        plt.plot(range(self.max_iter), self.iter_val, 'red', label='Value')
        plt.title('The convergence curve of the RF-CNS')
        plt.legend()
        plt.xlabel('Iterations')
        plt.ylabel('Values and weights')
        plt.show()

    def print_result(self):
        print('The value of the best individual:', self.ind_best_val)
        print('The weight of the best individual:', self.ind_best_wgh)
        print('The best individual:', self.ind_best)

    def run(self):
        self.initialize()
        self.train_model()
        for i in range(self.max_iter):
            self.neighborhood_search_1()
            self.neighborhood_search_2()
            self.neighborhood_competition()
            self.selection()
            self.update_model()
            self.update_best()
            self.iter_val[i] = self.ind_best_val
            self.iter_wgh[i] = self.ind_best_wgh
            print('Iteration,Value,Weight', (i, self.iter_val[i], self.iter_wgh[i]))
