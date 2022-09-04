# An example for quadratic knapsack problems

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
from problems import problem_quadratic_knapsack
from cns import CNS


def main():
    max_iter = 500  # Iterations
    var_num = 50  # Problem dimensions
    pop_size = 100  # Population size
    w = 200  # Knapsack size

    cns_new = CNS(max_iter, var_num, pop_size, problem_quadratic_knapsack, w)
    cns_new.run()
    cns_new.show_result()
    cns_new.print_result()


if __name__ == '__main__':
    main()
