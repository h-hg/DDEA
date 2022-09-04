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

import sys
sys.path.append('..')
from sampling import sample
from saeaprg import SAEAPRG

problems = ['madelon', 'isolet', 'CNAE-9', 'ORL', 'COIL20', 'warpPIE10P', 'lung', 'lymphoma', 'relathe', 'DLBCL', 'leukemia', 'Gutenberg', 'arcene']
#
#
for pro in problems:
    for i in range(1):
        A = sample(pro)
        data = A[:, :-1]
        fitness = A[:, -1]
        algorithm = SAEAPRG(data, fitness, pro)
        algorithm.optimize()
print('Terminated!')
