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
from sklearn.cluster import KMeans
from numba import jit

@jit
# 计算a和b之间的汉明距离
def Hamming_distance(a, b):
    dis = np.nonzero(a - b)[0]
    return dis.size

@jit
def Gaussianfun(center, data_point, sigma): # 高斯函数
    return np.exp(-0.5 * np.power(Hamming_distance(center, data_point) / sigma, 2))

@jit
def discretize_centers(centers):
    for i in range(centers.shape[0]):
        for j in range(centers.shape[1]):
            if centers[i, j] <= 0.5:
                centers[i,j] = 0
            else:
                centers[i,j] = 1
    return centers

@jit
def calsigma(hidden_shape, centers):
    max = 0.0
    num = 0
    total = 0.0
    for i in range(hidden_shape-1):
        for j in range(i+1, hidden_shape):
            # 汉明距离
            dis = Hamming_distance(centers[i], centers[j])
            total = total + dis
            num += 1
            if dis >max:
                max = dis
    sigma = 2*total/(num+0.001)
    return sigma

@jit
def calculate_interpolation_matrix(X, hidden_shape, centers, sigma):
    G = np.zeros((X.shape[0], hidden_shape))
    for data_point_arg, data_point in enumerate(X):
        for center_arg, center in enumerate(centers):
            G[data_point_arg,center_arg] = Gaussianfun(center, data_point, sigma)
    return G


class RBFN(object):

    def __init__(self, input_shape, hidden_shape):

        self.input_shape = input_shape
        self.hidden_shape = hidden_shape
        self.centers = None
        self.sigm = None
        self.weights = None
        self.bias = None

    def fit(self, X, Y):
        km = KMeans(n_clusters=self.hidden_shape).fit(X)
        self.centers = km.cluster_centers_
        self.centers = discretize_centers(self.centers)
        self.sigma = calsigma(self.hidden_shape, self.centers)
        G = calculate_interpolation_matrix(X, self.hidden_shape, self.centers, self.sigma)
        temp = np.ones((len(X)))
        temp = np.column_stack((G, temp))
        temp = np.dot(np.linalg.pinv(temp), Y)
        self.weights = temp[:self.hidden_shape]
        self.bias = temp[self.hidden_shape]


    def predict(self, X):
        G = calculate_interpolation_matrix(X, self.hidden_shape, self.centers, self.sigma)
        predictions = np.dot(G, self.weights) + self.bias
        return predictions
