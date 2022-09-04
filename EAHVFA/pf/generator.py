import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv


def fast_nondominated_sort(fits):
    """快速非支配排序算法
    Params:
        fits: 适应度， narray数组
    Return:
        ranks: 每个个体所对应的等级，一维数组
    """
    nPop = fits.shape[0]
    nF = fits.shape[1]  # 目标函数的个数
    ranks = np.zeros(nPop, dtype=np.int32)
    nPs = np.zeros(nPop)  # 每个个体p被支配解的个数
    sPs = []  # 每个个体支配的解的集合，把索引放进去
    for i in range(nPop):
        iSet = []  # 解i的支配解集
        for j in range(nPop):
            if i == j:
                continue
            isDom1 = fits[i] <= fits[j]
            isDom2 = fits[i] < fits[j]
            # 是否支配该解-> i支配j
            if sum(isDom1) == nF and sum(isDom2) >= 1:
                iSet.append(j)
            # 是否被支配-> i被j支配
            if sum(~isDom2) == nF and sum(~isDom1) >= 1:
                nPs[i] += 1
        sPs.append(iSet)  # 添加i支配的解的索引
    r = 0  # 当前等级为 0， 等级越低越好
    indices = np.arange(nPop)
    while sum(nPs==0) != 0:
        rIdices = indices[nPs==0]  # 当前被支配数为0的索引
        ranks[rIdices] = r
        for rIdx in rIdices:
            iSet = sPs[rIdx]
            nPs[iSet] -= 1
        nPs[rIdices] = -1  # 当前等级的被支配数设置为负数
        r += 1
    return ranks


def zdt1():
    x1 = np.arange(0, 1.02, 0.02)
    f1 = np.around(x1, 2)
    f2 = 1 - np.sqrt(f1)
    plt.scatter(f1, f2)
    plt.show()

    save_file = 'pf_zdt1.csv'
    with open(save_file, 'a+') as csvFile:
        writer = csv.writer(csvFile)
        # 写入多行用writerows
        for i in range(len(f1)):
            writer.writerow([f1[i], np.round(f2[i], 4)])


def zdt2():
    x1 = np.arange(0, 1.02, 0.02)
    f1 = np.around(x1, 2)
    f2 = 1 - np.power(f1, 2)
    plt.scatter(f1, f2)
    plt.show()

    save_file = 'pf_zdt2.csv'
    with open(save_file, 'a+') as csvFile:
        writer = csv.writer(csvFile)
        # 写入多行用writerows
        for i in range(len(f1)):
            writer.writerow([f1[i], np.round(f2[i], 4)])



def zdt3():
    x1 = np.arange(0, 1.02, 0.02)
    f1 = np.around(x1, 2)
    f2 = 1 - np.sqrt(f1) - f1 * np.sin(10*np.pi*f1)
    f =  []
    for i in range(len(f1)):
        f.append([f1[i], np.round(f2[i], 4)])
    ranks = fast_nondominated_sort(np.array(f))
    pf = []
    for j in range(len(ranks)):
        if ranks[j] == 0:
            pf.append(f[j])
    plt.scatter(np.array(pf)[:,0], np.array(pf)[:,1])
    plt.show()

    save_file = 'pf_zdt3.csv'
    with open(save_file, 'a+') as csvFile:
        writer = csv.writer(csvFile)
        # 写入多行用writerows
        for i in range(len(pf)):
            writer.writerow(pf[i])

def zdt4():
    x1 = np.arange(0, 1.02, 0.02)
    f1 = np.around(x1, 2)
    f2 = 1 - np.sqrt(f1)
    plt.scatter(f1, f2)
    plt.show()

    save_file = 'pf_zdt4.csv'
    with open(save_file, 'a+') as csvFile:
        writer = csv.writer(csvFile)
        # 写入多行用writerows
        for i in range(len(f1)):
            writer.writerow([f1[i], np.round(f2[i], 4)])

def zdt6():
    x1 = np.arange(0, 1.02, 0.02)
    x1 = np.around(x1, 2)
    t1 = np.sin(6*np.pi*x1)
    t2 = np.power(t1, 6)
    f1 = 1 - np.exp(-4*x1)*t2
    f2 = 1 - np.power(f1, 2)
    plt.scatter(f1, f2)
    plt.show()

    save_file = 'pf_zdt6.csv'
    with open(save_file, 'a+') as csvFile:
        writer = csv.writer(csvFile)
        # 写入多行用writerows
        for i in range(len(f1)):
            writer.writerow([f1[i], np.round(f2[i], 4)])

def dtlz1():
    x1 = np.arange(0, 1.02, 0.02)
    x2 = np.arange(0, 1.02, 0.02)
    # g_x = 100 * (8 + 8 * (0.25 - np.cos(-10*np.pi)))
    g_x = 0
    pf = []
    for i in range(len(x1)):
        t1 = x1[i]
        for j in range(len(x2)):
            t2 = x2[j]
            f1 = 0.5 * t1 * t2 * (1 + g_x)
            f2 = 0.5 * t1 * (1 - t2) * (1 + g_x)
            f3 = 0.5 * (1 - t1) * (1 + g_x)
            s = f1 + f2 + f3
            if np.round(s, 2) == 0.5:
                pf.append([f1, f2, f3])
    pf = np.array(pf)
    pf = np.round(pf, 4)
    pf = np.unique(pf, axis=0)
    print(pf.shape)
    save_file = 'pf_dtlz1.csv'
    with open(save_file, 'a+') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(pf)


def dtlz2():
    x1 = np.arange(0, 1.02, 0.02)
    x2 = np.arange(0, 1.02, 0.02)
    # g_x = 100 * (8 + 8 * (0.25 - np.cos(-10*np.pi)))
    g_x = 0
    pf = []
    for i in range(len(x1)):
        t1 = x1[i]
        for j in range(len(x2)):
            t2 = x2[j]
            f1 = np.cos(0.5 * np.pi * t1) * np.cos(0.5 * np.pi * t2)
            f2 = np.cos(0.5 * np.pi * t1) * np.sin(0.5 * np.pi * t2)
            f3 = np.sin(0.5 * np.pi * t1)
            s = f1**2 + f2**2 + f3**2
            if np.round(s, 2) == 1:
                pf.append([f1, f2, f3])
    pf = np.array(pf)
    pf = np.round(pf, 4)
    pf = np.unique(pf, axis=0)
    print(pf.shape)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(pf[:,0], pf[:,1], pf[:,2])
    plt.show()
    save_file = 'pf_dtlz2.csv'
    with open(save_file, 'a+') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(pf)

def dtlz3():
    x1 = np.arange(0, 1.02, 0.02)
    x2 = np.arange(0, 1.02, 0.02)
    # g_x = 100 * (8 + 8 * (0.25 - np.cos(-10*np.pi)))
    g_x = 0
    pf = []
    for i in range(len(x1)):
        t1 = x1[i]
        for j in range(len(x2)):
            t2 = x2[j]
            f1 = np.cos(0.5 * np.pi * t1) * np.cos(0.5 * np.pi * t2)
            f2 = np.cos(0.5 * np.pi * t1) * np.sin(0.5 * np.pi * t2)
            f3 = np.sin(0.5 * np.pi * t1)
            s = f1**2 + f2**2 + f3**2
            if np.round(s, 2) == 1:
                pf.append([f1, f2, f3])
    pf = np.array(pf)
    pf = np.round(pf, 4)
    pf = np.unique(pf, axis=0)
    print(pf.shape)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(pf[:,0], pf[:,1], pf[:,2])
    plt.show()
    save_file = 'pf_dtlz3.csv'
    with open(save_file, 'a+') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(pf)

def dtlz4():
    x1 = np.arange(0, 1.02, 0.02)
    x2 = np.arange(0, 1.02, 0.02)
    pf = []
    for i in range(len(x1)):
        t1 = np.power(x1[i], 100)
        for j in range(len(x2)):
            t2 = np.power(x2[j], 100)
            f1 = np.cos(0.5 * np.pi * t1) * np.cos(0.5 * np.pi * t2)
            f2 = np.cos(0.5 * np.pi * t1) * np.sin(0.5 * np.pi * t2)
            f3 = np.sin(0.5 * np.pi * t1)
            s = f1**2 + f2**2 + f3**2
            if np.round(s, 2) == 1:
                pf.append([f1, f2, f3])
    pf = np.array(pf)
    pf = np.round(pf, 4)
    pf = np.unique(pf, axis=0)
    print(pf.shape)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(pf[:,0], pf[:,1], pf[:,2])
    plt.show()
    save_file = 'pf_dtlz4.csv'
    with open(save_file, 'a+') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(pf)

def dtlz5():
    x1 = np.arange(0, 1.02, 0.02)
    x2 = np.arange(0, 1.02, 0.02)
    pf = []
    for i in range(len(x1)):
        t1 = x1[i]
        theta1 = t1
        for j in range(len(x2)):
            t2 = x2[j]
            theta2 = 0.25 * np.pi
            f1 = np.cos(0.5 * np.pi * theta1) * np.cos(0.5 * np.pi * theta2)
            f2 = np.cos(0.5 * np.pi * theta1) * np.sin(0.5 * np.pi * theta2)
            f3 = np.sin(0.5 * np.pi * theta1)
            s = f1**2 + f2**2 + f3**2
            if np.round(s, 2) == 1:
                pf.append([f1, f2, f3])
    pf = np.array(pf)
    pf = np.round(pf, 4)
    pf = np.unique(pf, axis=0)
    print(pf.shape)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(pf[:,0], pf[:,1], pf[:,2])
    plt.show()
    save_file = 'pf_dtlz5.csv'
    with open(save_file, 'a+') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(pf)

def dtlz6():
    x1 = np.arange(0, 1.02, 0.02)
    x2 = np.arange(0, 1.02, 0.02)
    pf = []
    for i in range(len(x1)):
        t1 = x1[i]
        theta1 = t1
        for j in range(len(x2)):
            t2 = x2[j]
            theta2 = 0.25 * np.pi
            f1 = np.cos(0.5 * np.pi * theta1) * np.cos(0.5 * np.pi * theta2)
            f2 = np.cos(0.5 * np.pi * theta1) * np.sin(0.5 * np.pi * theta2)
            f3 = np.sin(0.5 * np.pi * theta1)
            s = f1**2 + f2**2 + f3**2
            if np.round(s, 2) == 1:
                pf.append([f1, f2, f3])
    pf = np.array(pf)
    pf = np.round(pf, 4)
    pf = np.unique(pf, axis=0)
    print(pf.shape)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(pf[:,0], pf[:,1], pf[:,2])
    plt.show()
    save_file = 'pf_dtlz6.csv'
    with open(save_file, 'a+') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(pf)

def dtlz7():
    x1 = np.arange(0, 1.02, 0.02)
    x2 = np.arange(0, 1.02, 0.02)
    pf = []
    for i in range(len(x1)):
        f1 = x1[i]
        for j in range(len(x2)):
            f2 = x2[j]
            f3 = 2 * (3 - ((f1/2)*(1+np.sin(3*np.pi*f1)) + (f2/2)*(1+np.sin(3*np.pi*f2))))
            pf.append([f1, f2, f3])
    pf = np.array(pf)
    pf = np.round(pf, 4)
    pf = np.unique(pf, axis=0)
    ranks = fast_nondominated_sort(pf)
    rpf = []
    for j in range(len(ranks)):
        if ranks[j] == 0:
            rpf.append(pf[j])
    rpf = np.array(rpf)
    print(rpf.shape)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(rpf[:,0], rpf[:,1], rpf[:,2])
    plt.show()
    save_file = 'pf_dtlz7.csv'
    with open(save_file, 'a+') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(pf)


def uf1():
    x1 = np.arange(0, 1.02, 0.02)
    x1 = np.round(x1, 2)
    x = [x1.tolist()]
    for i in range(2, 11):
        x2 = np.sin(6*np.pi*x1+(i*np.pi)/10)
        x2 = np.round(x2, 2)
        x = x + [list(x2)]
    x = np.array(x)
    x = x.T
    n = x.shape[1]
    print(x.shape)
    for ind in x:
        c1 = 0
        c2 = 0
        s1 = 0
        s2 = 0
        for i in range(1, len(ind), 2):
            t1 = ind[i] - np.sin(6 * np.pi * ind[0] + (i + 1) * np.pi / n)
            s1 = s1 + np.power(t1, 2)
            c1 = c1 + 1
        # i表示所有奇数的下标
        for j in range(2, len(ind), 2):
            t2 = ind[j] - np.sin(6 * np.pi * ind[0] + (j + 1) * np.pi / n)
            s2 = s2 + np.power(t2, 2)
            c2 = c2 + 1
        f1 = ind[0] + (2 / c2) * s2
        f2 = 1 - np.sqrt(ind[0]) + (2 / c1) * s1
        s = np.sqrt(f1)+f2
        if np.round(s, 2) == 1.0:
            print((ind*100) % 2)
            flag = (ind*100) % 2
            sf = np.sum(flag)
            sf = np.round(sf, 2)
            print(sf)




# uf1()
# dtlz4()

uf1()