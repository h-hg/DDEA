import numpy as np
from pymoo.problems.many.wfg import WFG1, WFG2, WFG5

## benchmark problems, inclusing ZDT, DTLZ, UF, WFG
## We can define two different types of noises, multiplicative and additive noises.

def ZDT1(ind, t):
    # control the contamination level
    l = 0.2
    # 真实的函数值
    f1 = ind[0]
    s = 0
    for i in range(1, len(ind)):
        s = s + ind[i]
    g_x = 1 + 9 * s / (len(ind) - 1)
    h_x = 1 - np.sqrt(f1 / g_x)
    f2 = h_x * g_x
    # additive noises
    # e1 = np.random.normal(0, np.abs(1.0 * l), size=t)
    # e2 = np.random.normal(0, np.abs(10.0 * l), size=t)
    # multiplicative noises
    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    # print(e2)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)
    return r1, r2

def ZDT2(ind, t):
    # control the contamination level
    l = 0.2
    # 真实的函数值
    f1 = ind[0]
    s = 0
    for i in range(1, len(ind)):
        s = s + ind[i]
    g_x = 1 + 9 * s / (len(ind) - 1)
    h_x = 1 - np.square((f1 / g_x))
    # h_x = 1 - np.sqrt(f1 / g_x)
    f2 = h_x * g_x
    # additive noises
    # e1 = np.random.normal(0, 1.0 * l, size=t)
    # e2 = np.random.normal(0, 10.0 * l, size=t)
    # multiplicative noises
    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)
    return r1, r2

def ZDT3(ind, t):
    # 噪音的扰乱因
    l = 0.2
    # 真实的函数值
    f1 = ind[0]
    s = 0
    for i in range(1, len(ind)):
        s = s + ind[i]
    g_x = 1 + 9 * s / (len(ind) - 1)
    h_x = 1 - np.sqrt(f1 / g_x) - (f1 / g_x) * np.sin(10*np.pi*f1)
    f2 = h_x * g_x
    # 添加噪音
    # e1 = np.random.normal(0, 1.0 * l, size=t)
    # e2 = np.random.normal(0, 10 * l, size=t)
    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    # print(e2)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)
    return r1, r2

def ZDT4(ind, t):
    # 噪音的扰乱因子
    l = 0.2
    f1 = ind[0]
    s = 0
    for i in range(1, len(ind)):
        s = s + (np.square(ind[i]) - 10 * np.cos(4*np.pi*ind[i]))
    g_x = 1 + 10 * (len(ind) -1) + s
    h_x = 1 - np.sqrt(f1 / g_x)
    f2 = h_x * g_x
    # thita1 = 1.0 * l
    # e1 = np.random.normal(0, thita1, size=t)
    # if len(ind) == 10:
    #     e2 = np.random.normal(0, 316 * l, size=t)
    # elif len(ind) == 30:
    #     e2 = np.random.normal(0, 1000 * l, size=t)
    # elif len(ind) == 50:
    #     e2 = np.random.normal(0, 1700 * l, size=t)
    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)
    return r1, r2


def ZDT6(ind, t):
    # 噪音的扰乱因子
    l = 0.2
    f1 = 1 - np.exp((-4*ind[0]))*np.power(np.sin(6*np.pi*ind[0]), 6)
    s = 0
    for i in range(1, len(ind)):
        s = s + ind[i]
    g_x = 1 + 9 * np.power((s / (len(ind) - 1)), 0.25)
    h_x = 1 - np.square((f1 / g_x))
    f2 = h_x * g_x
    # 添加噪音
    # e1 = np.random.normal(0, 0.8 * l, size=t)
    # e2 = np.random.normal(0, 10 * l, size=t)
    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    # print(e2)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)
    return r1, r2


def DTLZ1(ind, t):
    s = 0
    l = 0.2
    for i in range(2, len(ind)):
        s = s + np.square(ind[i] - 0.5) - np.cos(20*np.pi*(ind[i] - 0.5))
    g_x = 100 * (len(ind) - 2 + s)
    f1 = 0.5 * ind[0] * ind[1] * (1 + g_x)
    f2 = 0.5 * ind[0] *(1-ind[1]) * (1 + g_x)
    f3 = 0.5 * (1-ind[0]) * (1 + g_x)
    # thita1 = 500 * l
    # thita2 = 1000 * l
    # thita3 = 2500 * l
    #
    # if len(ind) == 10:
    #     e1 = np.random.normal(0, thita1, size=t)
    #     e2 = np.random.normal(0, thita1, size=t)
    #     e3 = np.random.normal(0, thita1, size=t)
    # elif len(ind) == 30:
    #     e1 = np.random.normal(0, thita2, size=t)
    #     e2 = np.random.normal(0, thita2, size=t)
    #     e3 = np.random.normal(0, thita2, size=t)
    # elif len(ind) == 50:
    #     e1 = np.random.normal(0, thita3, size=t)
    #     e2 = np.random.normal(0, thita3, size=t)
    #     e3 = np.random.normal(0, thita3, size=t)
    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    e3 = np.random.normal(0, np.abs(f3 * l), size=t)

    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r3 = f3 + np.sum(e3) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)
    r3 = np.round(r3, 2)
    return r1, r2, r3


def DTLZ2(ind, t):
    l = 0.2
    s = 0
    for i in range(2, len(ind)):
        s = s + np.square(ind[i] - 0.5)
    f1 = np.cos(0.5*np.pi*ind[0])*np.cos(0.5*np.pi*ind[1])*(1+s)
    f2 = np.cos(0.5*np.pi*ind[0])*np.sin(0.5*np.pi*ind[1])*(1+s)
    f3 = np.sin(0.5*np.pi*ind[0])*(1+s)

    # if len(ind) == 10:
    #     e1 = np.random.normal(0, 3 * l, size=t)
    #     e2 = np.random.normal(0, 3 * l, size=t)
    #     e3 = np.random.normal(0, 3 * l, size=t)
    # elif len(ind) == 30:
    #     e1 = np.random.normal(0, 8 * l, size=t)
    #     e2 = np.random.normal(0, 8 * l, size=t)
    #     e3 = np.random.normal(0, 8 * l, size=t)
    # elif len(ind) == 50:
    #     e1 = np.random.normal(0, 13 * l, size=t)
    #     e2 = np.random.normal(0, 13 * l, size=t)
    #     e3 = np.random.normal(0, 13 * l, size=t)
    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    e3 = np.random.normal(0, np.abs(f3 * l), size=t)

    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r3 = f3 + np.sum(e3) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)
    r3 = np.round(r3, 2)
    return r1, r2, r3


def DTLZ3(ind, t):
    l = 0.2
    s = 0
    for i in range(2, len(ind)):
        s = s + np.square(ind[i] - 0.5) - np.cos(20 * np.pi * (ind[i] - 0.5))
    g_x = 100 * (len(ind) - 2 + s)
    f1 = (1+g_x)*np.cos(0.5*np.pi*ind[0])*np.cos(0.5*np.pi*ind[1])
    f2 = (1+g_x)*np.cos(0.5*np.pi*ind[0])*np.sin(0.5*np.pi*ind[1])
    f3 = (1+g_x)*np.sin(0.5*np.pi*ind[0])

    # if len(ind) == 10:
    #     e1 = np.random.normal(0, 1200 * l, size=t)
    #     e2 = np.random.normal(0, 1200 * l, size=t)
    #     e3 = np.random.normal(0, 1200 * l, size=t)
    # elif len(ind) == 30:
    #     e1 = np.random.normal(0, 3500 * l, size=t)
    #     e2 = np.random.normal(0, 3500 * l, size=t)
    #     e3 = np.random.normal(0, 3500 * l, size=t)
    # elif len(ind) == 50:
    #     e1 = np.random.normal(0, 6000 * l, size=t)
    #     e2 = np.random.normal(0, 6000 * l, size=t)
    #     e3 = np.random.normal(0, 6000 * l, size=t)
    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    e3 = np.random.normal(0, np.abs(f3 * l), size=t)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r3 = f3 + np.sum(e3) / t
    r1 = np.round(r1, 3)
    r2 = np.round(r2, 3)
    r3 = np.round(r3, 3)
    return r1, r2, r3


def DTLZ4(ind, t):
    l = 0.2
    s = 0
    for i in range(2, len(ind)):
        s = s + np.square((ind[i] - 0.5))
    f1 = np.cos(0.5*np.pi*np.power(ind[0], 100))*np.cos(0.5*np.pi*np.power(ind[1], 100))*(1+s)
    f2 = np.cos(0.5*np.pi*np.power(ind[0], 100))*np.sin(0.5*np.pi*np.power(ind[1], 100))*(1+s)
    f3 = np.sin(0.5*np.pi*np.power(ind[0], 100))*(1+s)

    # if len(ind) == 10:
    #     e1 = np.random.normal(0, 2.5 * l, size=t)
    #     e2 = np.random.normal(0, 2.5 * l, size=t)
    #     e3 = np.random.normal(0, 0.001 * l, size=t)
    # elif len(ind) == 30:
    #     e1 = np.random.normal(0, 6 * l, size=t)
    #     e2 = np.random.normal(0, 6 * l, size=t)
    #     e3 = np.random.normal(0, 6 * l, size=t)
    # elif len(ind) == 50:
    #     e1 = np.random.normal(0, 8 * l, size=t)
    #     e2 = np.random.normal(0, 8 * l, size=t)
    #     e3 = np.random.normal(0, 8 * l, size=t)

    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    e3 = np.random.normal(0, np.abs(f3 * l), size=t)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r3 = f3 + np.sum(e3) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)
    r3 = np.round(r3, 2)
    return r1, r2, r3

def DTLZ5(ind, t):
    l = 0.2
    s = 0

    for i in range(2, len(ind)):
        s = s + np.square((ind[i] - 0.5))
    theta = np.pi / (4 * (1 + s)) * (1 + 2 * s * ind[1])
    f1 = np.cos(0.5 * np.pi * ind[0]) * np.cos(0.5 * np.pi * theta) * (1 + s)
    f2 = np.cos(0.5 * np.pi * ind[0]) * np.sin(0.5 * np.pi * theta) * (1 + s)
    f3 = np.sin(0.5 * np.pi * ind[0]) * (1 + s)
    #
    # if len(ind) == 10:
    #     e1 = np.random.normal(0, 3.0 * l, size=t)
    #     e2 = np.random.normal(0, 3.0 * l, size=t)
    #     e3 = np.random.normal(0, 3.5 * l, size=t)
    # elif len(ind) == 30:
    #     e1 = np.random.normal(0, 5.0 * l, size=t)
    #     e2 = np.random.normal(0, 5.0 * l, size=t)
    #     e3 = np.random.normal(0, 5.5 * l, size=t)
    # elif len(ind) == 50:
    #     e1 = np.random.normal(0, 7.0 * l, size=t)
    #     e2 = np.random.normal(0, 7.0 * l, size=t)
    #     e3 = np.random.normal(0, 7.5 * l, size=t)
    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    e3 = np.random.normal(0, np.abs(f3 * l), size=t)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r3 = f3 + np.sum(e3) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)
    r3 = np.round(r3, 2)


    return r1, r2, r3

def DTLZ6(ind, t):
    l = 0.2
    s = 0
    #
    for i in range(2, len(ind)):
        s = s + np.power(ind[i], 0.1)
    theta = np.pi / (4 * (1 + s)) * (1 + 2 * s * ind[1])
    f1 = np.cos(0.5 * np.pi * ind[0]) * np.cos(0.5 * np.pi * theta) * (1 + s)
    f2 = np.cos(0.5 * np.pi * ind[0]) * np.sin(0.5 * np.pi * theta) * (1 + s)
    f3 = np.sin(0.5 * np.pi * ind[0]) * (1 + s)
    #
    # if len(ind) == 10:
    #     e1 = np.random.normal(0, 10 * l, size=t)
    #     e2 = np.random.normal(0, 10 * l, size=t)
    #     e3 = np.random.normal(0, 10 * l, size=t)
    # elif len(ind) == 30:
    #     e1 = np.random.normal(0, 30.0 * l, size=t)
    #     e2 = np.random.normal(0, 30.0 * l, size=t)
    #     e3 = np.random.normal(0, 30.5 * l, size=t)
    # elif len(ind) == 50:
    #     e1 = np.random.normal(0, 50.0 * l, size=t)
    #     e2 = np.random.normal(0, 50.0 * l, size=t)
    #     e3 = np.random.normal(0, 50.5 * l, size=t)
    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    e3 = np.random.normal(0, np.abs(f3 * l), size=t)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r3 = f3 + np.sum(e3) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)
    r3 = np.round(r3, 2)


    return r1, r2, r3


def DTLZ7(ind, t):
    l = 0.2
    s = 0
    #
    for i in range(2, len(ind)):
        s = s + ind[i]
    g = 1 + 9 / (len(ind) - 2) * s
    f1 = ind[0]
    f2 = ind[1]
    h = 3 - (f1 / (1 + g)) * (1 + np.sin(3 * np.pi * f1)) - (f2 / (1 + g)) * (1 + np.sin(3 * np.pi * f2))
    f3 = (1 + g) * h
    #
    # if len(ind) == 10:
    #     e1 = np.random.normal(0, 30.0 * l, size=t)
    #     e2 = np.random.normal(0, 30.0 * l, size=t)
    #     e3 = np.random.normal(0, 30.0 * l, size=t)
    # elif len(ind) == 30:
    #     e1 = np.random.normal(0, 30.0 * l, size=t)
    #     e2 = np.random.normal(0, 30.0 * l, size=t)
    #     e3 = np.random.normal(0, 30.0 * l, size=t)
    # elif len(ind) == 50:
    #     e1 = np.random.normal(0, 30.0 * l, size=t)
    #     e2 = np.random.normal(0, 30.0 * l, size=t)
    #     e3 = np.random.normal(0, 30.0 * l, size=t)
    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    e3 = np.random.normal(0, np.abs(f3 * l), size=t)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r3 = f3 + np.sum(e3) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)
    r3 = np.round(r3, 2)


    return r1, r2, r3



def uf1(ind, t=1):
    l = 0.2
    n = len(ind)
    c1 = 0
    c2 = 0
    s1 = 0
    s2 = 0
    # i表示所有偶数的下标
    for i in range(1, len(ind), 2):
        t1 = ind[i] - np.sin(6*np.pi*ind[0]+(i+1)*np.pi/n)
        s1 = s1 + np.power(t1, 2)
        c1 = c1 + 1
    # i表示所有奇数的下标
    for j in range(2, len(ind), 2):
        t2 = ind[j] - np.sin(6 * np.pi * ind[0] + (j + 1) * np.pi / n)
        s2 = s2 + np.power(t2, 2)
        c2 = c2 + 1
    f1 = ind[0] + (2/c2) * s2
    f2 = 1 - np.sqrt(ind[0]) + (2/c1) * s1

    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)

    return r1, r2


def uf2(ind, t=1):
    l = 0.2
    n = len(ind)
    c1 = 0
    c2 = 0
    s1 = 0
    s2 = 0
    # i表示所有偶数的下标
    for i in range(1, len(ind), 2):
        t1 = ind[i] - (0.3*np.power(ind[0], 2)*np.cos(24*np.pi*ind[0]+4*(i+1)*np.pi/n)+0.6*ind[0])*np.sin(6*np.pi*ind[0]+(i+1)*np.pi/n)
        s1 = s1 + np.power(t1, 2)
        c1 = c1 + 1

    # i表示所有奇数的下标
    for j in range(2, len(ind), 2):
        t2 = ind[j] - (0.3*np.power(ind[0], 2)*np.cos(24*np.pi*ind[0]+4*(j+1)*np.pi/n)+0.6*ind[0])*np.cos(6*np.pi*ind[0]+(j+1)*np.pi/n)
        s2 = s2 + np.power(t2, 2)
        c2 = c2 + 1

    f1 = ind[0] + (2 / c2) * s2
    f2 = 1 - np.sqrt(ind[0]) + (2 / c1) * s1

    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)

    return r1, r2


def uf4(ind, t=1):
    l = 0.2
    n = len(ind)
    c1 = 0
    c2 = 0
    s1 = 0
    s2 = 0
    # i表示所有偶数的下标
    for i in range(1, len(ind), 2):
        y = ind[i] - np.sin(6*np.pi*ind[0]+(i+1)*np.pi/n)
        h = np.abs(y) / (1+np.exp(2*np.abs(y)))
        s1 = s1 + h
        c1 = c1 + 1

    # i表示所有奇数的下标
    for j in range(2, len(ind), 2):
        y = ind[j] - np.sin(6 * np.pi * ind[0] + (j + 1) * np.pi / n)
        h = np.abs(y) / (1 + np.exp(2 * np.abs(y)))
        s2 = s2 + h
        c2 = c2 + 1

    f1 = ind[0] + (2 / c2) * s2
    f2 = 1 - np.power(ind[0], 2) + (2 / c1) * s1

    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)

    return r1, r2



def uf6(ind, t=1):
    l = 0.2
    n = len(ind)
    c1 = 0
    c2 = 0
    s1 = 0
    s2 = 0
    p1 = 1
    p2 = 1
    # i表示所有偶数的下标
    for i in range(1, len(ind), 2):
        y = ind[i] - np.sin(6*np.pi*ind[0]+(i+1)*np.pi/n)
        s1 = s1 + np.power(y, 2)
        p1 = p1 * np.cos(20*y*np.pi/np.sqrt(i+1))
        c1 = c1 + 1

    # i表示所有奇数的下标
    for j in range(2, len(ind), 2):
        y = ind[j] - np.sin(6 * np.pi * ind[0] + (j + 1) * np.pi / n)
        s2 = s2 + np.power(y, 2)
        p2 = p2 * np.cos(20 * y * np.pi / np.sqrt(j + 1))
        c2 = c2 + 1

    # odd
    f1 = ind[0] + np.max([0, 0.7*np.sin(4*np.pi*ind[0])]) + (2 / c2) * (4*s2 - 2*p2 + 2)
    # even
    f2 = 1 - ind[0] + np.max([0, 0.7*np.sin(4*np.pi*ind[0])]) + (2 / c1) * (4*s1 - 2*p1 + 2)

    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)

    return r1, r2

def uf8(ind, t=1):
    l = 0.2
    s1 = 0
    s2 = 0
    s3 = 0
    c1 = 0
    c2 = 0
    c3 = 0
    n = len(ind)

    # J3
    for i in range(2, len(ind), 3):
        tmp = ind[i] - 2 * ind[1] * np.sin(2*np.pi*ind[0] + (i+1)*np.pi/n )
        s3 = s3 + np.power(tmp, 2)
        c3 = c3 + 1
    # J2
    for i in range(4, len(ind), 3):
        tmp = ind[i] - 2 * ind[1] * np.sin(2*np.pi*ind[0] + (i+1)*np.pi/n )
        s3 = s3 + np.power(tmp, 2)
        c2 = c2 + 1
    # J1
    for i in range(3, len(ind), 3):
        tmp = ind[i] - 2 * ind[1] * np.sin(2*np.pi*ind[0] + (i+1)*np.pi/n )
        s1 = s1 + np.power(tmp, 2)
        c1 = c1 + 1

    f1 = np.cos(0.5 * np.pi * ind[0]) * np.cos(0.5 * np.pi * ind[1]) + (2 / c1) * s1
    f2 = np.cos(0.5 * np.pi * ind[0]) * np.sin(0.5 * np.pi * ind[1]) + (2 / c2) * s2
    f3 = np.sin(0.5 * np.pi * ind[0]) + (2 / c3) * s3

    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    e3 = np.random.normal(0, np.abs(f3 * l), size=t)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r3 = f3 + np.sum(e3) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)
    r3 = np.round(r3, 2)

    return r1, r2, r3



def uf9(ind, t=1):
    l = 0.2
    s1 = 0
    s2 = 0
    s3 = 0
    c1 = 0
    c2 = 0
    c3 = 0
    n = len(ind)
    # J3
    for i in range(2, len(ind), 3):
        tmp = ind[i] - 2 * ind[1] * np.sin(2*np.pi*ind[0] + (i+1)*np.pi/n )
        s3 = s3 + np.power(tmp, 2)
        c3 = c3 + 1
    # J2
    for i in range(4, len(ind), 3):
        tmp = ind[i] - 2 * ind[1] * np.sin(2*np.pi*ind[0] + (i+1)*np.pi/n )
        s3 = s3 + np.power(tmp, 2)
        c2 = c2 + 1
    # J1
    for i in range(3, len(ind), 3):
        tmp = ind[i] - 2 * ind[1] * np.sin(2*np.pi*ind[0] + (i+1)*np.pi/n )
        s1 = s1 + np.power(tmp, 2)
        c1 = c1 + 1

    f1 = 0.5 * (np.max([0, 1.1*(1-4*np.power(2*ind[0]-1 ,2))]) + 2*ind[0]) * ind[1] + (2 / c1) * s1
    f2 = 0.5 * (np.max([0, 1.1*(1-4*np.power(2*ind[0]-1 ,2))]) - 2*ind[0] + 2) * ind[1] + (2 / c2) * s2
    f3 = 1 - ind[1] + (2 / c3) * s3

    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    e3 = np.random.normal(0, np.abs(f3 * l), size=t)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r3 = f3 + np.sum(e3) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)
    r3 = np.round(r3, 2)

    return r1, r2, r3



def wfg1(ind, t=1):
    l = 0.2
    wfg = WFG1(n_var=len(ind), n_obj=2)
    y = wfg.evaluate(np.array([ind]))[0]
    f1 = y[0]
    f2 = y[1]
    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)
    return r1, r2

def wfg2(ind, t=1):
    l = 0.2
    wfg = WFG2(n_var=len(ind), n_obj=2)
    y = wfg.evaluate(np.array([ind]))[0]
    f1 = y[0]
    f2 = y[1]
    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)
    return r1, r2


def wfg5(ind, t=1):
    l = 0.2
    wfg = WFG5(n_var=len(ind), n_obj=2)
    y = wfg.evaluate(np.array([ind]))[0]
    f1 = y[0]
    f2 = y[1]
    e1 = np.random.normal(0, np.abs(f1 * l), size=t)
    e2 = np.random.normal(0, np.abs(f2 * l), size=t)
    r1 = f1 + np.sum(e1) / t
    r2 = f2 + np.sum(e2) / t
    r1 = np.round(r1, 2)
    r2 = np.round(r2, 2)
    return r1, r2



def init_pop(l=-5, u=5.01, interval=0.02, size=(300, 10)):
    x2 = np.arange(l, u, interval)
    x1 = np.arange(0, 1, interval)
    pop = np.zeros(shape=size)
    for i in range(pop.shape[0]):
        r1 = np.random.randint(0, x1.shape[0])
        pop[i, 0] = np.round(x1[r1, ], 2)
        for j in range(1, pop.shape[1]):
            r2 = np.random.randint(0, x2.shape[0])
            pop[i, j] = np.round(x2[r2, ], 2)
    return pop


def evaluate_true(pop, p_namme, t):
    for i in range(len(pop)):
        if p_namme == 'zdt1':
            pop[i].fitness.values = ZDT1(pop[i], t)
        elif p_namme == 'zdt2':
            pop[i].fitness.values = ZDT2(pop[i], t)
        elif p_namme == 'zdt3':
            pop[i].fitness.values = ZDT3(pop[i], t)
        elif p_namme == 'zdt4':
            pop[i].fitness.values = ZDT4(pop[i], t)
        elif p_namme == 'zdt6':
            pop[i].fitness.values = ZDT6(pop[i], t)
        elif p_namme == 'dtlz1':
            pop[i].fitness.values = DTLZ1(pop[i], t)
        elif p_namme == 'dtlz2':
            pop[i].fitness.values = DTLZ2(pop[i], t)
        elif p_namme == 'dtlz3':
            pop[i].fitness.values = DTLZ3(pop[i], t)
        elif p_namme == 'dtlz4':
            pop[i].fitness.values = DTLZ4(pop[i], t)
        elif p_namme == 'dtlz5':
            pop[i].fitness.values = DTLZ5(pop[i], t)
        elif p_namme == 'dtlz6':
            pop[i].fitness.values = DTLZ6(pop[i], t)
        elif p_namme == 'dtlz7':
            pop[i].fitness.values = DTLZ7(pop[i], t)
        elif p_namme == 'uf1':
            pop[i].fitness.values = uf1(pop[i], t)
        elif p_namme == 'uf2':
            pop[i].fitness.values = uf2(pop[i], t)
        elif p_namme == 'uf4':
            pop[i].fitness.values = uf4(pop[i], t)
        elif p_namme == 'uf6':
            pop[i].fitness.values = uf6(pop[i], t)
        elif p_namme == 'uf8':
            pop[i].fitness.values = uf8(pop[i], t)
        elif p_namme == 'uf9':
            pop[i].fitness.values = uf9(pop[i], t)
        elif p_namme == 'wfg1':
            pop[i].fitness.values = wfg1(pop[i], t)
        elif p_namme == 'wfg2':
            pop[i].fitness.values = wfg2(pop[i], t)
        elif p_namme == 'wfg5':
            pop[i].fitness.values = wfg5(pop[i], t)




def evaluate_ntspea(pop, p_namme, t):
    fes = 0
    for i in range(len(pop)):
        # if pop[i].fe < 10:
        if p_namme == 'zdt1':
            r1, r2 = ZDT1(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'zdt2':
            r1, r2 = ZDT2(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'zdt3':
            r1, r2 = ZDT3(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'zdt4':
            r1, r2 = ZDT4(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'zdt6':
            r1, r2 = ZDT6(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'dtlz1':
            r1, r2, r3 = DTLZ1(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'dtlz2':
            r1, r2, r3 = DTLZ2(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'dtlz3':
            r1, r2, r3 = DTLZ3(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'dtlz4':
            r1, r2, r3 = DTLZ4(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'dtlz5':
            r1, r2, r3 = DTLZ5(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'dtlz6':
            r1, r2, r3 = DTLZ6(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'dtlz7':
            r1, r2, r3 = DTLZ7(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'uf1':
            r1, r2 = uf1(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'uf2':
            r1, r2 = uf2(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'uf4':
            r1, r2 = uf4(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'uf6':
            r1, r2 = uf6(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'uf8':
            r1, r2, r3 = uf8(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'uf9':
            r1, r2, r3 = uf9(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'wfg1':
            r1, r2 = wfg1(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'wfg2':
            r1, r2 = wfg2(pop[i], t=t)
            fes = fes + t
        elif p_namme == 'wfg5':
            r1, r2 = wfg5(pop[i], t=t)
            fes = fes + t

        pop[i].fe = pop[i].fe + 1
        if len(pop[i].fitness.values) == 0:
            if p_namme[0:3] == 'zdt' or p_namme[0:3] == 'wfg' or (p_namme[0:3]  in ['uf1', 'uf2', 'uf4', 'uf6']):
                pop[i].fitness.values = [r1, r2]
            else:
                pop[i].fitness.values = [r1, r2, r3]
        else:
            temp = pop[i].fitness.values
            t1 = np.round((temp[0] + r1) / 2, 3)
            t2 = np.round((temp[1] + r2) / 2, 3)
            if p_namme[0:3] == 'dtl' or p_namme[0:3]  == 'uf8' or p_namme[0:3]  == 'uf9':
                t3 = np.round((temp[2] + r3) / 2, 3)
                pop[i].fitness.values = [t1, t2, t3]
            else:
                pop[i].fitness.values = [t1, t2]
    return fes



def get_fitness(ind, p_namme, t):
    if p_namme == 'zdt1':
        return ZDT1(ind, t=t)
    elif p_namme == 'zdt2':
        return ZDT2(ind, t=t)
    elif p_namme == 'zdt3':
        return ZDT3(ind, t=t)
    elif p_namme == 'zdt4':
        return ZDT4(ind, t=t)
    elif p_namme == 'zdt6':
        return ZDT6(ind, t=t)
    elif p_namme == 'dtlz1':
        return DTLZ1(ind, t=t)
    elif p_namme == 'dtlz2':
        return DTLZ2(ind, t=t)
    elif p_namme == 'dtlz3':
        return DTLZ3(ind, t=t)
    elif p_namme == 'dtlz4':
        return DTLZ4(ind, t=t)
    elif p_namme == 'dtlz5':
        return DTLZ5(ind, t=t)
    elif p_namme == 'dtlz6':
        return DTLZ6(ind, t=t)
    elif p_namme == 'dtlz7':
        return DTLZ7(ind, t=t)
    elif p_namme == 'uf1':
        return uf1(ind, t=t)
    elif p_namme == 'uf2':
        return uf2(ind, t=t)
    elif p_namme == 'uf4':
        return uf4(ind, t=t)
    elif p_namme == 'uf6':
        return uf6(ind, t=t)
    elif p_namme == 'uf8':
        return uf8(ind, t=t)
    elif p_namme == 'uf9':
        return uf9(ind, t=t)
    elif p_namme == 'wfg1':
        return wfg1(ind, t=t)
    elif p_namme == 'wfg2':
        return wfg2(ind, t=t)
    elif p_namme == 'wfg5':
        return wfg5(ind, t=t)




def dominance(p, q):
    flag = False
    g = p.fitness.values
    f = q.fitness.values
    obj_count = len(p.fitness.values)
    if obj_count == 2:
        if g[0] <= f[0] and g[1] <= f[1]:
             if g[0] < f[0] or g[1] < f[1]:
                    flag = True
    elif obj_count == 3:
            if g[0] <= f[0] and g[1] <= f[1] and g[2] <= f[2]:
                if g[0] < f[0] or g[1] < f[1] or g[2] < f[2]:
                    flag = True
    return flag