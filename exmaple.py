# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 13:18:32 2022

@author: LENOVO
"""
import Q_p
import math
from Q_p_plot import plot

m = 0
s = 0.4
p = 3
K = 7
Z_k = Q_p.Z_N_group(p, K)


def A(x):
    return math.exp(-Q_p.norm_p(x, p))


plot(A, Z_k, As="kernel")
plot(A, Z_k, As="convolution")
plot(A, Z_k, As="fractal", s=0.5, m=0)
plot(A, Z_k, As="function")
