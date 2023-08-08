# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 16:12:35 2023

@author: LENOVO
"""

import time
import Q_p
import test_functions as test
import math
import numpy as np
from Q_p_plot import Q_p_plot

p = 3
K = 8  # max 8 for 1 min memory
Z_k = Q_p.Z_N_group(p, K)

# Parameters
b_ee = 1.5
b_ei = 1.35
b_ie = 1.35
b_ii = 1.8
sigma_ee = 40/10
sigma_ei = 60/10
sigma_ie = 60/10
sigma_ii = 30/10
# Kernels


def W_ee(x):
    return b_ee*math.exp(sigma_ee)-b_ee*math.exp(Q_p.norm_p(x, p)*sigma_ee)


start = time.time()

kernel = Q_p_plot(W_ee, Z_k)

kernel.plot_conv()

# kernel.plot_vect()

end = time.time()
print("Time ejecution solver", end-start)
