# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 10:40:12 2023

@author: LENOVO
"""
import Q_p
from Q_p_plot import Q_p_plot


p = 5
k = 2
Z_k = Q_p.Z_N_group(p, k)

# 1 dimension


def f(x):
    return Q_p.norm_p(x, p)


# object plot
plot = Q_p_plot(f, Z_k)

# As kernel
image = plot.plot_conv(ylabels=False, xlabels=False,
                       with_title=True, title="",
                       with_tree=True,
                       size=(8, 8), cmap=None)

# As function

image = plot.plot_vect(ylabels=False, xlabels=False,
                       with_title=False, title="",
                       with_tree=True,
                       size=(8, 8), cmap=None)


# 2 dimension
def f(x, t):
    return Q_p.norm_p(x, p)*t


# new object plot
plot = Q_p_plot(f, Z_k)

# Plot
image = plot.plot2var(t_min=-1, t_max=1, delta_t=0.5,
                      ylabels=False, xlabels=True,
                      with_title=False, title="",
                      with_tree=True,
                      size=(8, 8), cmap=None)
