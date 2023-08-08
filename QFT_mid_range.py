# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 10:40:12 2023

@author: LENOVO
"""
import Q_p
from Q_p_plot import Q_p_plot
import numpy as np
import matplotlib.pyplot as plt
import math
import cmath

p = 23
L = 2
Z_L = Q_p.Z_N_group(p, L)
Z_L_list = [a for a in Z_L]

# 2 dimension
alpha = 1

i = complex(0, 1)


# def psi_2(x, t):
#     if x == 0:
#         return 0

#     else:
#         a_1 = p**(-L/2)*Q_p.norm_p(x, p)**(-1)
#         a_2 = cmath.exp(-i*(t*Q_p.norm_p(x, p)**(-alpha)))
#         a_3 = cmath.exp(-i*(t*Q_p.norm_p(x, p)**(-alpha)*p**(alpha)))
#         return a_1*(a_2-a_3)

M = 25


def psi_2(x, t):
    if x == 0:
        return 0
    else:
        a_1 = p**(-L/2)*Q_p.norm_p(x, p)**(-1)
        a_2 = cmath.exp(-i*(t*Q_p.norm_p(x, p)**(-alpha)*p**(alpha)))
        a_3 = 0
        for m in range(M):
            a_3 += (p**(-m)) * cmath.exp(-i *
                                         (t*Q_p.norm_p(x, p)**(-alpha)*p**(-alpha*m)))
        return a_1*((1-p**(-1))*a_3-a_2)


b = Z_L_list[1]
a = Z_L_list[-1]


def psi_mid(x, t):
    return abs(psi_2(x-a, t)+psi_2(x-b, t))**2


# new object plot
plot = Q_p_plot(psi_mid, Z_L)

# Plot
t_min = 0
t_max = 3
delta_t = 0.05
image = plot.plot2var(t_min=t_min, t_max=t_max, delta_t=delta_t,
                      ylabels=False, xlabels=False,
                      with_title=False, title="",
                      with_tree=True,
                      size=(8, 8), cmap=None)

plt.show()
# Plot one point
x_0 = Z_L_list[10]


def psi_mid_point(t):
    return psi_mid(x_0, t)


times = np.arange(t_min, t_max, delta_t)

position = np.array([psi_mid_point(t) for t in times])

image_x_0 = plt.plot(times, position, lw=2)
plt.xlabel("Time")
plt.grid(True)
plt.show()

# Save image
# imag_heat_map = image.get_figure()
# imag_heat_map.savefig('simulation_1_state_r-2.png',
#                       dpi="figure",
#                       bbox_inches="tight")
