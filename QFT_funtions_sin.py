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

p = 3
L = 4
Z_L = Q_p.Z_N_group(p, L)

# 2 dimension
alpha = 1


def psi_2(x, t):
    if x == 0:
        return 0

    else:
        a_1 = (8/2**0.5)*p**(-L)*Q_p.norm_p(x, p)**(-2)
        a_2 = t*Q_p.norm_p(x, p)**(-alpha)*((1-p**(alpha))/2)
        return a_1*math.sin(a_2)**2


b = 1


def psi_mid(x, t):
    return psi_2(x, t)+psi_2(x-b, t)


# new object plot
plot = Q_p_plot(psi_mid, Z_L)

# Plot
t_min = 0
t_max = 1
delta_t = 0.005
image = plot.plot2var(t_min=t_min, t_max=t_max, delta_t=delta_t,
                      ylabels=False, xlabels=False,
                      with_title=False, title="",
                      with_tree=True,
                      size=(8, 8), cmap=None)

plt.show()
# Plot one point
x_0 = p**(L-1)


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
