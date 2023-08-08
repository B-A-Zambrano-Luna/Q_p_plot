# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 10:40:12 2023

@author: LENOVO
"""
import Q_p
from Q_p_plot import Q_p_plot
import numpy as np


p = 3
L = 3
Z_L = Q_p.Z_N_group(p, L)

# 2 dimension
alpha = 1


def f(z, t):
    if z != 0:
        a_0 = ((1-p**(alpha))**2)*(t**2)*(p**(-L))
        a_1 = (p**(4*L*(1+alpha)))\
            * (Q_p.norm_p(z, p)**(2*(1+alpha)))
        a_3 = np.sinc(t*p**(2*L*alpha)*((1-p**alpha)/2)
                      * Q_p.norm_p(z, p)**(-alpha))
    elif z == 0:
        return 0

    return (a_0/a_1)*(a_3**2)


# new object plot
plot = Q_p_plot(f, Z_L)

# Plot
image = plot.plot2var(t_min=0, t_max=1, delta_t=0.05,
                      ylabels=False, xlabels=True,
                      with_title=False, title="",
                      with_tree=True,
                      size=(8, 8), cmap=None)

# Save image
# imag_heat_map = image.get_figure()
# imag_heat_map.savefig('simulation_1_state_r-2.png',
#                       dpi="figure",
#                       bbox_inches="tight")
