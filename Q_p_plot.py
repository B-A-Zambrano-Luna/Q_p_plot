# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 16:01:20 2023

@author: LENOVO
"""

import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch
import pylab
from p_adic_aux_function import \
    vectorize_function,\
    p_adic_Convolution_matrix
import Q_p


class Q_p_plot(object):
    def __init__(self, f, Z_k):
        self.f = f
        self.Z_k = Z_k
        self.conv = []
        self.vect = []
        self.var2 = []

    def convolution(self):
        self.conv = p_adic_Convolution_matrix(self.f, self.Z_k)

    def vector(self):
        self.vect = np.array([vectorize_function(self.f, self.Z_k)])

    def varibles2(self, t_min, t_max, delta_t):
        dt = t_min
        def f_t(x): return self.f(x, dt)
        self.var2 = np.array([vectorize_function(f_t, self.Z_k)])
        dt = dt + delta_t
        while dt <= t_max:
            def f_t(x): return self.f(x, dt)
            f_dt_x = np.array([vectorize_function(f_t, self.Z_k)])
            self.var2 = np.concatenate([self.var2, f_dt_x], axis=0)
            dt = dt + delta_t
        self.var2 = np.rot90(self.var2, k=1, axes=(0, 1))

    def plot_conv(self, ylabels=True, xlabels=True,
                  with_title=True, title="",
                  with_tree=False,
                  size=(8, 8), cmap="magma"):
        if self.conv == []:
            self.convolution()
        p = self.Z_k.get_prime()
        if with_tree:
            """ Generate distance matrix. """
            D = np.zeros([len(self.Z_k), len(self.Z_k)])
            for i in self.Z_k:
                for j in self.Z_k:
                    D[i, j] = Q_p.norm_p(i-j, p)
            condensedD = squareform(D)
            """ Compute and plot first dendrogram. """
            # First Tree
            first_map_position = [0.1, 0.25, 0.2, 0.5]
            fig_tree = pylab.figure(figsize=size)
            ax1 = fig_tree.add_axes(first_map_position)
            Y1 = sch.linkage(condensedD)
            Z1 = sch.dendrogram(Y1, orientation='left')
            ax1.set_xticks([])
            ax1.set_xticklabels([], rotation='vertical')
            ax1.set_yticks([])
            # ax1.set_ylabel("G_"+str(self.Z_k.get_radio()), loc="center")
            ax1.set_ylabel("")
            # Second Tree
            second_map_position = [0.3, 0.05, 0.56, 0.2]
            ax2 = fig_tree.add_axes(second_map_position)
            Z2 = sch.dendrogram(Y1, orientation='bottom')
            ax2.set_xticks([])
            ax2.set_xticklabels([], rotation='vertical')
            ax2.set_yticks([])
            # ax2.set_ylabel("G_"+str(self.Z_k.get_radio()), loc="center")
            fig_tree.add_axes([0.3, 0.25, 0.7, 0.5])
        Z_k_list = list(self.Z_k)
        ax = sns.heatmap(self.conv, cmap=cmap)
        if with_tree:
            ax.set(xlabel="",
                   ylabel="")
        else:
            ax.set(xlabel="p-adic integers",
                   ylabel="p-adic integers")
        if with_title and title == "":
            plt.title("Kernel ", fontsize=13)
        elif with_title and title != "":
            plt.title(title)
        else:
            plt.title("")
        if xlabels and ylabels:
            x_ticks = ax.xaxis.get_ticklocs()
            y_ticks = ax.yaxis.get_ticklocs()
            step = round(len(self.Z_k)/len(x_ticks))
            ax.xaxis.set_ticklabels([Z_k_list[a*step]
                                     for a in range(len(x_ticks))],
                                    rotation=90)
            ax.yaxis.set_ticklabels([Z_k_list[a*step]
                                     for a in range(len(y_ticks))],
                                    rotation=0)
        elif not xlabels and ylabels:
            ax.xaxis.set_ticklabels([])
            y_ticks = ax.yaxis.get_ticklocs()
            step = round(len(self.Z_k)/len(y_ticks))
            ax.yaxis.set_ticklabels([Z_k_list[a*step]
                                     for a in range(len(y_ticks))],
                                    rotation=0)
        elif not ylabels and xlabels:
            x_ticks = ax.xaxis.get_ticklocs()
            step = round(len(self.Z_k)/len(x_ticks))
            ax.xaxis.set_ticklabels([Z_k_list[a*step]
                                     for a in range(len(x_ticks))],
                                    rotation=0)
            ax.yaxis.set_ticklabels([])

        elif not ylabels and not xlabels:
            ax.yaxis.set_ticklabels([])
            ax.xaxis.set_ticklabels([])

        if with_tree:
            return fig_tree
        else:
            return ax

    def plot_vect(self, ylabels=True, xlabels=True,
                  with_title=True, title="",
                  with_tree=False,
                  size=(8, 8), cmap="magma"):
        if self.vect == []:
            self.vector()
        p = self.Z_k.get_prime()
        if with_tree:
            """ Generate distance matrix. """
            D = np.zeros([len(self.Z_k), len(self.Z_k)])
            for i in self.Z_k:
                for j in self.Z_k:
                    D[i, j] = Q_p.norm_p(i-j, p)
            condensedD = squareform(D)
            """ Compute and plot first dendrogram. """
            orientation = "top"
            first_map_position = [0.3, 0.16, 0.477, 0.13]
            add_axis = [0.3, 0.1, 0.6, 0.06]
            fig_tree = pylab.figure(figsize=size)
            ax1 = fig_tree.add_axes(first_map_position)
            Y1 = sch.linkage(condensedD)

            Z1 = sch.dendrogram(Y1, orientation=orientation)
            ax1.set_xticks([])
            ax1.set_xticklabels([], rotation='vertical')
            ax1.set_yticks([])
            # ax1.set_ylabel("G_"+str(self.Z_k.get_radio()), loc="center")
            ax1.set_ylabel("")
            if with_title:
                ax1.set_title("Position", loc="center", fontsize=13)
            else:
                ax1.set_title("", loc="center")
            fig_tree.add_axes(add_axis)
            ax = sns.heatmap(self.vect, cmap=cmap)
        else:
            fig, ax = plt.subplots(figsize=(12, 2))
            ax = sns.heatmap(self.vect, ax=ax, cmap=cmap)

        if with_title and title == "":
            plt.title("Function ")
        elif with_title and title != "":
            plt.title(title)
        else:
            plt.title("")
        if with_tree:
            ax.set(xlabel="",
                   ylabel="")
        else:
            ax.set(xlabel="p-adic integers",
                   ylabel="")
            ax.yaxis.set_ticks([])
        if xlabels:
            Z_k_list = list(self.Z_k)
            ax.xaxis.set_ticklabels(Z_k_list,
                                    rotation=0)
        else:
            ax.xaxis.set_ticklabels([])

        if with_tree:
            return fig_tree
        else:
            return ax

    def plot2var(self, t_min=0, t_max=1, delta_t=0.05,
                 ylabels=True, xlabels=True,
                 with_title=True, title="",
                 with_tree=False,
                 size=(8, 8), cmap="magma"):
        if self.var2 == []:
            self.varibles2(t_min, t_max, delta_t)
        p = self.Z_k.get_prime()

        if with_tree:
            """ Generate distance matrix. """
            D = np.zeros([len(self.Z_k), len(self.Z_k)])
            for i in self.Z_k:
                for j in self.Z_k:
                    D[i, j] = Q_p.norm_p(i-j, p)
            condensedD = squareform(D)
            """ Compute and plot first dendrogram. """
            first_map_position = [0.1, 0.1, 0.2, 0.5]
            fig_tree = pylab.figure(figsize=size)
            ax1 = fig_tree.add_axes(first_map_position)
            Y1 = sch.linkage(condensedD)
            Z1 = sch.dendrogram(Y1, orientation='left')
            ax1.set_xticks([])
            ax1.set_xticklabels([], rotation='vertical')
            ax1.set_yticks([])
            ax1.set_ylabel("Position", fontsize=18)
            fig_tree.add_axes([0.3, 0.1, 1.0, 0.5])
        ax = sns.heatmap(self.var2, cmap=cmap)
        if with_tree:
            ax.set_xlabel("Time", fontsize=18)
            ax.set_ylabel("")
        else:
            ax.set(xlabel="Time",
                   ylabel="p-adic integers")
        if with_title and title == "":
            plt.title("Function")
        elif with_title and title != "":
            plt.title(title)
        else:
            plt.title("")
        x_ticks = ax.xaxis.get_ticklocs()
        if xlabels:
            t_step = int((t_max-t_min)/delta_t)
            step = t_step/(len(x_ticks)-1)
            ax.xaxis.set_ticklabels([round(t*step*delta_t+t_min, 2)
                                     for t in range(len(x_ticks))])
        else:
            ax.xaxis.set_ticklabels([])

        if ylabels:
            y_ticks = ax.yaxis.get_ticklocs()
            Z_k_list = list(self.Z_k)
            step = round(len(self.Z_k)/len(y_ticks))
            y_labels = [Z_k_list[a*step] for a in range(len(y_ticks))]
            y_labels.reverse()
            ax.yaxis.set_ticklabels(y_labels,
                                    rotation=0)
        else:
            ax.yaxis.set_ticklabels([], rotation=0)

        if with_tree:
            return fig_tree
        else:
            return ax
