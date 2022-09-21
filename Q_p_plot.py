# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 12:17:25 2022

@author: LENOVO
"""
from Q_p import norm_p
from Q_p_as_fractal import Christiakov_emmending
import numpy as np
from scipy.spatial.distance import squareform
from matplotlib import pyplot as plt
import pylab
import scipy.cluster.hierarchy as sch


def p_adic_matrix(g, Z_k):
    p = Z_k.get_prime()
    k = Z_k.get_radio()
    vect_g = []
    index = np.array([i for i in Z_k], dtype=int)
    Dic_index = dict()
    for a in index:
        Dic_index[a] = g(a)
    vect_g = np.array(vect_g)
    matrix_g = np.zeros((p**k, p**k))
    num_row = 0
    for j in index:
        New_Values = (index - j) % p**k
        matrix_g[num_row] = np.array([Dic_index[a] for a in New_Values])
        num_row += 1
    return matrix_g


def p_adic_Convolution_matrix(g, Z_k):
    p = Z_k.get_prime()
    k = Z_k.get_radio()
    vect_g = []
    index = np.array([i for i in Z_k], dtype=int)
    Dic_index = dict()
    for a in index:
        Dic_index[a] = g(a)
    vect_g = np.array(vect_g)
    matrix_g = np.zeros((p**k, p**k))
    num_row = 0
    for j in index:
        New_Values = (index - j) % p**k
        matrix_g[num_row] = np.array([Dic_index[a] for a in New_Values])
        num_row += 1
    return p**(-k)*matrix_g


def plot(A, Z_k, As="function",
         all_tree=True,
         size=(8, 8),
         s=1/5,
         m=0,
         size_points=0.5):
    """


    Parameters
    ----------
    A : TYPE
        DESCRIPTION.
    Z_k : TYPE
        DESCRIPTION.
    As : TYPE, (optional) Have to be "function",
                            "kernel", "convolution",
                            "fractal"
        DESCRIPTION. The default is "function".
    all_tree : TYPE, optional
        DESCRIPTION. The default is True.
    size : TYPE, optional
        DESCRIPTION. The default is (8, 8).
    s : TYPE, optional
        DESCRIPTION. The default is 1/5.
    m : TYPE, optional
        DESCRIPTION. The default is 0.
    size_points : TYPE, optional
        DESCRIPTION. The default is 0.5.

    Returns
    -------
    None.

    """
    if As == "function":
        """ Generate distance matrix. """
        K = Z_k.get_radio()
        p = Z_k.get_prime()
        Number = len(Z_k)
        Y = np.array([b for b in Z_k])
        D = np.zeros([Number, Number])
        C_A = np.zeros([Number, Number])
        if A == 0:
            for i in range(Number):
                for j in range(Number):
                    D[i, j] = norm_p(Y[i]-Y[j], p)
        else:
            for i in range(Number):
                for j in range(Number):
                    D[i, j] = norm_p(Y[i]-Y[j], p)
                    C_A[i, j] = A(Y[j])
        condensedD = squareform(D)

        """ Tree levels"""
        with_labels = True
        if K > 3:
            with_labels = False

        # labels
        if with_labels:
            if all_tree:
                level_max = K
            else:
                level_max = 0
            x_levels = []
            x_sticks = []
            for i in range(0, K):
                x_levels.append(p**(-i))
                x_sticks.append("Level " + str(i))
            x_levels.append(0)
            x_sticks.append("Level "+str(K))
        else:
            x_levels = []
            x_sticks = []
        """ y_axis"""
        if with_labels:
            Yticks = range(Number)
            Yticklabels = Y
            # position heat maps A and B
            first_map_position = [0.04, 0.1, 0.2, 0.6]
            second_map_position = [0.3, 0.74, 0.6, 0.2]
            # position heat maps X_0 and U
            position_tree_input = [0.3, 0.3, 0.6, 0.25]
        else:
            # labels
            Yticks = []
            Yticklabels = []
            # position heat maps A and B
            first_map_position = [0.1, 0.1, 0.2, 0.6]
            second_map_position = [0.3, 0.702, 0.6, 0.2]
            # position heat maps X_0 and U
            position_tree_input = [0.3, 0.25, 0.6, 0.25]
        """  Heat map of the A"""
        fig = pylab.figure(figsize=size)

        # Compute and plot second dendrogram.
        ax2 = fig.add_axes(position_tree_input)
        Y2 = sch.linkage(condensedD, method='single')
        Z2 = sch.dendrogram(Y2)
        ax2.set_yticks(x_levels)
        ax2.set_yticklabels(x_sticks, rotation='horizontal')
        ax2.set_xticklabels(Yticklabels, minor=False)
        plt.title("As Function")

        # Function A
        axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.15])
        idx2 = Z2['leaves']
        im = axmatrix.matshow(C_A, aspect='auto', origin='lower')
        axmatrix.set_xticks([])
        axmatrix.set_yticks([])

        # Plot colorbar.
        axcolor = fig.add_axes([0.98, 0.1, 0.05, 0.45])
        pylab.colorbar(im, cax=axcolor, format='%.4f')

        # Plot axes
        axmatrix.set_xticks(Yticks)
        axmatrix.set_xticklabels(Yticklabels, minor=False)
        axmatrix.xaxis.set_label_position('bottom')
        axmatrix.xaxis.tick_bottom()

        pylab.xticks(rotation=-90, fontsize=8)

    elif As == "kernel":
        """ Generate distance matrix. """
        K = Z_k.get_radio()
        p = Z_k.get_prime()
        Number = len(Z_k)
        Y = np.array([b for b in Z_k])
        D = np.zeros([Number, Number])
        """ Generate distance matrix. """
        if A == 0:
            C_A = np.zeros([Number, Number])
            for i in range(Number):
                for j in range(Number):
                    D[i, j] = norm_p(Y[i]-Y[j], p)
        else:
            C_A = p_adic_matrix(A, Z_k)
            for i in range(Number):
                for j in range(Number):
                    D[i, j] = norm_p(Y[i]-Y[j], p)
        condensedD = squareform(D)

        """ Tree levels"""
        with_labels = True
        if K > 3:
            with_labels = False

        # labels
        if with_labels:
            if all_tree:
                level_max = K
            else:
                level_max = 0
            x_levels = []
            x_sticks = []
            for i in range(0, K):
                x_levels.append(p**(-i))
                x_sticks.append("Level " + str(i))
            x_levels.append(0)
            x_sticks.append("Level "+str(K))
        else:
            x_levels = []
            x_sticks = []
        """ y_axis"""
        if with_labels:
            Yticks = range(Number)
            Yticklabels = Y
            # position heat maps A and B
            first_map_position = [0.04, 0.1, 0.2, 0.6]
            second_map_position = [0.3, 0.74, 0.6, 0.2]
            # position heat maps X_0 and U
            position_tree_input = [0.3, 0.3, 0.6, 0.25]
        else:
            # labels
            Yticks = []
            Yticklabels = []
            # position heat maps A and B
            first_map_position = [0.1, 0.1, 0.2, 0.6]
            second_map_position = [0.3, 0.702, 0.6, 0.2]
            # position heat maps X_0 and U
            position_tree_input = [0.3, 0.25, 0.6, 0.25]
        # Compute and plot first dendrogram.
        fig = pylab.figure(figsize=size)
        ax1 = fig.add_axes(first_map_position)
        Y1 = sch.linkage(condensedD, method='centroid')
        Z1 = sch.dendrogram(Y1, orientation='left')
        ax1.set_xticks(x_levels)
        ax1.set_xticklabels(x_sticks, rotation='vertical')
        ax1.set_yticks([])

        # Compute and plot second dendrogram.
        ax2 = fig.add_axes(second_map_position)
        Y2 = sch.linkage(condensedD, method='single')
        Z2 = sch.dendrogram(Y2)
        ax2.set_yticks(x_levels)
        ax2.set_yticklabels(x_sticks, rotation='horizontal')
        ax2.set_xticklabels(Yticklabels, minor=False)
        plt.title("As Kernel")

        # Plot distance matrix.
        axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.6])
        idx1 = Z1['leaves']
        idx2 = Z2['leaves']
        im = axmatrix.matshow(C_A, aspect='auto', origin='lower')
        axmatrix.set_xticks([])
        axmatrix.set_yticks([])

        # Plot colorbar.
        axcolor = fig.add_axes([0.98, 0.1, 0.05, 0.6])
        pylab.colorbar(im, cax=axcolor, format='%.4f')
        # fig.show()

        # Plot axes
        axmatrix.set_xticks(Yticks)
        axmatrix.set_xticklabels(Yticklabels, minor=False)
        axmatrix.xaxis.set_label_position('bottom')
        axmatrix.xaxis.tick_bottom()

        pylab.xticks(rotation=-90, fontsize=8)

        # axmatrix.set_yticks(range(40))
        axmatrix.set_yticks(Yticks)
        axmatrix.set_yticklabels(Yticklabels, minor=False)
        axmatrix.yaxis.set_label_position('right')
    elif As == "convolution":
        """ Generate distance matrix. """
        K = Z_k.get_radio()
        p = Z_k.get_prime()
        Number = len(Z_k)
        Y = np.array([b for b in Z_k])
        D = np.zeros([Number, Number])
        """ Generate distance matrix. """
        if A == 0:
            C_A = np.zeros([Number, Number])
            for i in range(Number):
                for j in range(Number):
                    D[i, j] = norm_p(Y[i]-Y[j], p)
        else:
            C_A = p_adic_Convolution_matrix(A, Z_k)
            for i in range(Number):
                for j in range(Number):
                    D[i, j] = norm_p(Y[i]-Y[j], p)
        condensedD = squareform(D)

        """ Tree levels"""
        with_labels = True
        if K > 3:
            with_labels = False

        # labels
        if with_labels:
            if all_tree:
                level_max = K
            else:
                level_max = 0
            x_levels = []
            x_sticks = []
            for i in range(0, K):
                x_levels.append(p**(-i))
                x_sticks.append("Level " + str(i))
            x_levels.append(0)
            x_sticks.append("Level "+str(K))
        else:
            x_levels = []
            x_sticks = []
        """ y_axis"""
        if with_labels:
            Yticks = range(Number)
            Yticklabels = Y
            # position heat maps A and B
            first_map_position = [0.04, 0.1, 0.2, 0.6]
            second_map_position = [0.3, 0.74, 0.6, 0.2]
            # position heat maps X_0 and U
            position_tree_input = [0.3, 0.3, 0.6, 0.25]
        else:
            # labels
            Yticks = []
            Yticklabels = []
            # position heat maps A and B
            first_map_position = [0.1, 0.1, 0.2, 0.6]
            second_map_position = [0.3, 0.702, 0.6, 0.2]
            # position heat maps X_0 and U
            position_tree_input = [0.3, 0.25, 0.6, 0.25]
        # Compute and plot first dendrogram.
        fig = pylab.figure(figsize=size)
        ax1 = fig.add_axes(first_map_position)
        Y1 = sch.linkage(condensedD, method='centroid')
        Z1 = sch.dendrogram(Y1, orientation='left')
        ax1.set_xticks(x_levels)
        ax1.set_xticklabels(x_sticks, rotation='vertical')
        ax1.set_yticks([])

        # Compute and plot second dendrogram.
        ax2 = fig.add_axes(second_map_position)
        Y2 = sch.linkage(condensedD, method='single')
        Z2 = sch.dendrogram(Y2)
        ax2.set_yticks(x_levels)
        ax2.set_yticklabels(x_sticks, rotation='horizontal')
        ax2.set_xticklabels(Yticklabels, minor=False)
        plt.title("As Convolution")

        # Plot distance matrix.
        axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.6])
        idx1 = Z1['leaves']
        idx2 = Z2['leaves']
        im = axmatrix.matshow(C_A, aspect='auto', origin='lower')
        axmatrix.set_xticks([])
        axmatrix.set_yticks([])

        # Plot colorbar.
        axcolor = fig.add_axes([0.98, 0.1, 0.05, 0.6])
        pylab.colorbar(im, cax=axcolor, format='%.4f')
        # fig.show()

        # Plot axes
        axmatrix.set_xticks(Yticks)
        axmatrix.set_xticklabels(Yticklabels, minor=False)
        axmatrix.xaxis.set_label_position('bottom')
        axmatrix.xaxis.tick_bottom()

        pylab.xticks(rotation=-90, fontsize=8)

        # axmatrix.set_yticks(range(40))
        axmatrix.set_yticks(Yticks)
        axmatrix.set_yticklabels(Yticklabels, minor=False)
        axmatrix.yaxis.set_label_position('right')
    elif As == "fractal":
        x_position, y_position = Christiakov_emmending(Z_k, m, s)
        if A == 0:
            A_list = np.zeros(len(Z_k))
        else:
            A_list = []
            for a in Z_k:
                A_list.append(A(a))
        # """Show function A"""
        fig_tree = pylab.figure(figsize=size)
        im = plt.scatter(x_position,
                         y_position,
                         c=A_list,
                         s=size_points)
        plt.title("As function")
        """ Plot colorbar. """
        axcolor = fig_tree.add_axes([0.95, 0.1, 0.05, 0.6])
        pylab.colorbar(im, cax=axcolor)
