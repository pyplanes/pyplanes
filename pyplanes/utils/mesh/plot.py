#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# plot.py
#
# This file is part of pyplanes, a software distributed under the MIT license.
# For any question, please contact mathieu@matael.org.
#
# Copyright (c) 2018 The pyplanes authors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#

import matplotlib.tri as mtri


def plot_triangular_mesh(ax, mesh, hide=None, labels=[]):
    """
    Plots a triangular mesh on a given axis

    Parameters
    ----------
    ax: plt.Axis
        axis on which to plot
    mesh: Mesh/MeshPart
        mesh-like object to plot
    hide: iterable, optional
        Can contain 'nodes', 'edges' or 'shapes' not to show the ids of this class
    labels: iterable
        Contain "nodes', 'edges' or 'shapes' to show the label corresponding to the object
        besides its id (not that if the object is in hide, it won't show here)
    """

    hide = hide if hide is not None else []

    triangles = mtri.Triangulation(
        mesh.nodes[:,0], mesh.nodes[:,1],
        mesh.shapes
    )
    ax.triplot(triangles, lw=.5)

    if not 'nodes' in hide:
        for i_n, n in enumerate(mesh.nodes):
            text = "{}({})".format(i_n, mesh.nodes_labels[i_n]) if 'nodes' in labels else str(i_n)
            ax.text(*n, text, ha='center', color='r')

    if not 'edges' in hide:
        for i_e, e in enumerate(mesh.edges):
            center = mesh.nodes[e,:].mean(axis=0)
            text = "{}({})".format(i_e, mesh.edges_labels[i_e]) if 'edges' in labels else str(i_e)
            ax.text(*center, text, ha='center', color='b')

    if not 'shapes' in hide:
        for i_s, s in enumerate(mesh.shapes):
            center = mesh.nodes[s,:].mean(axis=0)
            text = "{}({})".format(i_s, mesh.shapes_labels[i_s]) if 'shapes' in labels else str(i_s)
            ax.text(*center, text, ha='center', color='k')
