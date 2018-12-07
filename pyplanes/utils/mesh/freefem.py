#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# freefem.py
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

import numpy as np

from pyplanes.mesh import Mesh


def freefem_reader(freefem_mshfile):
    """
    FreeFEM++ mesh reader (return a Mesh instance)

    If the file can't be read, the exception is caught and thrown anew by the function

    Parameters
    ----------
    freefem_mshfile: str
        path to file
    """

    try:
        fh = open(freefem_mshfile)
    except FileNotFoundError:
        raise FileNotFoundError('The mesh file {} could not be found.'.format(freefem_mshfile))
    nb_nodes, nb_shapes, nb_edges = list(map(int, fh.readline().split(' ')))

    i, line_len  = 0, None
    nodes, nodes_labels = [], []
    while i < nb_nodes:
        split_line = fh.readline().split(' ')

        if line_len is None:
            line_len = len(split_line)
        elif line_len != len(split_line):
            raise ValueError('Number of coordinates differ between nodes (had seen {}, now {})'.format(line_len, len(split_line)))
        # The last element is always a label (0 if unset)
        nodes.append(np.array(split_line[:-1], dtype=np.float))
        nodes_labels.append(int(split_line[-1]))
        i += 1

    i, line_len = 0, None
    shapes, shapes_labels = [], []
    while i < nb_shapes:
        split_line = fh.readline().split(' ')

        if line_len is None:
            line_len = len(split_line)
        elif line_len != len(split_line):
            raise ValueError('Number of nodes differ between shapes (had seen {}, now {})'.format(line_len, len(split_line)))
        # The last element is always a label (0 if unset)
        # The nodes are indexed from 1 on by FF++, so we offset by -1 here
        shapes.append(np.array(split_line[:-1], dtype=np.float)-1)
        shapes_labels.append(int(split_line[-1]))
        i += 1

    i, line_len = 0, None
    edges, edges_labels = [], []
    while i < nb_edges:
        split_line = fh.readline().split(' ')

        if line_len is None:
            line_len = len(split_line)
        elif line_len != len(split_line):
            raise ValueError('Number of nodes differ between shapes (had seen {}, now {})'.format(line_len, len(split_line)))
        # The last element is always a label (0 if unset)
        # The nodes are indexed from 1 on by FF++, so we offset by -1 here
        edges.append(np.array(split_line[:-1], dtype=np.float)-1)
        edges_labels.append(int(split_line[-1]))
        i += 1

    fh.close()
    return Mesh(np.array(nodes), nodes_labels, np.array(edges, dtype=int), edges_labels, np.array(shapes, dtype=int), shapes_labels)
