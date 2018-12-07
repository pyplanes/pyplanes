#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# meshpart.py
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
from functools import reduce

from .mesh import Mesh
from .nodelistpart import NodeListPart


class MeshPart(Mesh):
    """
    Holds a part of mesh, initialised with the complete mesh and a separation criterion

    Attributes
    ----------
    _full_mesh: Mesh-like
        original mesh-like instance
    _limits: dict
        criteria given to separated the mesh

    Parameters
    ----------
    base_mesh: Mesh-like
        mesh to use as a starting point
    **limits:
        set of criteria to use for separation
    """

    def __init__(self, base_mesh, **limits):
        if isinstance(base_mesh, MeshPart):
            super().__init__(base_mesh=base_mesh._full_mesh)
            self._full_mesh = base_mesh._full_mesh
            self._limits = base_mesh._limits
        else:
            super().__init__(base_mesh=base_mesh)
            self._full_mesh = base_mesh
            self._limits = limits
        self.__build_lists(self._limits)

    def __build_lists(self, limits):
        """Builds the lists of nodes/edges and shapes on which to iterate

        This method operates in place and modifies the attributes of the instance

        Parameters
        ----------
        limits: dict
            criteria used for separation (see Notes)

        Raises
        ------
        ValueError:
            if the criteria proposed are invalid

        Notes
        -----

        Criteria that can be specified are:
        labels=list(int)
            only the shapes whose label is in the criterion's list will be kept
        whole_mesh=True
            keep to full mesh (wrapper mode)
        """

        limits_types = list(limits.keys())
        if 'labels' in limits_types:
            labels = list(limits['labels'])
            shapes_idxs = [i_s for i_s, s in enumerate(self.shapes) if self.shapes_labels[i_s] in labels]
            nodes_idxs, edges_idxs = self.__find_nodes_edges(shapes_idxs)
            self.__filter_lists(nodes_idxs, edges_idxs, shapes_idxs)
        elif 'whole_mesh' in limits_types and limits['whole_mesh']:
            self.valid_nodes_idxs = list(range(len(self.nodes)))
            self.all_nodes = self.nodes
            self.nodes = NodeListPart(self.all_nodes, self.valid_nodes_idxs)
        else:
            raise ValueError('Unknown limits: {}'.format(limits))

    def __filter_lists(self, nodes_idxs, edges_idxs, shapes_idxs):
        """
        Given a list of nodes, edges and shapes ids, prunes the lists held by the instance
        in place

        Update the counts in ``self.nb`` as well.

        Parameters
        ----------
        nodes_idxs, edges_idxs, shapes_idxs: list
            lists of ids for the nodes, edges and shapes
        """
        self.shapes = [self._full_mesh.shapes[_] for _ in shapes_idxs]
        self.shapes_labels = [self._full_mesh.shapes_labels[_] for _ in shapes_idxs]
        self.edges = [self._full_mesh.edges[_] for _ in edges_idxs]
        self.edges_labels = [self._full_mesh.edges_labels[_] for _ in edges_idxs]

        self.all_nodes, self.all_nodes_labels = self.nodes, self.nodes_labels
        self.valid_nodes_idxs = nodes_idxs
        self.nodes = NodeListPart(self.all_nodes, self.valid_nodes_idxs)
        self.nodes_labels = NodeListPart(self.all_nodes_labels, self.valid_nodes_idxs)

        self.nb = {
            'nodes': len(self.valid_nodes_idxs),
            'edges': len(self.edges),
            'shapes': len(self.shapes)
        }

    def __find_nodes_edges(self, shapes_idxs):
        """
        Given a set of shapes ids, find all connected nodes and edges

        Parameters
        ----------
        shapes_idxs: list
            ids of the selected shapes

        Returns
        -------
        nodes_idxs, edges_idxs: tuple(list)
            returns the tuple of lists of ids for nodes and edges
        """

        nodes_idxs = list(reduce(
            lambda p, c: p.union(set(c)),
            [self._full_mesh.shapes[_] for _ in shapes_idxs],
            set()))
        nodes_boollist = np.zeros(len(self._full_mesh.nodes,), dtype=bool)
        for idx in nodes_idxs:
            nodes_boollist[idx] = True
        edges_idxs = [i_e for i_e, e in enumerate(self._full_mesh.edges) if nodes_boollist[e].all()]

        return nodes_idxs, edges_idxs
