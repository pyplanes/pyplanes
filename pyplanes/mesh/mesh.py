#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# mesh.py
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

from pyplanes.mesh import Node, Edge, Shape


class Mesh(object):
    """ Holds a mesh and provides access primitives for nodes, edges and shapes
    The mesh is defined as lists and can be used either as a storage (through [] access)
    or an object provider (through <object_name>() methods).

    nodes - list of nodes  coords (index-identified)
    nodes_labels - list of labels for nodes (index-matched with the previous one)
    edges - list of edges (refering to nodes)
    edges_labels - list of labels for edges (index-matched)
    shapes - list of shapes (refering to nodes)
    shapes_labels - list of labels for shapes (index-matched)
    """

    def __init__(self,
                 nodes=None, nodes_labels=None,
                 edges=None, edges_labels=None,
                 shapes=None, shapes_labels=None):

        if nodes_labels is not None and not len(nodes) == len(nodes_labels):
            raise ValueError('Nodes and nodes labels lists must have the same length')
        if edges_labels is not None and not len(edges) == len(edges_labels):
            raise ValueError('Edges and edges labels lists must have the same length')
        if shapes_labels is not None and not len(shapes) == len(shapes_labels):
            raise ValueError('Shapes and shapes labels lists must have the same length')

        self.nodes = nodes
        self.nodes_labels = nodes_labels
        self.edges = edges
        self.edges_labels = edges_labels
        self.shapes = shapes
        self.shapes_labels = shapes_labels
        self.nb = {
            'nodes': len(self.nodes),
            'edges': len(self.edges),
            'shapes': len(self.shapes)
        }

    def node(self, idx):
        return Node(self.nodes[idx], self.nodes_labels[idx])

    def edge(self, idx):
        return Edge(list(map(Node, self.edges[idx])), self.edges_labels[idx])

    def shape(self, idx):
        return Shape(list(map(Node, self.shapes[idx])), self.shapes_labels[idx])
