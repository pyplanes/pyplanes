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

    Parameters
    ----------

    nodes : list coords
        list of all nodes coords (index-identified)
    nodes_labels : list[str]
        list of labels for nodes (index-matched with the previous one)
    edges : list[list[int]]
        list of nodes linked by an edge (refering to nodes: id_n1 , id_n2)
    edges_labels : list[str]
        list of labels for edges (index-matched)
    shapes : list[list[int]]
        list of shapes (refering to nodes)
    shapes_labels : list[str]
        list of labels for shapes (index-matched)

    Methods
    -------
    iter_nodes(), iter_edges(), iter_shapes():
        return a iterator of all the nodes, edges or shapes in the mesh
    node(idx), edge(idx), shape(idx):
        return the node/edge/shape whose identifier is idx as a Node/Edge or Shape instance
    add_midpoints():
        Create a second-order mesh by adding a node at the middle of every edge.
    """

    def __init__(self,
                 nodes=None, nodes_labels=None,
                 edges=None, edges_labels=None,
                 shapes=None, shapes_labels=None, base_mesh=None):

        if base_mesh is not None:
            self.nodes = base_mesh.nodes
            self.nodes_labels = base_mesh.nodes_labels
            self.edges = base_mesh.edges
            self.edges_labels = base_mesh.edges_labels
            self.shapes = base_mesh.shapes
            self.shapes_labels = base_mesh.shapes_labels
            self.nb = base_mesh.nb

        else:
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

    def iter_nodes(self):
        return (_ for _ in self.nodes)

    def iter_edges(self):
        return (_ for _ in self.edges)

    def iter_shapes(self):
        return (_ for _ in self.shapes)

    def node(self, idx):
        """ Access primitive for a Node

        Parameters
        ----------
        idx : int
            ID of the node of interest

        Returns
        -------
        n : Node
            Correspnding node instance

        Raises
        ------
        IndexError
            if idx is out of range
        """
        return Node(self.nodes[idx], self.nodes_labels[idx])

    def edge(self, idx):
        """ Access primitive for an Edge

        Parameters
        ----------
        idx : int
            ID of the edge of interest

        Returns
        -------
        n : Edge
            Correspnding Edge instance

        Raises
        ------
        IndexError
            if idx is out of range
        """
        return Edge(list(map(Node, self.edges[idx])), self.edges_labels[idx])

    def shape(self, idx):
        """ Access primitive for a Shape

        Parameters
        ----------
        idx : int
            ID of the shape of interest

        Returns
        -------
        n : Shape
            Correspnding Shape instance

        Raises
        ------
        IndexError
            if idx is out of range
        """
        return Shape(list(map(Node, self.shapes[idx])), self.shapes_labels[idx])

    def add_midpoints(self):
        """Adds a node in the middle of each edge to make it a second-order mesh."""

        Nshapes, Nshapes_labels = [], self.shapes_labels[:]
        Nedges, Nedges_labels = [], []
        Nnodes, Nnodes_labels = self.nodes.tolist(), self.nodes_labels[:]

        # build a backref dict for edge labels
        backref_edges = {tuple(self.edges[i_e]): label for i_e, label in enumerate(self.edges_labels)}

        new_nodes_dict = {}
        last_node_id = len(Nnodes) - 1
        for shape_nodes in self.shapes:
            Nshapes.append([])
            for i_node, node in enumerate(shape_nodes):
                id_node2 = 0 if (i_node+1) == len(shape_nodes) else i_node+1
                n2 = (node, id_node2)
                n2_inv = tuple(n2[::-1])

                if new_nodes_dict.get(n2) is None and new_nodes_dict.get(n2_inv) is None:
                    Nnodes.append(self.nodes[n2,:].mean(axis=0).tolist())
                    Nnodes_labels.append(0)
                    last_node_id += 1
                    mid_id = last_node_id

                    new_nodes_dict[n2] = mid_id
                    new_nodes_dict[n2_inv] = mid_id

                    # get the new edges label
                    if backref_edges.get(n2) is not None:
                        edge_label = backref_edges[n2]
                    elif backref_edges.get(n2_inv) is not None:
                        edge_label = backref_edges[n2]
                    else:
                        edge_label = [0]

                    # register the new edges
                    Nedges.append([n2[0], mid_id])
                    Nedges.append([mid_id, n2[1]])
                    Nedges_labels += [edge_label, edge_label]
                else:
                    mid_id = new_nodes_dict.get(n2, new_nodes_dict[n2_inv])

                Nshapes[-1] += [node, mid_id]

        Nnodes = np.array(Nnodes)
        return Mesh(Nnodes, Nnodes_labels, Nedges, Nedges_labels, Nshapes, Nshapes_labels)
