#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# nodelistpart.py
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


class NodeListPart(list):
    """
    Proxy class to a list to handle a node list with consistent IDs for a pruned mesh

    Iterates only on the valid nodes but allows to ``__getitem__`` on the full list so the
    IDs are correct.

    Attributes
    ----------
    all_nodes: array
        complete list of nodes as specified in a Mesh
    valid_nodes_idxs: array
        all IDs of valid nodes

    Parameters
    ----------
    all_nodes: array
        complete list of nodes as specified in a Mesh
    valid_nodes_idxs: array
        all IDs of valid nodes
    """

    def __init__(self, all_nodes, valid_nodes_idxs):
        self.all_nodes = all_nodes
        self.valid_nodes_idxs = valid_nodes_idxs
        self.__nb_valid = len(self.valid_nodes_idxs)
        self.__n = 0

    def __len__(self):
        return self.__nb_valid

    def __iter__(self):
        self.__n = 0
        return self

    def __next__(self):
        if self.__n < self.__nb_valid:
            self.__n += 1
            return self.all_nodes[self.valid_nodes_idxs[self.__n-1]]
        else:
            raise StopIteration

    def __getitem__(self, *a, **kw):
        return self.all_nodes.__getitem__(*a, **kw)

    def __setitem__(self, *a, **kw):
        return self.all_nodes.__setitem__(*a, **kw)
