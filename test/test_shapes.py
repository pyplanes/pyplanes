#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# shapes.py
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

import pytest
import numpy as np

from pyplanes.mesh import Shape, Node


@pytest.fixture
def nodes():
    return [
        Node(np.array([0,0])),
        Node(np.array([1,0])),
        Node(np.array([0,1]))
    ]


def test_sanity(nodes):
    """Tests that it's not possible to define a 2 nodes shape"""

    with pytest.raises(ValueError):
        Shape(nodes[:-1])


def test_instanciation(nodes):
    s = Shape(nodes)
    for i in range(3):
        assert nodes[i] == s.nodes[i]
    assert (np.array([[0,0],[1,0],[0,1]])==s.coords).all()


def test_normal_positive(nodes):
    s = Shape(nodes)
    assert s.normal == np.float(1)


def test_normal_negative(nodes):
    s = Shape(nodes[::-1])
    assert s.normal == np.float(-1)


def test_centre(nodes):
    centre = Shape(nodes).centre
    assert (np.array([1/3, 1/3])==centre).all()


def test_addmidpoint():
    nodes = [(0,0), (1,0), (1,0), (1,1)]
    nodes_labels = [[0], [0], [0], [0]]
    shapes = [[1,2,4], [1,4,3]]
    shapes_labels = [[0], [0]]
    edges = [[1,2], [2,4], [4,3], [3,1]]
    edges_labels = [[1], [2], [3], [4]]

