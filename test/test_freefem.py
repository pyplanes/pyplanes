#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# test_freefem.py
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

from pyplanes.utils.mesh import freefem_reader


def test_unreadablefile():

    inexistent_path = 'files/inexistent'
    with pytest.raises(FileNotFoundError):
        freefem_reader(inexistent_path)


def test_reading():
    mesh = freefem_reader('test/files/test.msh')
    expected_nodes = list(map(np.array, [
        [0.0, 0.0], [0.5, 0.0], [1.0, 0.0],
        [0.0, 0.5], [0.5, 0.5], [1.0, 0.5],
        [0.0, 1.0], [0.5, 1.0], [1.0, 1.0]
    ]))
    expected_nodes_labels = [1, 0, 0, 0, 1, 0, 0, 0, 1]

    for i in range(len(expected_nodes)):
        assert (mesh.nodes[i] == expected_nodes[i]).all()
    assert mesh.nodes_labels == expected_nodes_labels
    assert mesh.nb['nodes'] == 9
    assert mesh.nb['shapes'] == 8
    assert mesh.nb['edges'] == 16
