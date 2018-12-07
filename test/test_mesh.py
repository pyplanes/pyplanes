#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# test_mesh.py
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
from pyplanes.mesh import MeshPart


def test_wrong_limits():
    full_mesh = freefem_reader('test/files/test_parts.msh')
    with pytest.raises(ValueError):
        part = MeshPart(full_mesh, weird_name=None)


def test_meshparts():
    full_mesh = freefem_reader('test/files/test_parts.msh')
    part = MeshPart(full_mesh, labels=(1,))

    for i_s in [4,5,6,7]:
        assert (full_mesh.shapes[i_s] == part.shapes[i_s-4]).all()

    for ii, i_e in enumerate([3,4,5,6,9,10,11,14,15]):
        assert (full_mesh.edges[i_e]==part.edges[ii]).all()

    nodes = list(part.iter_nodes())
    for ii, i_n in enumerate([3,4,5,6,7,8]):
        assert (full_mesh.nodes[i_n]==nodes[ii]).all()

