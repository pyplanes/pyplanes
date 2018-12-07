#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# 2d_solve.py
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

from pyplanes import Solver
from pyplanes.media import Air

@pytest.fixture
def mesh_filename():
    return 'test/files/square.msh'


def test_2d_linear(mesh_filename):
    labels = { 0: {'medium': 'air' }, }
    S = Solver(mesh_filename, labels)
    A, b, x = S.solve(10, output_matrices=True)
    assert A.size

