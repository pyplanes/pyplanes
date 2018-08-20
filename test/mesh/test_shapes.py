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

import unittest
import numpy as np

from pyplanes.mesh import Shape, Node


class ShapeTestCase(unittest.TestCase):

    def setUp(self):
        self.nodes = [
            Node(np.array([0,0])),
            Node(np.array([1,0])),
            Node(np.array([0,1]))
        ]

    def test_sanity(self):
        """Tests that it's not possible to define a 2 nodes shape"""

        self.assertRaises(ValueError, Shape, self.nodes[:-1])

    def test_instanciation(self):
        s = Shape(self.nodes)
        for i in range(3):
            self.assertEqual(self.nodes[i], s.nodes[i])

    def test_normal_positive(self):
        s = Shape(self.nodes)
        self.assertEqual(s.compute_normal(), np.float(1))
        self.assertEqual(s.normal, np.float(1))

    def test_normal_negative(self):
        s = Shape(self.nodes[::-1])
        self.assertEqual(s.compute_normal(), np.float(-1))
        self.assertEqual(s.normal, np.float(-1))
