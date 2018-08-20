#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# edge.py
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
from numpy.lib.scimath import sqrt


class Edge(object):
    """ Holds an edge (2 nodes) and its labels
    n1, n2 - connected nodes
    labels - list of labels
    """

    def __init__(self, n1, n2, labels=None):

        assert n1.dim == n2.dim, 'Both the Node objects passed to Edge must have the same dimension'

        self.n1 = n1
        self.n2 = n2
        self.labels = [] if labels is None else labels

    def length(self):
        return sqrt(self.n1.coords**2 - self.n2.coords**2)
