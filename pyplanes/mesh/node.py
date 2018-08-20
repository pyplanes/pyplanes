#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# node.py
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


class Node(object):
    """ Stores a node's coordinates and labels

    coords -- coordinates (stored as a numpy array)
    labels -- list of labels
    dim -- number of coords (default: 2)
    """

    def __init__(self, coords, labels=None, dim=2):

        self.coords = np.array(coords)
        self.labels = [] if labels is None else labels
        self.dim = dim
