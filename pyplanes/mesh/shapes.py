#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# shape.py
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


class Shape(object):
    """
    Holds a shape

    Parameters
    ----------
    nodes: node list[len>2]
        list of nodes. must be correctly winded, it *won't* be checked
    labels: list
        list of labels
    """

    def __init__(self, nodes, labels=None):

        if len(nodes) < 3:
            raise ValueError("Less than 3 nodes don't make a surface")

        self.nodes = nodes
        self.coords = np.array([_.coords for _ in nodes])
        self.labels = [] if labels is None else labels

        self._normal = None
        self._centre = None

    @property
    def normal(self):
        """
        Compute and store a normalized surface normal

        Return
        ------
        normal: -1 or 1
            the shape's normal
        """

        if self._normal is None:
            raw_normal = np.cross(
                self.coords[0] - self.coords[1],
                self.coords[0] - self.coords[2]
            )
            self._normal = raw_normal/np.abs(raw_normal)
        return self._normal

    @property
    def centre(self):
        """
        Compute and store the centre of the shape

        Defaults to an average of the coordinates.

        Returns
        -------

        centre: array of floats
            Position of the centre
        """
        if self._centre is None:
            self._centre = self.coords.mean(axis=0)
        return self._centre
