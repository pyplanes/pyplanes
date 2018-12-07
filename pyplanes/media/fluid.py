#! /usr/bin/env python
# -*- coding:utf8 -*-
#
# fluid.py
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
# This file is part of pymls, a software distributed under the MIT license.
# For any question, please contact one of the authors cited below.
#

from .medium import Medium


class Fluid(Medium):
    """ Represent a fluid medium

    Attributes
    ----------

    rho : float
        density
    c : float
        sound speed
    """

    MEDIUM_TYPE = 'fluid'
    MODEL = MEDIUM_TYPE
    EXPECTED_FIELDS = [
        ('rho', float),  # Density
        ('c', float),  # Sound speed
    ]

    def __init__(self):
        """ Nothing special todo for fluid """
        super().__init__()

        self.rho = None
        self.c = None

    def update_frequency(self, f):
        """ For a fluid, does nothing."""
        super().update_frequency(f)
