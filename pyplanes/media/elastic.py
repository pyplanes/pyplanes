#! /usr/bin/env python
# -*- coding:utf8 -*-
#
# elastic.py
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

from .medium import Medium


class Elastic(Medium):
    """Represent an elastic solid

    Attributes
    ----------
    E:
        Young's modulus
    nu:
        Poisson's Ratio
    rho:
        density
    eta:
        loss factor
    """

    MEDIUM_TYPE = 'elastic'
    MODEL = MEDIUM_TYPE
    EXPECTED_PARAMS = [
        ('E', float),  # Young's modulus
        ('nu', float),  # Poisson ratio
        ('rho', float),  # Density
        ('eta', float),  # loss factor
    ]
    OPT_PARAMS = [
        ('lambda_', complex),
        ('mu', complex)
    ]

    def __init__(self):
        super().__init__()

        self.E = None
        self.rho = None
        self.nu = None
        self.eta = None
        self.lambda_ = None
        self.mu = None

    def from_dict(self, *a, **kw):
        super().from_dict(*a, **kw)

    def _compute_missing(self):

        if self.lambda_ is None:
            self.lambda_ = (1+1j*self.eta)*(self.E*self.nu)/((1+self.nu)*(1-2*self.nu))
        if self.mu is None:
            self.mu = (1+1j*self.eta)*(self.E)/(2*(1+self.nu))

    def update_frequency(self, f):
        self.f = f
