#! /usr/bin/env python3
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


import numpy as np

from .femmedium import FEMMedium
from pyplanes.media import Fluid as MediaFluid


class Fluid(FEMMedium, MediaFluid):
    """FEM representation of a fluid medium

    *Include weak-form*
    """

    def __init__(self):
        self.rho, self.c = None, None
        FEMMedium.__init__(self, 'p')
        MediaFluid.__init__(self)

    def update_frequency(self, f):
        FEMMedium.update_frequency(self, f)
        MediaFluid.update_frequency(self, f)

    def compute_multipliers(self, field):
        omega = 2*np.pi*self.f

        return 1/(self.rho*omega**2),\
            -1/(self.rho*self.c**2)
