#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# air.py
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

from pyplanes.media import Air as MediaAir
from .fluid import Fluid


class Air(Fluid, MediaAir):

    def __init__(self):
        Fluid.__init__(self)
        MediaAir.__init__(self)
        self.rho = MediaAir.rho
        self.c = MediaAir.c

    def update_frequency(self, f):
        Fluid.update_frequency(self, f)
