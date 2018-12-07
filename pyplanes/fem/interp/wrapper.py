#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# wrapper.py
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

class InterpWrapper(object):
    """Template for FEM interpolation methods.

    Attributes
    ----------
    __I: InterpolationMethod subclass
        Actual interpolator class (all methods defined in InterpolationMethod are mirrored)
    """

    def __init__(self, interp, *a, **kw):
        self.__I = interp(*a, **kw)

    def get_unit(self, *a, **kw):
        return self.__I.get_unit(*a, **kw)

    def get_matrices(self, *a, **kw):
        return self.__I.get_matrices(*a, **kw)

