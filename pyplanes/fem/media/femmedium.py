#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# femmedium.py
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

from enum import IntEnum


class FEMMedium(object):
    """Holds the FEM representation of a given medium

    Attributes
    ----------
    fields: IntEnum
        all the fields considered to write the weak form (usually between 1 and 3)
    multipliers: dict
        for each field, the list of factors appearing in front of the matrices in the
        discretised weak form
    """

    def __init__(self, fields):
        self.f = -1
        self.multipliers = {}
        self.fields = IntEnum('fields', fields, start=0)

    def get_multipliers(self, field):
        """Return the multipliers to the elementary matrices for the medium"""
        if self.multipliers.get(field) is None:
            self.multipliers[field] = self.compute_multipliers(field)
        return self.multipliers.get(field)

    def compute_multipliers(self, field):
        """Compute the multipliers using the frequency specified in ``self.f``"""
        pass

    def update_frequency(self, f):
        """ Update the frequency and resets the multpliers

        Parameters
        ----------

        f: complex
            frequency
        """
        self.f = f
        self.multipliers = {}
