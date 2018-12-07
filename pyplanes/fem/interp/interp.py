#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# interp.py
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


class InterpolationMethod(object):
    """Template for FEM interpolation methods.

    Attributes
    ----------
    mesh: Mesh-like instance
        The mesh on which to operate

    Parameters
    ----------
    mesh: Mesh-like instance
        The mesh on which to operate
    """

    def __init__(self, mesh):
        self._mesh = mesh

    @classmethod
    def get_linsys_size(self):
        """Compute the size of the maximal linear system for 1 variable"""
        pass

    def get_unit(self, length):
        """Integrate the shape functions on a length ``length``.

        To be overridden upon inheritance from an Interpolation Class

        Parameters
        ----------
        length: float
            length of the edge

        Returns
        -------
        I: float
            Integral value (real valued)
        """
        pass

    def get_matrices(self, nodes_coords, get_C=False):
        """Computes the elementary matrices over the element described by nodes

        Parameters
        ----------
        nodes_coords: array[Nx{2 or 3}]
            list of the nodes' coords
        get_C: bool, False
            Choose to return the crossed matrix C (unimplemented)

        Returns
        -------
        H, Q, [C]: np.array
            Elementary matrices, H with derivatives, Q without, C mixed (derived/underived)
        """
        pass
