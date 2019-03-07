#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# triangles.py
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

from .interp import InterpolationMethod


class TR3(InterpolationMethod):
    """First-Order Polynomial Interpolation on 3-points triangles"""

    def __init__(self, mesh):
        super().__init__(mesh)

    @classmethod
    def get_linsys_size(self, mesh):
        return len(mesh.nodes)

    def get_unit(self, length):
        """Integrate the shape functions on a length ``length``.

        Parameters
        ----------
        length: float
            length of the edge

        Returns
        -------
        I: float
            Integral value (real valued)
        """

        _PG1 = 0.339981043584856
        _PG2 = 0.861136311594053
        _WG1 = 0.652145154862546
        _WG2 = 0.347854845137454
        gauss_points = (
            (_PG1, _WG1),
            (-_PG1, _WG1),
            (_PG2, _WG2),
            (-_PG2, _WG2),
        )

        I = np.zeros((2,))
        for xi, w in gauss_points:
            I[0] += (1+xi/2)*w
            I[1] += (1-xi/2)*w
        I *= length/2
        return I

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
        dN = np.array([
            [-1, 1, 0],
            [-1, 0, 1]])
        J = dN @ nodes_coords
        Jinv = np.linalg.inv(J)
        Jdet = np.linalg.det(J)

        _PG1 = 0.445948490915965
        _PG2 = 0.091576213509771
        _WG1 = 0.111690794839005
        _WG2 = 0.054975871827661
        gauss_points = (
            (1-2*_PG1, _PG1, _WG1),
            (_PG1, _PG1, _WG1),
            (_PG1, 1-2*_PG1, _WG1),
            (_PG2, _PG2, _WG2),
            (1-2*_PG2, _PG2, _WG2),
            (_PG2, 1-2*_PG2, _WG2)
        )

        H, Q = np.zeros((3,3)), np.zeros((3,3))
        for xi, eta, w in gauss_points:
            N = np.array([1-xi-eta, xi, eta]).reshape(1,3)
            Bd = Jinv @ dN
            H += Bd.T @ Bd * Jdet * w
            Q += N.T @ N * Jdet * w
        return H, Q


class TR6(InterpolationMethod):
    """Second-Order Polynomial Interpolation on 6-points triangles"""

    def __init__(self, mesh):
        super().__init__(mesh)

    @classmethod
    def get_linsys_size(self, mesh):
        return len(mesh.nodes)

    def get_unit(self, length):
        """Integrate the shape functions on a length ``length``.

        Parameters
        ----------
        length: float
            length of the edge

        Returns
        -------
        I: float
            Integral value (real valued)
        """

        _PG1 = 0.339981043584856
        _PG2 = 0.861136311594053
        _WG1 = 0.652145154862546
        _WG2 = 0.347854845137454
        gauss_points = (
            (_PG1, _WG1),
            (-_PG1, _WG1),
            (_PG2, _WG2),
            (-_PG2, _WG2),
        )

        I = np.zeros((3,))
        for xi, w in gauss_points:
            I[0] += (-(1-xi)*xi/2) * w
            I[1] += ((1+xi)*xi/2) * w
            I[2] += ((1-xi)*(1+xi)) * w
        I *= length/2
        return I

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
        _PG1 = 0.445948490915965
        _PG2 = 0.091576213509771
        _WG1 = 0.111690794839005
        _WG2 = 0.054975871827661
        gauss_points = (
            (1-2*_PG1, _PG1, _WG1),
            (_PG1, _PG1, _WG1),
            (_PG1, 1-2*_PG1, _WG1),
            (_PG2, _PG2, _WG2),
            (1-2*_PG2, _PG2, _WG2),
            (_PG2, 1-2*_PG2, _WG2)
        )

        H, Q = np.zeros((6,6)), np.zeros((6,6))
        for xi, eta, w in gauss_points:
            zeta = 1 - xi - eta

            N = np.array([
                -zeta*(1-2*zeta),
                4*xi*zeta,
                -xi*(1-2*xi),
                4*xi*eta,
                -eta*(1-2*eta),
                4*eta*zeta]).reshape(1,3)

            dN = np.array([
                [1-4*zeta, 4*(zeta-xi), -1+4*xi, 4*eta, 0, -4*eta],
                [1-4*zeta, -4*xi, 0, 4*xi, -1+4*eta, 4*(zeta-eta)],
            ])
            J = dN @ node_coords
            Jinv = np.linalg.inv(J)
            Jdet = np.linalg.det(J)

            Bd = Jinv @ dN
            H += Bd.T @ Bd * Jdet * w
            Q += N.T @ N * Jdet * w
        return H, Q

