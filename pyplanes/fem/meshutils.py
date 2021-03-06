#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# meshutils.py
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

from pyplanes.mesh import MeshPart


class FEMMesh(MeshPart):
    """Holds a mesh and utility methods for querying FEM-aware features

    Parameters
    ----------
    base_mesh: pyplanes.mesh.Mesh
        mesh instance upon which the FEMMesh is based
    """

    def __init__(self, base_mesh):
        super().__init__(base_mesh=base_mesh)
