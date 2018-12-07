#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# fem.py
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

from pyplanes.core.method import Method
from pyplanes.fem.interp import InterpWrapper
from pyplanes.fem.meshutils import FEMMesh
from pyplanes.mesh import Mesh, MeshPart


class FEMDomain(Method, InterpWrapper):
    """
    Represents a FEM modelled domain, with its part of Mesh

    FEM domains link a Mesh, a Medium, and InterpolationMethod and a local linear system.

    Attributes
    ----------
    mesh: FEMMesh
        mesh of the domain (wrapper around a ``pyplanes.mesh.MeshPart`` object)
    medium: FEMMedium
        FEM aware medium class that spreads the domain
    nb_dof: int
        number of degrees of freedom of the underlying linear system
    node_map: array
        maps a node number from the global mesh to the local linear system. An offset must
        be applied to get the other fields (check ``FEMDomain.to_dof()``)
    bc_edges: list
        list of edges on which a non trivial boundary condition is applied

    Parameters
    ----------
    interp: InterpolationMethod's subclass
        interpolation method to be used in the domain
    mesh: Mesh/MeshPart
        mesh on which the domain spans
    media: dict
        global media dictionnary (to get an update medium everytime)
    labels: dict
        global labels definition dictionnary

    Methods
    -------
    to_dof(self, node, field)
        Convert a mesh node ID into a position in the linear system. The available
        ``fields`` are found from ``self.medium``
    _identify_medium(self)
        Parses the mesh and its label to identify the medium filling the current domain
    _identify_bc(self)
        Parses the edges and labels to identify the one requiring attention
    assemble(self, f)
        Assemble the matrices of the domain for the frequency f. The result is a ijv
        representation
    apply_boundary_conditions(self, f)
        Goes over the boundary conditions to create a forcing vector or mark lines for removal
    """

    def __init__(self, interp, mesh, media, labels):
        Method.__init__(self)
        if isinstance(mesh, Mesh):
            self.mesh = FEMMesh(MeshPart(mesh, whole_mesh=True))
        else:
            self.mesh = FEMMesh(mesh)
        InterpWrapper.__init__(self, interp, self.mesh)

        self._media = media
        self._labels = labels

        self.medium = self._identify_medium()

        self.nb_dof = len(self.medium.fields)*interp.get_linsys_size(self.mesh)

        self.nodes_map = np.zeros((len(self.mesh.nodes.all_nodes),), dtype=int)
        for i_n, n in enumerate(self.mesh.valid_nodes_idxs):
            self.nodes_map[n] = i_n

        self.bc_edges = self._identify_bc()

    def to_dof(self, node, field):
        """
        Convert a mesh node ID into a position in the domain's linear system

        Parameters
        ----------
        node: int
            mesh node ID
        field: str, Enum Item
            field of interest (for multiple fields domains mainly, required.)

        Returns
        -------
        pos: int
            position in the linear system
        """

        return self.medium.fields[field]*self.nb_dof+self.nodes_map[node]

    def _identify_medium(self):
        """
        Parses the mesh and its label to identify the medium filling the current domain

        Returns
        -------
        m: FEMMedium
            the medium spanning the entire domain

        Raises
        ------
        ValueError:
            if there's mode than one medium on the domain or if the medium's label doesn't
            correspond to any medium in the global list.
        """

        labels_list = self.mesh.nodes_labels + self.mesh.edges_labels + self.mesh.shapes_labels
        media_seen = set([_ for _ in labels_list if _ in self._labels.keys() and self._labels[_].get('medium') is not None])
        if len(media_seen) != 1:
            raise ValueError('Unable to process a domain with {} media'.format(len(media_seen)))
        else:
            return self._media[list(media_seen)[0]]

    def _identify_bc(self):
        """
        Parses the edges and labels to identify the one requiring attention

        This method doesn't retain interface edges : these are treated at the global level

        Returns
        -------
        edges: list
            all edges to be considered when applying boundary conditions

        Notes
        -----
        Placeholder implementation for now, requires modification.
        """

        return list(range(len(self.mesh.edges)))

    def assemble(self, f):
        """
        Assemble the matrices of the domain for a given frequency.

        The result is a ijv representation returned as a tuple.

        Parameters
        ----------
        f: complex
            frequency of the analysis

        Returns
        -------
        Ai, Aj, Av:
            arrays defining the local linear system's sparse matrix

        Notes
        -----
        This method stores the result of the computation in the instance as well with 3
        attributes ``Ai``, ``Aj``, ``Av``
        """

        # initialise the matrices
        Ai, Aj, Av = [], [], []

        for elem_id, node_ids in enumerate(self.mesh.shapes):
            coords = self.mesh.nodes[node_ids,:]
            dof_ids = self.to_dof(node_ids, 'p')

            Hmul, Qmul = self.medium.get_multipliers('p')
            H, Q = self.get_matrices(coords)

            element_matrix = Hmul*H + Qmul*Q

            idx_h, idx_v = np.meshgrid(dof_ids, dof_ids, indexing='ij')
            Ai.append(idx_h.reshape(idx_h.size))
            Aj.append(idx_v.reshape(idx_h.size))
            Av.append(element_matrix.reshape(idx_h.size))

        self.Ai = np.concatenate(Ai)
        self.Aj = np.concatenate(Aj)
        self.Av = np.concatenate(Av)
        return self.Ai, self.Aj, self.Av

    def apply_boundary_conditions(self, f):
        """
        Goes over the boundary conditions to create a forcing vector or mark lines for removal

        Outputs a i(j)v-like vectorset. The value vector might contain ``None`` if a
        specific dof must be removed.

        Parameters
        ----------
        f: complex
            frequency of the analysis

        Returns
        -------
        bi, bv: array
            sparse-ready representation of the local forcing vector. Might contain
            ``None`` in the value array to mark a dof for removal.

        Notes
        -----
        For now, this method uses hard coded labels and **this must be changed**.
        """

        omega = 2*np.pi*f

        # put v=None to *remove the dof*
        bi, bv = [], []

        # TODO/ DON'T USE HARD CODED LABELS. EVER.
        for i_e in self.bc_edges:
            label = self.mesh.edges_labels[i_e]
            nodes = self.mesh.nodes[self.mesh.edges[i_e]]
            edge_length = np.linalg.norm(nodes[0] - nodes[1], 2)
            if label == 1:
                dof_ids = self.to_dof(self.mesh.edges[i_e], 'p')
                bi.append(dof_ids)
                bv.append(self.get_unit(edge_length)/(1j*omega))

        self.bi = np.concatenate(bi)
        self.bv = np.concatenate(bv)
        return self.bi, self.bv
