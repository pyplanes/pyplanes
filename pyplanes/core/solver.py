#! /usr/bin/env python3
# -*- coding:utf8 -*-
#
# solver.py
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

import os
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from pyplanes.mesh import MeshPart
from pyplanes.fem.media import Air, from_yaml
from pyplanes.fem.interp import TR3
from pyplanes.utils.mesh import freefem_reader
from pyplanes.fem import FEMDomain


class Solver:
    """
    Solver class, orchestrates the creation and solving of the linear system

    Attributes
    ----------
    meshfile: str
        path to the mesh definition file
    mesh: pyplanes.mesh.Mesh
        mesh instance based on the provided mesh file
    labels: dict
        hashmap binding the mesh labels to their meaning, for elements and nodes
        holds the media used in the resolution. The keys are the labels indexes in which
        media are used.
    domains: list
        list of identified domains (part of the mesh with a method and a medium)
    domains_dof: list(tuple)
        for each domain of ``domains`` contains a 2-uple (first dof, number of dofs)
        allowing to identify which zone of the linear system belong to which method
    max_nb_dof: int
        maximal size of the linear system
    interfaces: list
        list of the identified inter-method interfaces

    Parameters
    ----------
    meshfile: str
        location of the meshfile to use. If the file doesn't exit, an IOError with be thrown
    labels: dict
        labels definition

    Methods
    -------
    solve(self, f, output_matrices=False)
        solve the linear system and (depending on ``output_matrices``) exports the
        matrices
    _identify_domains()
        used to populate the ``domains`` attributes. Based on the labels and the mesh,
        split the problems into subdomains
    _identify_media()
        extract all media paths from the ``labels`` attribute and load them

    Notes
    -----
    The current implementation automatically read the meshfile as a FreeFEM++ mesh.

    """

    def __init__(self, meshfile, labels):
        if not os.path.exists(meshfile):
            raise IOError('Unable to locate file {}'.format(meshfile))

        self.meshfile = meshfile
        self.mesh = freefem_reader(meshfile)

        self.labels = labels
        self.media = self._identify_media()

        self.domains = self._identify_domains()
        dof_per_domain = [_.nb_dof for _ in self.domains]
        self.domains_dof = [(sum(dof_per_domain[:i_nb]), nb) for i_nb, nb in enumerate(dof_per_domain)]

        # compute the maximum size of the global matrix
        self.max_nb_dof = sum(dof_per_domain)

    def _identify_media(self):
        """Parses the labels hashmap to find all media used and load them

        Returns
        -------
        media_dict: dict
            a {label: Medium subclass' instance} dictionnary

        Notes
        -----

        The string 'air' can be used to use the embedded air medium
        """

        mediapaths = {k: v['medium'] for k, v in self.labels.items() if v.get('medium') is not None}

        media_dict = {}
        for label, path in mediapaths.items():
            if path.lower() == 'air':
                media_dict[label] = Air()
            else:
                media_dict[label] = from_yaml(path)
        return media_dict

    def _identify_domains(self):
        """Parses the labels list and the mesh to identify domains

        Returns
        -------
        domains: list
            A list of domains

        Note
        ----

        **Must be reimplemented**
        """

        domains = [FEMDomain(TR3, MeshPart(self.mesh, labels=(0,)), self.media, self.labels)]
        return domains

    def __apply_dof_remove(self, belem_remove_idxs):
        pass

    def solve(self, f, output_matrices=False):
        """Assembles and solves the linear system by combining the domains at a given
        frequency f.

        Parameters
        ----------
        f: complex
            frequency of the analysis
        output_matrices: bool (optional)
            If true, the matrix and forcing vector of the linear system are returned along
            with the solution.

        Returns
        -------
        x: np.array
            solution vector
        A, b, x: scipy.sparse.coo_matrix, scipy.sparse.coo_matrix, np.array
            The matrices from the Ax=b linear system if the parameter ``output_matrices``
            is set to True.
        """

        # update the frequency in media
        for m in self.media.values():
            m.update_frequency(f)

        Ai, Aj, Av = [], [], []
        bi, bv = [], []
        for i_d, d in enumerate(self.domains):
            Ai_d, Aj_d, Av_d = d.assemble(f)
            bi_d, bv_d = d.apply_boundary_conditions(f)

            Ai.append(self.domains_dof[i_d][0]+Ai_d)
            Aj.append(self.domains_dof[i_d][0]+Aj_d)
            Av.append(Av_d)

            bi.append(self.domains_dof[i_d][0]+bi_d)
            bv.append(bv_d)

        # initialise the global matrix
        Ai = np.concatenate(Ai)
        Aj = np.concatenate(Aj)
        Av = np.concatenate(Av)
        self.A = coo_matrix((Av,(Ai, Aj)), shape=(self.max_nb_dof, self.max_nb_dof), dtype=np.complex)

        bi = np.concatenate(bi)
        bv = np.concatenate(bv)

        belem_remove_idxs = []
        actual_bi, actual_bv = [], []
        for i_belem, belem in enumerate(bv):
            if belem is None:
                belem_remove_idxs.append(bi[i_belem])
            else:
                actual_bi.append(bi[i_belem])
                actual_bv.append(belem)
        self.__apply_dof_remove(belem_remove_idxs)

        self.b = coo_matrix((actual_bv,(actual_bi, np.zeros((len(actual_bi))))), shape=(self.max_nb_dof,1), dtype=np.complex)

        self.x = spsolve(self.A.tocsc(), self.b)
        if output_matrices:
            return self.A, self.b, self.x
        else:
            return self.x
