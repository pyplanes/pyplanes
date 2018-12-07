#! /usr/bin/env python
# -*- coding:utf8 -*-
#
# eqf.py
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
from numpy.lib.scimath import sqrt

from .fluid import Fluid
from .air import Air


class EqFluidJCA(Fluid):
    """ Represents an equivalent fluid medium with the Johnson-Champoux-Allard approach

    The Johnson-Champoux-Allard model is described in references [1]_, [2]_ and [3]_.

    Attributes
    ----------
    phi : float
        Porosity
    sigma : float
        Flow resistivity
    alpha : float
        Tortuosity
    Lambda_prime : float
        Thermal characteristic length
    Lambda : float
        Viscous characteristic length
    rho_1 : float
        Mass of solid per unit volume of aggregate
    nu : float
        poisson ratio
    E : float
        Young's modulus
    eta : float
        viscosity

    Notes
    -----

    Populates the instance with the following attributes (and a few others left undocumented)

    `rho_eq_til`
        Frequency dependent equivalent density
    `alpha_til`
        Frequency dependent Johnson's tortuosity
    `alpha_prime_til`
        Frequency dependent Champoux-Allard's tortuosity
    `K_eq_til`
        Frequency dependent equivalent compressibility
    `c_eq_til`
        Frequency dependent equivalent speed of sound

    References
    ----------

    .. [1] D. L. Johnson *et al.*, "Theory of dynamic permeability and tortuosity in fluid-saturated porous media", Journal of fluid mechanics 176(1), 1987.
    .. [2] Y. Champoux & J.-F. Allard, "Dynamic tortuosity and bulk modulus in air-saturated porous media", Journale of Applied Physics, 1975, DOI: 10.1063/1.349482.
    .. [3] J.-F. Allard & N. Atalla, "Propagation of Sound in porous media: modelling sound absorbing materials", Wiley, 2009, ISBN: 978-0-470-74661-5.
    """

    MEDIUM_TYPE = 'eqf'
    MODEL = 'fluid'
    EXPECTED_PARAMS = [
        ('phi', float),  # Porosity
        ('sigma', float),  # Flow resistivity
        ('alpha', float),  # Tortuosity
        ('Lambda_prime', float),  # Thermal characteristic length
        ('Lambda', float),  # Viscous characteristic length
        ('rho_1', float),  # Mass of solid per unit volume of aggregate
        ('nu', float),  # poisson ratio
        ('E', float),  # Young's modulus
        ('eta', float)  # viscosity
    ]

    def __init__(self):
        super().__init__()

        self.phi = None
        self.sigma = None
        self.alpha = None
        self.Lambda_prime = None
        self.Lambda = None
        self.rho_1 = None
        self.nu = None
        self.E = None
        self.N = None
        self.eta = None

    def _compute_missing(self):
        """ Computes the required constant parameters missing from the definition

        For a JCA equivalent fluid, the shear modulus `N` is computed."""
        self.N = self.E/(2*(1+self.nu))

    def update_frequency(self, f):
        """ Computes the JCA parameters (see Notes on the class).

        Parameters
        ----------

        f :
            Frequency of interest
        """

        omega = 2*np.pi*f

        #  Johnson et al model for rho_eq_til
        self.omega_0 = self.sigma*self.phi/(Air.rho*self.alpha)
        self.omega_infty = (self.sigma*self.phi*self.Lambda)**2/(4*Air.mu*Air.rho*self.alpha**2)
        self.F_JKD = sqrt(1+1j*omega/self.omega_infty)
        self.rho_eq_til = (Air.rho*self.alpha/self.phi)*(1+(self.omega_0/(1j*omega))*self.F_JKD)
        self.alpha_til = self.phi*self.rho_eq_til/Air.rho

        #  Champoux-Allard model for K_eq_til
        self.omega_prime_infty = (16*Air.nu_prime)/(self.Lambda_prime**2)
        self.F_prime_CA = sqrt(1+1j*omega/self.omega_prime_infty)
        self.alpha_prime_til = 1+self.omega_prime_infty*self.F_prime_CA/(2*1j*omega)
        self.K_eq_til = (Air.gamma*Air.P/self.phi)/(Air.gamma-(Air.gamma-1)/self.alpha_prime_til)

        self.c_eq_til = sqrt(self.K_eq_til/self.rho_eq_til)

        self.c = self.c_eq_til
        self.rho = self.rho_eq_til
