#! /usr/bin/env python
# -*- coding:utf8 -*-
#
# pem.py
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

from .eqf import EqFluidJCA
from .air import Air


class PEM(EqFluidJCA):
    """Implementation of a Biot-Johnson-Champoux-Allard model representing a poroelastic
    material.

    The Johnson-Champoux-Allard model is described in references [1]_, [2]_ and [3]_. The
    Biot model is detailled in references [4]_, [5]_ and [6]_. The book of Allard and
    Atalla (Ref. [7]_) is a useful resource.

    Note that the implementation is based on the work of Dazel et al. (Ref. [8]_) and the
    inner variables are named to match these of this contribution.

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
    delta_1, delta_2, delta_3:
        Wavenumbers for the three waves (3 being the elastic shear wave)
    mu_1, mu_2, mu_3
        Amplitude ratii (see Ref. [8]_)
    delta_eq
        Wavenumber linked to the JCA model
    N
        Young modulus


    References
    ----------

    .. [1] D. L. Johnson *et al.*, "Theory of dynamic permeability and tortuosity in fluid-saturated porous media", Journal of fluid mechanics 176(1), 1987.
    .. [2] Y. Champoux & J.-F. Allard, "Dynamic tortuosity and bulk modulus in air-saturated porous media", Journale of Applied Physics, 1975, DOI: 10.1063/1.349482.
    .. [3] J.-F. Allard & N. Atalla, "Propagation of Sound in porous media: modelling sound absorbing materials", Wiley, 2009, ISBN: 978-0-470-74661-5.
    .. [4] M.A. Biot, "Theory of Propagation of Elastic Waves in a Fluid-Saturated Porous Solid", J. Acoustical Soc. of America 28(2), 1956.
    .. [5] M.A. Biot & D. Willis, "The Elastic Coefficients of the Theory of Consolidation", Journal of Applied Mechanics, 1957.
    .. [6] M.A. Biot, "Mechanics of Deformation and Acoustic Propagation in Porous Media", Journal of Applied Physics 33(4), 1962.
    .. [7] J-F Allard & N. Atalla, "Propagation of sound in porous media: modelling sound absorbing materials", Wiley, 2009 (ISBN: 9780470746615)
    .. [8] O. Dazel *et al.*, "An alternative Biot’s displacement formulation for porous materials", J. Acoustical Soc. of America 121(6), 2007.
    """

    MEDIUM_TYPE = 'pem'
    MODEL = 'pem'
    EXPECTED_PARAMS = [
        ('phi', float),  # Porosity
        ('sigma', float),  # Flow resistivity
        ('alpha', float),  # Tortuosity
        ('Lambda_prime', float),  # Thermal characteristic length
        ('Lambda', float),  # Viscous characteristic length
        ('rho_1', float),  # Mass of solid per unit volume of aggregate
        ('nu', float),  # poisson ratio
        ('E', float),  # Young's modulus
        ('eta', float),  # loss factor
        ('loss_type', str),  # type of losses for the frame (anelastic or structural)
    ]

    # def from_dict(self, *a, **kw):
    #     super().from_dict(*a, **kw)
    #     self.__compute_missing()

    def _compute_missing(self):
        super()._compute_missing()

        self.lambda_ = self.E*self.nu/((1+self.nu)*(1-2*self.nu))

        self.rho_12 = -self.phi*Air.rho*(self.alpha-1)
        self.rho_11 = self.rho_1-self.rho_12
        self.rho_2 = self.phi*Air.rho
        self.rho_22 = self.rho_2-self.rho_12

    def update_frequency(self, f):
        super().update_frequency(f)
        omega = 2*np.pi*f

        self.rho_22_til = self.phi**2*self.rho_eq_til
        self.rho_12_til = self.rho_2-self.rho_22_til
        self.rho_11_til = self.rho_1-self.rho_12_til
        self.rho_til = self.rho_11_til-((self.rho_12_til**2)/self.rho_22_til)

        self.gamma_til = self.phi*(self.rho_12_til/self.rho_22_til-(1-self.phi)/self.phi)
        self.rho_s_til = self.rho_til+self.gamma_til**2*self.rho_eq_til

        if self.loss_type == 'anelastic':
            raise NotImplementedError("Anelastic losses not implemented yet")
            # *HAS* to be in update_frequency() if anelastic is considered
            # self.structural_loss=1+(b_hat*(1j*omega/beta_hat)**alpha_hat)/(1+(1j*omega/beta_hat)**alpha_hat)
        elif self.loss_type == 'structural':
            self.structural_loss = 1+1j*self.eta
        elif self.loss_type == 'none':
            self.structural_loss = 1
        else:
            raise ValueError('Unknown type of losses')

        self.N = self.E/(2*(1+self.nu))*self.structural_loss
        self.A_hat = (self.E*self.nu)/((1+self.nu)*(1-2*self.nu))*self.structural_loss
        self.P_hat = self.A_hat+2*self.N

        # Biot 1956 elastic coefficients
        self.R_til = self.K_eq_til*self.phi**2
        self.Q_til = ((1-self.phi)/self.phi)*self.R_til
        self.P_til = self.P_hat+self.Q_til**2/self.R_til

        delta_eq = omega*sqrt(self.rho_eq_til/self.K_eq_til)
        delta_s_1 = omega*sqrt(self.rho_til/self.P_hat)
        delta_s_2 = omega*sqrt(self.rho_s_til/self.P_hat)

        Psi = ((delta_s_2**2+delta_eq**2)**2-4*delta_eq**2*delta_s_1**2)
        sdelta_total = sqrt(Psi)

        delta_1 = sqrt(0.5*(delta_s_2**2+delta_eq**2+sdelta_total))
        delta_2 = sqrt(0.5*(delta_s_2**2+delta_eq**2-sdelta_total))
        delta_3 = omega*sqrt(self.rho_til/self.N)

        mu_1 = self.gamma_til*delta_eq**2/(delta_1**2-delta_eq**2)
        mu_2 = self.gamma_til*delta_eq**2/(delta_2**2-delta_eq**2)
        mu_3 = -self.gamma_til

        self.delta_s_1 = delta_s_1
        self.delta_s_2 = delta_s_2
        self.delta_1 = delta_1
        self.delta_2 = delta_2
        self.delta_3 = delta_3
        self.delta_eq = delta_eq
        self.sdelta_total = sdelta_total
        self.mu_1 = mu_1
        self.mu_2 = mu_2
        self.mu_3 = mu_3
