#!/usr/bin/env python

"""slater.py: Slater matrix and determinant manipulation."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import sys
import numpy as np
import untangle

from miniapp.orbital import lcao

class Slater():
    """Slater wavefunction class.

    Since the wavefunction is factorised into product of spin-determinants, many of the routines
    return a tuple of values -- the first element in the tuple corresponds to the alpha-spin
    part while the second element, the beta-spin part.
    """

    def __init__(self, system, input_file):
        """Class constructor. Construct the LCAO and initialise the MO coefficients.
        """
        # Keep hold of the atomic orbitals that form the LCAO
        self.aos = lcao.LCAO(system, input_file)

        # Pull out the level of the xml we want to be looping through
        xml = untangle.parse(input_file)
        molecular_orbitals = xml.Input.Wavefunction.MolecularOrbitals

        # Extract the number of molecular orbitals of each spin-type
        # Note that we use the number of MOs interchangeably with the number of electrons
        # of each spin type, else we end up with non-square Slater matrices, which doesn't really
        # make sense (I don't think?), else an MO isn't a single-particle orbitals (which
        # is contrary to their definition)
        self.num_mos = (int(molecular_orbitals['num_alpha']), int(molecular_orbitals['num_beta']))
        self.mo_coeffs = (
            np.zeros((self.aos.num_aos, self.num_mos[0]), dtype=float),
            np.zeros((self.aos.num_aos, self.num_mos[1]), dtype=float)
        )

        # Some counters
        imo_a = 0
        imo_b = 0

        # Extract the coefficients for each molecular orbital. The user is allowed to designate
        # each MO as being of a particular spin, but we always verify at the end that we've
        # parsed as many as we expect
        for molecular_orbital in molecular_orbitals.Coefficients:
            coeffs = np.fromstring(molecular_orbital.cdata, dtype=float, sep=',')
            if molecular_orbital['spin'] == "Alpha" or molecular_orbital['spin'] == "Both":
                self.mo_coeffs[0][:, imo_a] = coeffs
                imo_a += 1
            if molecular_orbital['spin'] == "Beta"  or molecular_orbital['spin'] == "Both":
                self.mo_coeffs[1][:, imo_b] = coeffs
                imo_b += 1
        if imo_a != self.num_mos[0] or imo_b != self.num_mos[1]:
            sys.exit("Discrepancy in number of molecular orbitals in input file.")
        
        # The MO coefficients should transform the AO overlap matrix into the
        # identity (since they're mutually orthonormal)
        for imo in range(self.num_mos[0]):
            mo_coeffs = self.mo_coeffs[0][:, imo]
            norm = mo_coeffs.T @ self.aos.overlap_matrix() @ mo_coeffs
            self.mo_coeffs[0][:, imo] /= np.sqrt(norm)
        for imo in range(self.num_mos[1]):
            mo_coeffs = self.mo_coeffs[1][:, imo]
            norm = mo_coeffs.T @ self.aos.overlap_matrix() @ mo_coeffs
            self.mo_coeffs[1][:, imo] /= np.sqrt(norm)

    def __str__(self):
        """String representation of the Slater wavefunction.
        """
        np.set_printoptions(formatter={'float': '{:7.4f}'.format})
        slater_str = "Slater Wavefunction\n"
        slater_str += self.aos.__str__()
        slater_str += "Alpha MO Coefficients\n   {0}\n".format(self.mo_coeffs[0])
        slater_str += "Beta MO Coefficients\n   {0}\n".format(self.mo_coeffs[1])
        return slater_str

    def density_matrix(self):
        r"""Compute the density matrices for the alpha- and beta-spin wavefunctions.
        .. math::
            P_{\mu\nu}^\omega = \sum_a^{N_\omega}C_{\mu a}^\omega C_{\nu a}^\omega\,,

        where :math:`P_{\mu\nu}^\omega` is an element of the spin-density matrix,
        :math:`N_\omega` the number of molecular orbitals in the spin state and
        :math:`\mathbf{C}` the matrix of molecular orbital coefficients.
        """
        return (
            self.mo_coeffs[0] @ self.mo_coeffs[0].T,
            self.mo_coeffs[1] @ self.mo_coeffs[1].T
        )

    def density(self, pos):
        r"""Evaluate the electronic density at a point in :math:`\mathbb{R}^3`.
        We do so through the density matrix:
        .. math::
            \rho^\omega(\mathbf{r}) =
            \sum_{\mu,\nu} P_{\mu\nu}^\omega\phi_\mu(\mathbf{r})\phi_\nu(\mathbf{r})\,,

        where :math:`\rho^\omega(\mathbf{r})` is the density from the spin-state at
        the specified point, and :math:`\phi(\mathbf{r})` an atomic orbital evaluated
        at that point. The density is returned as the sum of contributions from the
        spin-states.
        """
        (density_matrix_alpha, density_matrix_beta) = self.density_matrix()
        ao_vals = self.aos.evaluate(pos)
        density_alpha = ao_vals.T @ density_matrix_alpha @ ao_vals
        density_beta = ao_vals.T @ density_matrix_beta @ ao_vals
        return density_alpha + density_beta

    def density_field(self, xs, ys, zs):
        """Evaluate the density scalar field on a grid.
        """
        field = np.zeros(xs.shape, dtype=float)
        for x in range(xs.shape[0]):
            for y in range(ys.shape[0]):
                for z in range(zs.shape[0]):
                    field[x,y,z] = self.density(np.array([xs[x,y,z], ys[x,y,z], zs[x,y,z]]))
        return field
    
    def matrix(self, pos):
        """Evaluate the spin-Slater matrices for all electrons that are given as an argument.
        """
        # If there are no electrons in a spin state, the Slater matrix is just unity, else
        # we explicitly evaluate it
        if self.num_mos[0] > 0:
            # Evaluate the atomic orbitals for each electron of the spin-state. This gives
            # us an array of size (num_aos x num_elecs_spin)
            alpha_elec_aos = np.apply_along_axis(self.aos.evaluate, axis=1, arr=pos[0]).T
            # The Slater matrix can now be formed my matrix multiplication:
            # (num_mos_spin,num_aos) . (num_aos,num_elecs_spin) = (num_mos_spin,num_elecs_spin)
            alpha = np.dot(self.mo_coeffs[0], alpha_elec_aos)
        else:
            alpha = np.ones((1,1), dtype=float)
            
        # As above, but for beta-spin electrons
        if self.num_mos[1] > 0:
            beta_elec_aos  = np.apply_along_axis(self.aos.evaluate, axis=1, arr=pos[1]).T
            beta = np.dot(self.mo_coeffs[1], beta_elec_aos)
        else:
            beta = np.ones((1,1), dtype=float)

        return (alpha, beta)
            
    def laplacian(self, pos):
        """Evaluate the spin-Slater laplacians for all electrons that are given as an argument.
        """
        # If there are no electrons in a spin state, the laplacian of the wavefunction is zero
        # to force any kinetic energy contributions to be zero
        if self.num_mos[0] > 0:
            # Evaluate the second derivatives of the atomic orbitals for each electron of the
            # spin-state. This gives us an array of size (num_aos x num_elecs)
            alpha_elec_ao_lapls = np.apply_along_axis(self.aos.laplacian, axis=1, arr=pos[0]).T
            # The laplacian of the Slater wavefunction can now be formed my matrix multiplication:
            # (num_mos_spin,num_aos) . (num_aos,num_elecs) = (num_mos_spin,num_elecs)
            alpha = np.dot(self.mo_coeffs[0], alpha_elec_ao_lapls)
        else:
            alpha = np.zeros((1,1), dtype=float)

        # As above, but for beta-spin electrons
        if self.num_mos[1] > 0:
            beta_elec_ao_lapls  = np.apply_along_axis(self.aos.laplacian, axis=1, arr=pos[1]).T
            beta = np.dot(self.mo_coeffs[1], beta_elec_ao_lapls)
        else:
            beta = np.zeros((1,1), dtype=float)

        return (alpha, beta)