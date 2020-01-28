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
    return a tuple of values -- the zeroth element in the tuple corresponds to the alpha-spin
    part while the first element the beta-spin part.
    """

    def __init__(self, system, input_file):
        """We have no instance variables in this class, so just verify the class variables have
        been initialised properly before object creation.
        """
        # Keep hold of the atomic orbitals that form the LCAO
        self.aos = lcao.LCAO(system, input_file)

        # Pull out the level of the xml we want to be looping through
        xml = untangle.parse(input_file)
        molecular_orbitals = xml.Input.Wavefunction.MolecularOrbitals

        # Extract the number of molecular orbitals of each spin-type
        # Note that we use the number of MOs interchangeable with the number of electrons
        # of each spin type, else we end up with non-square Slater matrices, which doesn't really
        # make sense since the MOs are single-particle orbitals
        self.num_mos = (int(molecular_orbitals['num_alpha']), int(molecular_orbitals['num_beta']))
        self.mo_coeffs = (
            np.zeros((self.num_mos[0], self.aos.num_aos), dtype=float),
            np.zeros((self.num_mos[1], self.aos.num_aos), dtype=float)
        )

        imo_a = 0; imo_b = 0
        # Extract the coefficients for each molecular orbital. The user is allowed to designate
        # each MO as being of a particular spin, but we always verify at the end that we've
        # parsed as many as we expect
        for molecular_orbital in molecular_orbitals.Coefficients:
            coeffs = np.fromstring(molecular_orbital.cdata, dtype=float, sep=',')
            if molecular_orbital['spin'] == "Alpha" or molecular_orbital['spin'] == "Both":
                self.mo_coeffs[0][imo_a,:] = coeffs; imo_a += 1
            if molecular_orbital['spin'] == "Beta"  or molecular_orbital['spin'] == "Both":
                self.mo_coeffs[1][imo_b,:] = coeffs; imo_b += 1
                
        if imo_a != self.num_mos[0] or imo_b != self.num_mos[1]:
            sys.exit("Discrepancy in number of molecular orbitals in input file.")
            
    def __str__(self):
        """String representation of the Slater wavefunction.
        """
        np.set_printoptions(formatter={'float': '{:7.4f}'.format})
        slater_str = "Slater Wavefunction\n"
        slater_str += self.aos.__str__()
        slater_str += "Alpha MO Coefficients\n   {0}\n".format(self.mo_coeffs[0])
        slater_str += "Beta MO Coefficients\n   {0}\n".format(self.mo_coeffs[1])
        return slater_str
            
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

    ### TODO:
    ### These methods should go in the Walker class
    
    def inverse_explicit(self):
        """Explicitly invert the Slater matrices via brute force cubic-scaling method.
        Should only really be called at object creation.
        """
        if Slater.num_mos[0] > 0:
            self.inv_alpha = np.linalg.inv(self.alpha)
        
        if Slater.num_mos[1] > 0:
            self.inv_beta = np.linalg.inv(self.beta)

        return (self.inv_alpha, self.inv_beta)
            
    def update(self, new_pos, iel, spin):
        """Sherman-Morrison update of the appropriate spin-Slater matrix and
        corresponding determinant.
        """
        if ispin == "Alpha":
            new_mos = np.dot(Slater.mo_coeffs[0], np.transpose(Slater.aos.evaluate(new_pos)))
            delta_mos = new_mos - self.alpha[:,iel] 
        elif ispin == "Beta":
            new_mos = np.dot(Slater.mo_coeffs[1], np.transpose(Slater.aos.evaluate(new_pos)))
            delta_mos = new_mos - self.beta[:,iel]
        else:
            sys.exit("Don't understand the spin-state {0}".format(spin))
            
        print("Sherman-Morrison updating unsupported.")
