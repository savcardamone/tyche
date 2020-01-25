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

    # Class variables (static across all Slater wavefunctions)
    aos = None; mos = None; mo_coeffs = None
    
    @classmethod
    def initialise_parameters(cls, system, input_file):
        """Initialisation of class variables.
        The use case of this class is to construct a Slater object for each walker. If things
        like the LCAO and molecular orbital coefficients are instance variables, we end up 
        replicating a lot of constant data. So rather we factor the constants out the class 
        variables, and require the user call this method before trying to do object construction.
        """
        # Keep hold of the atomic orbitals that form the LCAO
        cls.aos = lcao.LCAO(system, input_file)

        # Pull out the level of the xml we want to be looping through
        xml = untangle.parse(input_file)
        molecular_orbitals = xml.Input.Wavefunction.MolecularOrbitals

        # Extract the number of molecular orbitals of each spin-type
        # Note that we use the number of MOs interchangeable with the number of electrons
        # of each spin type, else we end up with non-square Slater matrices, which doesn't really
        # make sense since the MOs are single-particle orbitals
        cls.num_mos = (int(molecular_orbitals['num_alpha']), int(molecular_orbitals['num_beta']))
        cls.mo_coeffs = (
            np.zeros((cls.num_mos[0], cls.aos.num_aos), dtype=float),
            np.zeros((cls.num_mos[1], cls.aos.num_aos), dtype=float)
        )

        imo_a = 0; imo_b = 0
        # Extract the coefficients for each molecular orbital. The user is allowed to designate
        # each MO as being of a particular spin, but we always verify at the end that we've
        # parsed as many as we expect
        for molecular_orbital in molecular_orbitals.Coefficients:
            coeffs = np.fromstring(molecular_orbital.cdata, dtype=float, sep=',')
            if molecular_orbital['spin'] == "Alpha" or molecular_orbital['spin'] == "Both":
                cls.mo_coeffs[0][imo_a,:] = coeffs; imo_a += 1
            if molecular_orbital['spin'] == "Beta"  or molecular_orbital['spin'] == "Both":
                cls.mo_coeffs[1][imo_b,:] = coeffs; imo_b += 1
                
        if imo_a != cls.num_mos[0] or imo_b != cls.num_mos[1]:
            sys.exit("Discrepancy in number of molecular orbitals in input file.")
        
    def __init__(self):
        """Class constructor.
        """
        # Verify that we've called the class variable initialisation routine before we start
        # trying to do object creation
        if None in (Slater.aos, Slater.num_mos, Slater.mo_coeffs):
            sys.exit("Slater class variables must be initialised before object construction.")

        # If there are no electrons of either spin state, then we force the respective Slater
        # matrix to be a scalar that's always unity. Else when we do the product of Slater
        # determinants we get gibberish out
        if Slater.num_mos[0] == 0:
            self.alpha = np.zeros((1, 1), dtype=float)
            self.inv_alpha = np.zeros((1, 1), dtype=float)
        else:
            self.alpha = np.zeros((Slater.num_mos[0], Slater.num_mos[0]), dtype=float)
            self.inv_alpha = np.zeros((Slater.num_mos[0], Slater.num_mos[0]), dtype=float)

        if Slater.num_mos[1] == 0:
            self.beta = np.zeros((1, 1), dtype=float)
            self.inv_beta = np.zeros((1, 1), dtype=float)
        else:
            self.beta = np.zeros((Slater.num_mos[1], Slater.num_mos[1]), dtype=float)
            self.inv_beta = np.zeros((Slater.num_mos[1], Slater.num_mos[1]), dtype=float)

        # Since we don't have a configuration yet, just have the determinants as unity
        self.alpha_det = 1.0; self.beta_det = 1.0        
        
    def __str__(self):
        """String representation of the Slater wavefunction.
        """
        np.set_printoptions(formatter={'float': '{:7.4f}'.format})
        walker_str = ""
        walker_str += "Slater Wavefunction\n"
        walker_str += " * Alpha MO Coefficients\n   {0}\n".format(self.mo_coeffs[0])
        walker_str += " * Beta MO Coefficients\n   {0}\n".format(self.mo_coeffs[1])
        walker_str += " * Alpha Slater Matrix\n   {0}\n".format(self.alpha)
        walker_str += " * Beta Slater Matrix\n   {0}\n".format(self.beta)
        walker_str += " * Alpha Inverse Slater Matrix\n   {0}\n".format(self.inv_alpha)
        walker_str += " * Beta Inverse Slater Matrix\n   {0}\n".format(self.inv_beta)
        walker_str += " * Alpha Determinant : {0:7.4f}\n".format(self.alpha_det)
        walker_str += " * Beta Determinant  : {0:7.4f}\n".format(self.beta_det)
        walker_str += " * Slater Determinant: {0:7.4f}\n".format(self.value())

    def value(self):
        """Return the product of spin-Slater determinants.
        """
        return self.alpha_det * self.beta_det
        
    def matrix(self, pos):
        """Evaluate the spin-Slater matrices for all electrons that are given as an argument.
        """
        # Just skip Slater matrix evaluation if there are no electrons of this spin-state
        if Slater.num_mos[0] > 0:
            # Evaluate the atomic orbitals for each electron of the spin-state. This gives
            # us an array of size (num_aos x num_elecs_spin)
            alpha_elec_aos = np.apply_along_axis(Slater.aos.evaluate, axis=1, arr=pos[0]).T
            # The Slater matrix can now be formed my matrix multiplication:
            # (num_mos_spin,num_aos) . (num_aos,num_elecs_spin) = (num_mos_spin,num_elecs_spin)
            self.alpha = np.dot(Slater.mo_coeffs[0], alpha_elec_aos)

        # As above, but for beta-spin electrons
        if Slater.num_mos[1] > 0:
            beta_elec_aos  = np.apply_along_axis(Slater.aos.evaluate, axis=1, arr=pos[1]).T
            self.beta = np.dot(Slater.mo_coeffs[1], beta_elec_aos)

        return (self.alpha, self.beta)
            
    def laplacian(self, pos):
        """Evaluate the spin-Slater laplacians for all electrons that are given as an argument.
        """
        # Just skip laplacian evaluation if there are no electrons of this spin-state
        if Slater.num_mos[0] > 0:
            # Evaluate the second derivatives of the atomic orbitals for each electron of the
            # spin-state. This gives us an array of size (num_aos x num_elecs)
            alpha_elec_ao_lapls = np.apply_along_axis(Slater.aos.laplacian, axis=1, arr=pos[0]).T
            # The laplacian of the Slater wavefunction can now be formed my matrix multiplication:
            # (num_mos_spin,num_aos) . (num_aos,num_elecs) = (num_mos_spin,num_elecs)
            alpha_laplacian = np.dot(Slater.mo_coeffs[0], alpha_elec_aos)

        # As above, but for beta-spin electrons
        if Slater.num_mos[1] > 0:
            beta_elec_ao_lapls  = np.apply_along_axis(Slater.aos.laplacian, axis=1, arr=pos[1]).T
            beta_laplacian = np.dot(Slater.mo_coeffs[1], beta_elec_aos)

        return (alpha_laplacian, beta_laplacian)

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
