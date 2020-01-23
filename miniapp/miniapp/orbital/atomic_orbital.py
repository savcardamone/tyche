#!/usr/bin/env python

"""atomic_orbital.py: Contraction of Gaussian primitives."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import sys
import numpy as np
from scipy.special import factorial2 as dfac
import untangle

class AtomicOrbital():
    """Container for an atomic orbital, along with evaluation routines.
    """

    def __init__(self, centre, coeffs, zetas, ang_mom, normalise_ao=True):
        """Class constructor.
        Initialise the atomic orbital's centre, primitive coefficients/exponents and
        angular momentum. Normalise the AO if the user wants.
        """        
        self.centre = centre
        self.coeffs = coeffs
        self.zetas = zetas
        self.ang_mom = ang_mom

        # Normalise the primitive expansion
        if normalise_ao == True:
            self.normalise()

    def __str__(self):
        """Object string representation.
        """
        np.set_printoptions(formatter={'float': '{:7.4f}'.format})
        return "Atomic Orbital with {0} primitives.\n" \
            "   Centre      : {1}\n" \
            "   Coefficients: {2}\n" \
            "   Zeta        : {3}\n" \
            "   Ang. Mom.   : {4}".format(len(self.coeffs), self.centre, self.coeffs, self.zetas, self.ang_mom)

    def normalise(self):
        """Normalise the atomic orbital so that it's self-overlap is unity. Taken verbatim from
        Fermann & Valeev document on integral evaluation, equation (2.25)
        """
        tot_ang_mom = np.sum(self.ang_mom)
        prefactor_part = (dfac(2*self.ang_mom[0]-1)*dfac(2*self.ang_mom[1]-1)*dfac(2*self.ang_mom[2]-1)) / np.power(2, tot_ang_mom)
        prefactor = 1 / np.sqrt(np.power(np.pi, 1.5) * prefactor_part)
        
        summation = 0.0
        for iprim in range(len(self.coeffs)):
            for jprim in range(len(self.coeffs)):
                summation += self.coeffs[iprim] * self.coeffs[jprim] / np.power(self.zetas[iprim] + self.zetas[jprim], tot_ang_mom + 1.5)

        N = prefactor / np.sqrt(summation)
        for iprim in range(len(self.coeffs)):
            self.coeffs[iprim] /= N
    
    def evaluate(self, pos):
        """Evaluate the atomic orbital at a given position.
        """
        # Vector and square distance between atomic orbital centre and the position
        # at which we're evaluating it
        dr = pos - self.centre; dr_sq = np.linalg.norm(dr)

        ao_val = 0.0
        # Reduction over primitive values at given position
        for iprim in range(len(self.coeffs)):
            ao_val += self.coeffs[iprim] * np.exp(-self.zetas[iprim] * dr_sq)

        return np.power(dr[0], self.ang_mom[0]) * np.power(dr[1], self.ang_mom[1]) * np.power(dr[2], self.ang_mom[2]) * ao_val

    def gradient(self, pos):
        sys.exit("Gradient is currently unsupported.")

    def laplacian(self, pos):
        sys.exit("Laplacian is currently unsupported.")
