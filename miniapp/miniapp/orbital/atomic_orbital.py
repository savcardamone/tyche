#!/usr/bin/env python

"""atomic_orbital.py: Contraction of Gaussian primitives."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import sys
import numpy as np
from scipy.special import factorial2 as dfac

class AtomicOrbital():
    """Container for an atomic orbital, along with evaluation routines.
    """

    def __init__(self, centre, coeffs, zetas, ang_mom):
        """Class constructor.
        Initialise the atomic orbital's centre, primitive coefficients/exponents and
        angular momentum. Normalise the AO if the user wants.
        """        
        if len(coeffs) == len(zetas):
            self.num_prims = len(coeffs)
        else:
            sys.exit("Length of AO coefficient and zetas arrays must be the same.")

        self.centre = centre
        self.coeffs = coeffs
        self.zetas = zetas
        self.ang_mom = ang_mom

    def __str__(self):
        """Object string representation.
        """
        np.set_printoptions(formatter={'float': '{:9.4f}'.format})
        ao_string = ""
        ao_string += "Atomic Orbital with {0} primitives.\n".format(self.num_prims)
        ao_string += "   Centre      : {0}\n".format(self.centre) 
        ao_string += "   Coefficients: {0}\n".format(self.coeffs) 
        ao_string += "   Zeta        : {0}\n".format(self.zetas) 
        ao_string += "   Ang. Mom.   : {0}\n".format(self.ang_mom) 
        return ao_string

    @staticmethod
    def overlap(orb_a, orb_b):
        """Compute the overlap between two orbital contractions.
        """
        acc = 0
        for aprim in range(orb_a.num_prims):
            for bprim in range(orb_b.num_prims):
                gamma = orb_a.zetas[aprim] + orb_b.zetas[bprim]
                zeta_prod = orb_a.zetas[aprim] * orb_b.zetas[bprim]
                coeff_prod = orb_a.coeffs[aprim] * orb_b.coeffs[bprim]
                _, sq_dist = orb_a.distance(orb_b.centre)
                norm = np.power(np.pi / gamma, 1.5)
                acc += norm * coeff_prod * np.exp(-zeta_prod * sq_dist / gamma)
        return acc

    def normalise(self):
        """Normalise the atomic orbital so that it's self-overlap is unity. Taken verbatim from
        Fermann & Valeev document on integral evaluation, equation (2.25)
        """
        tot_ang_mom = np.sum(self.ang_mom)
        prefactor_x = dfac(2*self.ang_mom[0] - 1)
        prefactor_y = dfac(2*self.ang_mom[1] - 1)
        prefactor_z = dfac(2*self.ang_mom[2] - 1)
        prefactor_part = prefactor_x * prefactor_y * prefactor_z / np.power(2, tot_ang_mom)
        prefactor = np.power(np.pi, 1.5) * prefactor_part

        summation = 0.0
        for iprim in range(len(self.coeffs)):
            for jprim in range(len(self.coeffs)):
                coeff_prod = self.coeffs[iprim] * self.coeffs[jprim]
                denom = np.power(self.zetas[iprim] + self.zetas[jprim], tot_ang_mom + 1.5)
                summation += coeff_prod / denom

        N = 1.0 / np.sqrt(prefactor * summation)
        for iprim in range(len(self.coeffs)):
            self.coeffs[iprim] *= N

    def distance(self, pos):
        """Return vector and square distance between atomic orbital centre and the position
        at which we're evaluating it.
        """
        dr = pos - self.centre
        for iaxis in range(3):
            if dr[iaxis] == 0: dr[iaxis] += np.finfo(float).eps
        dr_sq = np.inner(dr, dr)
        return dr, dr_sq

    def evaluate(self, pos):
        """Evaluate the atomic orbital at a given position.
        """
        # Get vector and square distance between atomic orbital centre and evaluation position
        dr, dr_sq = self.distance(pos)
        
        ao_val = 0.0
        # Reduction over primitive values at given position
        for iprim in range(len(self.coeffs)):
            ao_val += self.coeffs[iprim] * np.exp(-self.zetas[iprim] * dr_sq)

        prefactor_x = np.power(dr[0], self.ang_mom[0]) 
        prefactor_y = np.power(dr[1], self.ang_mom[1]) 
        prefactor_z = np.power(dr[2], self.ang_mom[2])
        
        return  prefactor_x * prefactor_y * prefactor_z * ao_val

    def gradient(self, pos):
        sys.exit("Gradient is currently unsupported.")

    def laplacian(self, pos):
        """Evaluate the laplacian of the atomic orbital at a given position.
        """
        # Get vector and square distance between atomic orbital centre and evaluation position
        dr, dr_sq = self.distance(pos)

        ### Primitive and derivative evaluation
        ao_val = 0.0; ddr_ao_val = 0.0; ddr2_ao_val = 0.0
        # Reduction over primitive values and derivatives at given position
        for iprim in range(len(self.coeffs)):
            primitive_val = self.coeffs[iprim] * np.exp(-self.zetas[iprim] * dr_sq)
            ao_val += primitive_val
            ddr_ao_val += self.zetas[iprim] * primitive_val
            ddr2_ao_val += self.zetas[iprim] * self.zetas[iprim] * primitive_val

        ### We need lots of powers of the displacement vector components to evaluate the Laplacian
        xl = dr[0] ** self.ang_mom[0]; xl_up = xl * dr[0] * dr[0]; xl_down = xl / (dr[0] * dr[0])
        ym = dr[1] ** self.ang_mom[1]; ym_up = ym * dr[1] * dr[1]; ym_down = ym / (dr[1] * dr[1])
        zn = dr[2] ** self.ang_mom[2]; zn_up = zn * dr[2] * dr[2]; zn_down = zn / (dr[2] * dr[2])
            
        ### First part of Laplacian evaluation -- Angular momentum reduced by 2
        ang_mom_down_x = self.ang_mom[0] * (self.ang_mom[0] - 1) * xl_down * ym * zn
        ang_mom_down_y = self.ang_mom[1] * (self.ang_mom[1] - 1) * xl * ym_down * zn
        ang_mom_down_z = self.ang_mom[2] * (self.ang_mom[2] - 1) * xl * ym * zn_down
        ang_mom_down = ang_mom_down_x * ang_mom_down_y * ang_mom_down_z * ao_val
 
        ### Second part of Laplacian evaluation -- Angular momentum increased by 2
        ang_mom_up_x = xl_up * ym * zn
        ang_mom_up_y = xl * ym_up * zn
        ang_mom_up_z = xl * ym * zn_up
        ang_mom_up = 4 * ang_mom_up_x * ang_mom_up_y * ang_mom_up_z * ddr2_ao_val

        ### Third part of Laplacian evaluation -- Angular momentum stays the same
        ang_mom_same = 4 * (np.sum(self.ang_mom) + 6) * xl * ym * zn * ddr_ao_val

        ### Final part of Laplacian evaluation -- throw everything together
        return ang_mom_down + ang_mom_up - ang_mom_same
