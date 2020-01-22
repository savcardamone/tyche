#!/usr/bin/env python

"""atomic_orbital.py: Contraction of Gaussian primitives."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import sys
import numpy as np
import scipy.special as sp
import untangle

class AtomicOrbital():
    """Container for an atomic orbital, along with evaluation routines.
    """

    def __init__(self, centre=None, coeffs=None, zetas=None, ang_mom=(0,0,0)):
        """Class constructor.
        Initialise the atomic orbital's centre and primitive coefficients/exponents.
        We also initialise the angular momentum of the orbital, although default it to being
        spherically symmetric.
        """

        def normalise(self):
            """Normalise the atomic orbital so that it's self-overlap is unity. Taken verbatim from
            Fermann & Valeev document on integral evaluation, equation (2.25)
            """

            tot_ang_mom = ang_mom[0] + ang_mom[1] + ang_mom[2]
            prefactor_part = (sp.factorial2(2*ang_mom[0]-1)*sp.factorial2(2*ang_mom[1]-1)*sp.factorial2(2*ang_mom[2]-1)) / np.power(2, tot_ang_mom)
            prefactor = 1 / np.sqrt(np.power(np.pi, 1.5) * prefactor_part)

            summation = 0.0
            for iprim in range(len(self.coeffs)):
                for jprim in range(len(self.coeffs)):
                    summation += self.coeffs[iprim] * self.coeffs[jprim] / np.power(self.zetas[iprim] + self.zetas[jprim], tot_ang_mom + 1.5)

            N = prefactor / np.sqrt(summation)
            for iprim in range(len(self.coeffs)):
                self.coeffs[iprim] /= N

        if (type(centre) is not np.ndarray) or (type(coeffs) is not np.ndarray) or (type(zetas) is not np.ndarray):
            sys.exit("Input arguments to AtomicOrbital constructor must be numpy arrays.")
        if len(coeffs) != len(zetas):
            sys.exit("Nunber of coefficients must match number of exponents.")
        if len(centre) != 3 or len(ang_mom) != 3:
            sys.exit("Atomic orbital must live in 3D space.")

        self.centre = centre
        self.coeffs = coeffs
        self.zetas = zetas
        self.ang_mom = ang_mom

        # Normalise the primitive expansion
        normalise(self)

    def __str__(self):
        """Object string representation.
        """

        np.set_printoptions(formatter={'float': '{:7.4f}'.format})
        return "Atomic Orbital with {0} primitives.\n" \
            "   Centre      : {1}\n" \
            "   Coefficients: {2}\n" \
            "   Zeta        : {3}\n" \
            "   Ang. Mom.   : {4}".format(len(self.coeffs), self.centre, self.coeffs, self.zetas, self.ang_mom)

    @classmethod
    def from_atom(cls, atom, basis_set_file=None):
        
        xml = untangle.parse(basis_set_file)
        atomic_basis_sets = xml.Wavefunction.Basis

        aos = []
        for atomic_basis_set in atomic_basis_sets:
            if atomic_basis_set['atom'] == atom.atom_type:
                for contraction in atomic_basis_set.Contraction:
                    aos.append(cls(
                        np.transpose(atom.pos[0]),
                        np.fromstring(contraction.Zeta.cdata, dtype=float, sep=','),
                        np.fromstring(contraction.Coeff.cdata, dtype=float, sep=','),
                        np.fromstring(contraction.AngMom.cdata, dtype=int, sep=',')
                    ))
                if len(aos) != int(atomic_basis_set['num_funcs']):
                    sys.exit("Number of AOs parsed ({0}) does not equal the number designated ({1})".format(len(aos), int(atomic_basis_set['num_funcs'])))
                    
        if not aos:
            sys.exit("Could not find atomic basis for atom type {0}".format(atom.atom_type))
            
        return aos
            
    def evaluate(self, pos):
        """Evaluate the atomic orbital at a given position.
        """

        # Vector and square distance between atomic orbital centre and the position
        # at which we're evaluating it
        dr = pos - self.centre; dr_sq = norm(dr)

        ao_val = 0.0
        # Reduction over primitive values at given position
        for iprim in range(len(self.coeffs)):
            ao_val += self.coeffs[iprim] * np.exp(-self.zetas[iprim] * dr_sq)

        return np.power(dr[0], self.ang_mom[0]) * np.power(dr[1], self.ang_mom[1]) * np.power(dr[2], self.ang_mom[2]) * ao_val

    def gradient(self, pos):
        print("Gradient is currently unsupported.")

    def laplacian(self, pos):
        print("Laplacian is currently unsupported.")

def evaluate(aos, pos):
    """Evaluate each atomic orbital in the input list at the given position. Return a vector
    containing all atomic orbital values at that position.
    """

    ao_vals = zeros((len(aos), 1))
    for iao in range(len(aos)):
        ao_vals[iao] = aos[iao].evaluate(pos)

    return ao_vals

def gradient(aos, pos):
    """Evaluate the gradient of each atomic orbital in the input list at the given position.
    Return an array containing the gradient of all atomic orbital values at that position.
    """

    ao_grads = zeros((len(aos), 3))
    for iao in range(len(aos)):
        ao_grads[iao,:] = aos[iao].gradient(pos)

    return ao_grads

def laplacian(aos, pos):
    """Evaluate the laplacian of each atomic orbital in the input list at the given position.
    Return a vector containing all atomic orbital laplacians at that position.
    """

    ao_lapls = zeros((len(aos), 1))
    for iao in range(len(aos)):
        ao_lapls[iao] = aos[iao].laplacian(pos)

    return ao_lapls
