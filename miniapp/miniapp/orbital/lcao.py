#!/usr/bin/env python

"""lcao.py: Linear combination of atomic orbitals."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import sys
import itertools
import numpy as np
import untangle

from miniapp.orbital import atomic_orbital as ao

class LCAO():
    """Wrapper around list of atomic orbitals which constitute the LCAO.
    """

    def __init__(self, system, input_file):
        """Class constructor.
        Build up the LCAO by iterating over atoms in the system and building up the
        appropriate AOs on each centre based on the input atomic basis set.
        """
        # Pull out the level of the xml we want to be looping through
        xml = untangle.parse(input_file)
        atomic_basis_sets = xml.Input.Wavefunction.Basis

        self.aos = []
        for atom in system.atoms:
            atomic_aos = []
            # Loop over atomic basis sets specified in the file until we find the appropriate one
            for atomic_basis_set in atomic_basis_sets:
                if atomic_basis_set['atom'] == atom.atom_type:

                    # Loop over AOs in this contraction and append each to the LCAO
                    # We explicitly construct the multiple atomic orbitals arising from l > 0
                    for contraction in atomic_basis_set.Contraction:

                        coeffs = np.fromstring(contraction.Coeff.cdata, dtype=float, sep=',')
                        zetas = np.fromstring(contraction.Zeta.cdata, dtype=float, sep=',')

                        if contraction['ang_mom'] == "s":
                            ang_moms = [np.array([0, 0, 0], dtype=int)]
                        elif contraction['ang_mom'] == "p":
                            ang_moms = [
                                np.array([1, 0, 0], dtype=int), np.array([0, 1, 0], dtype=int),
                                np.array([0, 0, 1], dtype=int)
                            ]
                        elif contraction['ang_mom'] == "d":
                            ang_moms = [
                                np.array([2, 0, 0], dtype=int), np.array([0, 2, 0], dtype=int),
                                np.array([0, 0, 2], dtype=int), np.array([1, 1, 0], dtype=int),
                                np.array([1, 0, 1], dtype=int), np.array([0, 1, 1], dtype=int)
                            ]
                        else:
                            sys.exit("Unrecognised contraction ang mom: {0}".format(
                                contraction['ang_mom']))

                        # Loop over all (-l <= m <= l) and construct the corresponding
                        # atomic orbital
                        for ang_mom in ang_moms:
                            atomic_aos.append(
                                ao.AtomicOrbital(atom.pos, coeffs, zetas, ang_mom)
                            )

            if not atomic_aos:
                sys.exit("Could not find atomic basis for atom type {0}".format(atom.atom_type))
            else:
                self.aos.append(atomic_aos)

        self.aos = list(itertools.chain.from_iterable(self.aos))
        self.num_aos = len(self.aos)

    def __str__(self):
        """Object string representation.
        """
        lcao_str = "Linear Combination of Atomic Orbitals\n"
        lcao_str += "Number of AOs: {0}\n".format(self.num_aos)
        for iao in self.aos:
            lcao_str += iao.__str__()

        return lcao_str

    def overlap_matrix(self):
        """Compute the overlap matrix for all AOs in the LCAO.
        """
        matrix = np.zeros((self.num_aos, self.num_aos), dtype=float)
        for iao in range(self.num_aos):
            for jao in range(self.num_aos):
                matrix[iao, jao] = ao.AtomicOrbital.overlap(self.aos[iao], self.aos[jao])
        return matrix

    def evaluate(self, pos):
        """Evaluate all atomic orbitals in the LCAO at a given position.
        """
        ao_vals = np.zeros((self.num_aos,), dtype=float)
        for iao in range(self.num_aos):
            ao_vals[iao] = self.aos[iao].evaluate(pos)

        return ao_vals

    def gradient(self, pos):
        """Evaluate the gradient of all atomic orbitals in the LCAO at a given position.
        """
        ao_grads = np.zeros((self.num_aos, 3), dtype=float)
        for iao in range(self.num_aos):
            ao_grads[iao, :] = self.aos[iao].gradient(pos)

        return ao_grads

    def laplacian(self, pos):
        """Evaluate the laplacian of all atomic orbitals in the LCAO at a given position.
        """
        ao_lapls = np.zeros((self.num_aos,), dtype=float)
        for iao in range(self.num_aos):
            ao_lapls[iao] = self.aos[iao].laplacian(pos)

        return ao_lapls
