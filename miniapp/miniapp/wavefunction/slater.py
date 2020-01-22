
import sys
import numpy as np
import itertools

from miniapp.orbital import atomic_orbital as ao

class Slater():
    """Slater wavefunction class. Keeps alpha and beta spin components of the
    wavefunction separate.
    """
    
    def __init__(self, pos, system, mo_coeffs=None):
        """Class constructor.
        """

        def build_aos(system):
            """Construct the LCAO from the atomic information of the system.
            """
            
            aos = []
            for iatom in range(system.num_atoms):
                aos.append(ao.AtomicOrbital.from_atom(system.atoms[iatom]))
            return list(itertools.chain.from_iterable(aos))
            
        # Keep hold of the atomic orbitals that form the LCAO
        self.aos = build_aos(system)

        num_mos = (0, 0)
        if mo_coeffs == None:
            print("Reading MO coefficients from file currently unsupported.")
            # Fetch from file
        else:
            num_mos = (mo_coeffs[0].shape[0], mo_coeffs[1].shape[0]) 
            self.mo_coeffs = mo_coeffs

        # Number of alpha and beta electrons
        num_alpha, num_beta = pos[0].shape[0], pos[1].shape[0]
        if num_alpha > num_mos[0] or num_beta > num_mos[1]:
            sys.exit("More electrons than there are molecular orbitals.")

        # Make sure we don't allocate an empty array -- initialise spin-Slater
        # matrices as unity if there are no electrons in the spin state, else
        # explicitly allocate
        if num_alpha == 0:
            self.alpha = np.ones((1, 1))
        else:
            self.alpha = np.zeros((num_alpha, num_mos[0]))
            
        if num_beta == 0:
            self.beta = np.ones((1, 1))
        else:
            self.beta = np.zeros((num_beta,  num_mos[1]))
            
        # Evaluate the spin-Slater determinants
        self.evaluate(pos)

    def __str__(self):
        """String representation of the Slater wavefunction.
        """

        np.set_printoptions(formatter={'float': '{:7.4f}'.format})
        return "Slater Wavefunction\n" \
            " * Alpha MO Coefficients\n   {0}\n * Beta MO Coefficients\n   {1}\n" \
            " * Alpha Electrons\n   {2}\n * Beta Electrons\n   {3}\n" \
            " * Alpha Slater Matrix\n   {4}\n * Beta Slater Matrix\n   {5}\n" \
            " * Alpha Determinant : {6:7.4f}\n * Beta Determinant  : {7:7.4f}\n * Slater Determinant: {8:7.4f}"\
            .format(self.mo_coeffs[0], self.mo_coeffs[1], self.pos[0], self.pos[1], \
                    self.alpha, self.beta, self.alpha_det, self.beta_det, self.value())

    def value(self):
        """Return the product of spin-Slater determinants.
        """
        return self.alpha_det * self.beta_det
        
    def evaluate(self, pos):
        """Evaluate the spin-Slater matrices and determinants for all electrons in
        the system. Note that this is brute force determinant calculating, and not
        Sherman-Morrison, so cubic scaling in number of electrons.
        """
        
        # Keep hold of the electron configuration giving rise to the wavefunction value
        self.pos = pos

        # Loop over alpha electrons, compute the LCAO at the electron's position
        # and evaluate the molecular orbitals
        for ialpha in range(pos[0].shape[0]):
            ao_vals = ao.evaluate(self.aos, self.pos[0][ialpha,:])
            self.alpha[ialpha,:] = np.dot(self.mo_coeffs[0], ao_vals)
        # Form the alpha-spin Slater determinant
        self.alpha_det = np.det(self.alpha)
            
        # Loop over beta electrons, compute the LCAO at the electron's position
        # and evaluate the molecular orbitals
        for ibeta in range(pos[1].shape[0]):
            ao_vals = ao.evaluate(self.aos, self.pos[1][ibeta,:])
            self.beta[ibeta,:] = np.dot(self.mo_coeffs[1], ao_vals)
        # Form the beta-spin Slater determinant
        self.beta_det = np.det(self.beta)

    def update(self, pos, iel):
        """Sherman-Morrison update of the appropriate spin-Slater matrix and
        corresponding determinant.
        """
        print("Sherman-Morrison updating unsupported.")
