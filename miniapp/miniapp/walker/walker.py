#!/usr/bin/env python

"""walker.py: QMC walker (or Monte Carlo sample, rather)."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import sys
import numpy as np
import itertools

class Walker():
    """Walker class. A realisation of a Markov chain which samples the electronic wavefunction.
    """
    
    # Class variables (static across all walkers)
    system = None; wavefunction = None
    
    @classmethod
    def initialise_parameters(cls, system, wavefunction):
        """Initialisation of class variables.
        We're going to need to system to evaluate nuclear contributions to the local energy,
        but this is constant data, so no need to replicate for each walker.
        """
        # Keep hold of the atomic orbitals that form the LCAO
        cls.system = system
        cls.wavefunction = wavefunction
        
    def __init__(self):
        """Class constructor.
        """
        # Verify that we've called the class variable initialisation routine before we start
        # trying to do object creation
        if None in (Walker.system, Walker.wavefunction):
            sys.exit("Walker class variables must be initialised before object construction.")
            
        def place_electron():
            """Place an electron within the van der Waals radius of a random atom,
            weighted by its nuclear charge.
            """
            random_charge = np.random.randint(low=0, high=Walker.system.tot_z)
            cumsum_charge = list(itertools.accumulate([atom.Z for atom in Walker.system.atoms]))
            atom_idx = [x > random_charge for x in cumsum_charge].index(True)
            # Return a Gaussian variate centred on the chosen atom
            return np.random.normal(
                Walker.system.atoms[atom_idx].pos, Walker.system.atoms[atom_idx].Rvdw, size=3
            )
        
        # Initialise the array to store the alpha and beta electrons
        self.pos = (
            np.zeros((Walker.wavefunction.num_mos[0],3), dtype=float),
            np.zeros((Walker.wavefunction.num_mos[1],3), dtype=float)
        )
        
        # Initialise alpha-spin electrons
        for ielec in range(Walker.wavefunction.num_mos[0]):
            self.pos[0][ielec,:] = place_electron()
        # Initialise beta-spin electrons
        for ielec in range(Walker.wavefunction.num_mos[1]):
            self.pos[1][ielec,:] = place_electron()

        self.matrix = Walker.wavefunction.matrix(self.pos)
        self.inverse = self.invert_wavefunction()
        self.dets = self.determinant_wavefunction()
        #self.local_energy()
            
    def __str__(self):
        """String representation of the Walker.
        """
        np.set_printoptions(formatter={'float': '{:7.4f}'.format})
        return_str = "Walker\n"
        return_str += "{0} alpha- and {1} beta-spin electrons\n".format(self.pos[0].shape[0], self.pos[1].shape[0])
        return_str += "Pot: {0:7.4f}Ha   Kin: {1:7.4f}Ha\n".format(self.potential(), self.kinetic())
        return_str += "Alpha Electron Positions\n"
        for ialpha in range(self.pos[0].shape[0]):
            return_str += "{0}\n".format(self.pos[0][ialpha,:]) 
        return_str += "Beta Electron Positions\n"
        for ibeta in range(self.pos[1].shape[0]):
            return_str += "{0}\n".format(self.pos[1][ibeta,:])
        return return_str
        
    def local_energy(self):
        """Compute the local energy of the walker's configuration.
        """
        return self.potential() + self.kinetic()
    
    def potential(self):
        """Compute the potential energy contribution to the walker's local energy.
        """
        def vne(self):
            """Compute the nuclear-electron attractive potential energy contribution to
            the walker's local energy.
            """
            vne_acc = 0.0
            electrons = self.pos[0] + self.pos[1]
            for iel in range(electrons.shape[0]):
                for atom in Walker.system.atoms:
                    rij = electrons[iel,:] - atom.pos
                    vne_acc -= atom.Z / np.sqrt(np.linalg.norm(rij))
            return vne_acc
                    
        def vee(self):
            """Compute the electron-electron repulsive potential energy contribution to
            the walker's local energy.
            """
            vee_acc = 0.0
            electrons = self.pos[0] + self.pos[1]
            for iel in range(electrons.shape[0]):
                for jel in range(iel+1,electrons.shape[0]):
                    rij = electrons[iel,:] - electrons[jel,:]
                    vee_acc += 1.0 / np.sqrt(np.linalg.norm(rij))
            return vee_acc
        
        return vne(self) + vee(self) + Walker.system.Vnn
    
    def kinetic(self):
        """Compute the total electron kinetic energy contribution to the walker's local energy.
        """
        laplacian_wfn = Walker.wavefunction.laplacian(self.pos)
        alpha_kin = np.einsum('ij,ij->i', laplacian_wfn[0], self.inverse[0].T)
        beta_kin  = np.einsum('ij,ij->i', laplacian_wfn[1], self.inverse[1].T)
        return np.sum(alpha_kin + beta_kin)
    
    def invert_wavefunction(self):
        return (np.linalg.inv(self.matrix[0]), np.linalg.inv(self.matrix[1]))

    def determinant_wavefunction(self):
        return (np.linalg.det(self.matrix[0]), np.linalg.det(self.matrix[1]))

    def wavefunction_value(self):
        return self.dets[0] * self.dets[1]
