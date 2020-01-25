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
    system = None
    
    @classmethod
    def initialise_parameters(cls, system):
        """Initialisation of class variables.
        We're going to need to system to evaluate nuclear contributions to the local energy,
        but this is constant data, so no need to replicate for each walker.
        """
        # Keep hold of the atomic orbitals that form the LCAO
        cls.system = system
        
    def __init__(self, wfn):
        """Class constructor.
        """
        # Verify that we've called the class variable initialisation routine before we start
        # trying to do object creation
        if Walker.system is None:
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
        
        # Give the walker a copy of the wavefunction
        self.wfn = wfn
        # Initialise the array to store the alpha and beta electrons
        self.pos = (
            np.zeros((self.wfn.num_mos[0],3), dtype=float),
            np.zeros((self.wfn.num_mos[1],3), dtype=float)
        )
        
        # Initialise alpha-spin electrons
        for ielec in range(self.wfn.num_mos[0]):
            self.pos[0][ielec,:] = place_electron()
        # Initialise beta-spin electrons
        for ielec in range(self.wfn.num_mos[1]):
            self.pos[1][ielec,:] = place_electron()
            
        # Evaluate the wavefunction for the walker's configuration
        self.wfn.evaluate(self.pos)
            
    def __str__(self):
        """String representation of the Walker.
        """
        np.set_printoptions(formatter={'float': '{:7.4f}'.format})
        return_str = ""
        return_str += "Walker with {0} alpha- and {1} beta-spin electrons\n".format(self.pos[0].shape[0], self.pos[1].shape[0])
        return_str += "Potential Energy: {0:7.4f} Ha   Kinetic Energy: {1:7.4f} Ha\n".format(self.potential(), self.kinetic())
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
        laplacian_wfn = self.wfn.aos.laplacian(self, self.pos)
        return 0.0
    
