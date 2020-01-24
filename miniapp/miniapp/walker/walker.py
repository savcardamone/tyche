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
            random_charge = np.random.randint(low=0, high=system.tot_z)
            cumsum_charge = itertools.accumulate([atom.z for atom in Walker.system.atoms])
            atom_idx = np.where(filter(lambda x: x < random_charge, cumsum_charge))[0][-1]
            
            # Return a Gaussian variate centred on the chosen atom
            return np.random.normal(Walker.system.atoms[idx].pos, Walker.system.atoms[idx].Rvdw, size=3)
        
        # Give the walker a copy of the wavefunction
        self.wfn = wfn
        # Initialise the array to store the alpha and beta electrons
        self.pos = (np.zeros((self.wfn.num_mos[0],3)), np.zeros(self.wfn.num_mos[1],3))
        
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
        pass
    
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
            pass
        
        def vee(self):
            """Compute the electron-electron repulsive potential energy contribution to
            the walker's local energy.
            """
            pass
        
        return vne(self) + vee(self) + Walker.system.Vnn
    
    def kinetic(self):
        """Compute the total electron kinetic energy contribution to the walker's local energy.
        """
        pass
    
