#!/usr/bin/env python

"""system.py: Atomic system description."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import numpy as np

class Atom():
    """Container for an atom. Largely a wrapper around parameters.
    """
    
    # Dictionary or nuclear charges ascribed to atom types
    nuclear_charge = {
        'H' : 1.0, 'He' : 2.0,
        'Li' : 3.0, 'Be' : 4.0, 'B' : 5.0, 'C' : 6.0, 'N' : 7.0, 'O' : 8.0, 'F' : 9.0, 'Ne' : 10.0
    }
    
    def __init__(self, config):
        """Class constructor.
        config input parameter is a tuple comprising atomic element and the
        atom's position vector.
        """
        
        self.atom_type = config[0]
        self.pos = np.zeros((1,3)) + np.array([config[1], config[2], config[3]])
        self.Z = Atom.nuclear_charge[self.atom_type]
        
    def __str__(self):
        """Object string representation.
        """
        
        np.set_printoptions(formatter={'float': '{:7.4f}'.format})
        return "Atom {0:2s} (Z = {1}):  {2}".format(self.atom_type, self.Z, self.pos)

class System():
    """Container for an atomic system.
    """
    
    def __init__(self, config=None):
        """Class constructor.
        config input parameter is a list of tuples, each comprising atomic element and the
        atom's position vector.
        """
        
        def nuclear_repulsion(self):
            """Compute the nuclear-nuclear repulsion energy arising in the system.
            """
            
            nuc_rep = 0.0
            # Just loop over all unique pairs of atoms and accumulate
            for iatom in range(self.num_atoms):
                for jatom in range(iatom+1, self.num_atoms):
                    Rij = self.atoms[iatom].pos - self.atoms[jatom].pos
                    nuc_rep += self.atoms[iatom].Z * self.atoms[jatom].Z / np.sqrt(np.linalg.norm(Rij))

            return nuc_rep
        
        if config == None:
            print("Reading system configuration from file unsupported.")
        # Read the system from the config argument and explicitly construct a list of
        # atom objects
        else:
            self.num_atoms = len(config)
            self.tot_z = 0
            self.atoms = []
            for iatom in range(self.num_atoms):
                self.atoms.append(Atom(config[iatom]))
                self.tot_z += self.atoms[iatom].Z
                
        # Compute the nuclear repulsion energy arising in the system
        self.Vnn = nuclear_repulsion(self)
                
    def __str__(self):
        """Object string representation.
        """

        sys_string = "Atomic system comprising {0} atom(s)\n".format(self.num_atoms)
        for iatom in range(self.num_atoms):
            sys_string += "   {0}\n".format(self.atoms[iatom].__str__())
        sys_string += "   Nuclear repulsion energy: {0:7.4f} a.u.".format(self.Vnn)

        return sys_string

