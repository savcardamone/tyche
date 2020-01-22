#!/usr/bin/env python

"""system.py: Atom container and associated manipulation methods."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import sys
import numpy as np

# Dictionary or nuclear charges ascribed to atom types
nuclear_charge = {
    'H' : 1.0, 'He' : 2.0,
    'Li' : 3.0, 'Be' : 4.0, 'B' : 5.0, 'C' : 6.0, 'N' : 7.0, 'O' : 8.0, 'F' : 9.0, 'Ne' : 10.0
}

class Atom():
    """Container for an atom. Largely a wrapper around parameters.
    Once we add some dynamics we'll need to expand this class -- or maybe factor it out completely
    and just have everything in System.
    """
    
    def __init__(self, atom_type, position):
        """Class constructor.
        Takes the atom type and the position vector of the atom.
        """ 
        self.atom_type = atom_type
        self.pos = position
        try:
            self.Z = nuclear_charge[self.atom_type]
        except KeyError:
            sys.exit("Couldn't find nuclear charge of atom type {0}".format(atom_type))
            
    def __str__(self):
        """Object string representation.
        """
        np.set_printoptions(formatter={'float': '{:7.4f}'.format})
        return "Atom {0:2s} (Z = {1}):  {2}".format(self.atom_type, self.Z, self.pos)
