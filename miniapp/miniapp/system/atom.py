#!/usr/bin/env python

"""system.py: Atom container and associated manipulation methods."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import sys
import numpy as np

# Dictionary of nuclear charges ascribed to atom types
nuclear_charge = {
    'H' : 1.0, 'He' : 2.0,
    'Li' : 3.0, 'Be' : 4.0, 'B' : 5.0, 'C' : 6.0, 'N' : 7.0, 'O' : 8.0, 'F' : 9.0, 'Ne' : 10.0
}
# Dictionary of van der Waals radii ascribed to atom types
vdw_radius = {
    'H' : 1.2, 'He' : 1.4,
    'Li' : 1.82, 'Be' : 1.53, 'B' : 1.92, 'C' : 1.70, 'N' : 1.55, 'O' : 1.52, 'F' : 1.47, 'Ne' : 1.54
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
            self.Rvdw = vdw_radius[self.atom_type]
        except KeyError:
            sys.exit("Couldn't locate atom type {0} during parameter search".format(atom_type))
            
    def __str__(self):
        """Object string representation.
        """
        np.set_printoptions(formatter={'float': '{:7.4f}'.format})
        return "Atom {0:2s} (Z = {1}, Rvdw = {2}):  {3}".format(self.atom_type, self.Z, self.Rvdw, self.pos)
