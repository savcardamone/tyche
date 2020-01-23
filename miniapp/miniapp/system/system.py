#!/usr/bin/env python

"""system.py: Atomic system description."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import numpy as np
import untangle
import xml.etree.cElementTree as et

from miniapp.system import atom

class System():
    """Container for an atomic system.
    """
    
    def __init__(self, input_file):
        """Class constructor.
        config input parameter is a list of tuples, each comprising atomic element and the
        atom's position vector.
        """
        
        system_root = untangle.parse(input_file).Input.System
        self.num_atoms = int(system_root['num_atoms'])

        self.atoms = []
        for xml_atom in system_root.Atom:
            self.atoms.append(atom.Atom(xml_atom['type'], np.fromstring(xml_atom['pos'], dtype=float, sep=',')))
            
        if not self.atoms:
            sys.exit("No atoms were specified in the input file.")

        # Compute the nuclear repulsion energy arising in the system
        self.Vnn = self.nuclear_repulsion()

    def __str__(self):
        """Object string representation.
        """

        sys_string = "Atomic system comprising {0} atom(s)\n".format(self.num_atoms)
        for iatom in range(self.num_atoms):
            sys_string += "   {0}\n".format(self.atoms[iatom].__str__())
        sys_string += "   Nuclear repulsion energy: {0:7.4f} a.u.".format(self.Vnn)

        return sys_string

    def as_xml(self):
        """Dump the System object to an output file in xml format. System will be to root.
        """
        sys_xml = et.Element("System", num_atoms="{0}".format(self.num_atoms))
        for atom in self.atoms:
            atom_xml = et.SubElement(sys_xml, "Atom", type="{0}".format(atom.atom_type), pos="{0}".format(atom.pos.tolist()).strip("[]"))

        return et.tostring(sys_xml)
            
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
