#!/usr/bin/env python

"""test_atom.py: Verify the Atom class functions as it's meant to."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import unittest as ut
import numpy as np

from miniapp.system import atom

class TestAtom(ut.TestCase):
    """Unit tests for the Atom class.
    """

    def test_1(self):
        """Make a Neon atom and verify the object construction.
        """
        atom_type = "Ne"
        position_vector = np.array([1.0, 1.0, 1.0])

        test_atom = atom.Atom(atom_type, position_vector)

        self.assertEqual(test_atom.atom_type, atom_type)
        self.assertEqual(test_atom.pos.tolist(), position_vector.tolist())
        self.assertEqual(test_atom.Z, atom.nuclear_charge[atom_type])
