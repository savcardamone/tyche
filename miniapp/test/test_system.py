#!/usr/bin/env python

"""test_system.py: Verify the System class functions as it's meant to."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import os
import unittest as ut
import numpy as np
import tempfile

from test import test_base
from miniapp.system import system, atom

class TestSystem(ut.TestCase):
    """Unit tests for the System class.
    """

    def test_1(self):
        """Extract the atomic system from a set of test input files. Dump the xmls from
        the objects that get created then try to read them back into new objects. 
        Verify the objects are equivalent. 
        """
        test_files = test_base.grab_input_files()
        
        for test_file in test_files:
            
            print("Working on input file: {0}".format(test_file))
            test_system = system.System("{0}".format(test_file))

            # Dump the System object to an xml file
            test_system_xml_string = test_system.as_xml()
            # Dump the xml into /tmp
            temp_file, temp_filename = tempfile.mkstemp()
            with os.fdopen(temp_file, 'w') as tmp:
                tmp.write("<Input>")
                tmp.write(str(test_system_xml_string, 'utf-8'))
                tmp.write("</Input>")
            # Read the xml file we've just dumped into another System object
            retest_system = system.System(temp_filename)
            
            # Verify System object equality
            self.assertEqual(test_system.num_atoms, retest_system.num_atoms)
            for a, b in zip(test_system.atoms, retest_system.atoms):
                self.assertEqual(a.atom_type, b.atom_type)
                self.assertEqual(a.pos.tolist(), b.pos.tolist())
                self.assertEqual(a.Z, b.Z)
