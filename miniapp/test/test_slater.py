#!/usr/bin/env python

"""test_slater.py: Verify the Slater class functions as it's meant to."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import numpy as np
import unittest as ut

from test import test_base
from miniapp.system import system
from miniapp.wavefunction import slater

class TestSlater(ut.TestCase):

    def test_1(self):
        """Build Slater wavefunctions for each of the test files.
        """
        test_files = test_base.grab_input_files()
        
        print("TODO: No numerical verification has been added to this yet.")
        for test_file in test_files:
            
            print("Working on input file: {0}".format(test_file))
            
            test_system = system.System("{0}".format(test_file))
            test_slater = slater.Slater(test_system, "{0}".format(test_file))

            test_walker = (np.array([[1.0, 1.0, 1.0]]), np.array([[1.0, 1.0, 1.0]]))
            test_slater.evaluate(test_walker)
            print(test_slater)
            
