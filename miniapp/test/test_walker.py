#!/usr/bin/env python

"""test_walker.py: Verify the Walker class functions as it's meant to."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import numpy as np
import unittest as ut

from test import test_base
from miniapp.system import system
from miniapp.wavefunction import slater
from miniapp.walker import walker

class TestWalker(ut.TestCase):

    def test_1(self):
        """Build a Walker for each of the test files.
        """
        test_files = test_base.grab_input_files()
        
        print("TODO: No numerical verification has been added to this yet.")
        for test_file in test_files:
            
            print("Working on input file: {0}".format(test_file))
            
            test_system = system.System("{0}".format(test_file))
            test_slater = slater.Slater(test_system, "{0}".format(test_file))

            walker.Walker.initialise_parameters(test_system, test_slater)
            test_walker = walker.Walker()
