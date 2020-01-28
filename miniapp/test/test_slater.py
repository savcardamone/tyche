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
            print(test_slater)

            test_walker = (
                np.random.uniform(-10.0, 10.0, size=(test_slater.num_mos[0],3)),
                np.random.uniform(-10.0, 10.0, size=(test_slater.num_mos[1],3))
            )

            (alpha_mat,beta_mat) = test_slater.matrix(test_walker)
            (alpha_lap,beta_lap) = test_slater.laplacian(test_walker)

            np.set_printoptions(formatter={'float': '{:11.4f}'.format})
            print("Alpha Electrons:\n{0}".format(test_walker[0]))
            print("Beta  Electrons:\n{0}".format(test_walker[1]))
            print("Alpha Slater Matrix:\n{0}".format(alpha_mat))
            print("Beta  Slater Matrix:\n{0}".format(beta_mat))
            print("Alpha Laplacian Matrix:\n{0}".format(alpha_lap))
            print("Beta  Laplacian Matrix:\n{0}".format(beta_lap))
            
