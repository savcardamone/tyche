#!/usr/bin/env python

"""test_lcao.py: Verify the LCAO class functions as it's meant to."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import os
import unittest as ut
import numpy as np

from test import test_base
from miniapp.system import system
from miniapp.orbital import lcao

class TestSystem(ut.TestCase):
    """Unit tests for the LCAO class.
    """

    def test_1(self):
        """Extract the atomic system from a set of test input files. Dump the xmls from
        the objects that get created then try to read them back into new objects. 
        Verify the objects are equivalent. 
        """
        test_files = test_base.grab_input_files()
        
        print("TODO: No numerical verification has been added to this yet.")
        for test_file in test_files:
            
            print("Working on input file: {0}".format(test_file))
            
            test_system = system.System("{0}".format(test_file))
            test_lcao = lcao.LCAO(test_system, "{0}".format(test_file))

            test_pos = np.array([1.0, 1.0, 1.0], dtype=float)
            test_lcao.evaluate(test_pos)
