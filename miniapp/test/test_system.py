#!/usr/bin/env python

"""test_system.py: Verify the System class functions as it's meant to."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import os
import unittest as ut
import numpy as np

from miniapp.system import system, atom

class TestSystem(ut.TestCase):
    """Unit tests for the System class.
    """

    def test_1(self):
        """Extract the atomic system from a set of test input files.
        """
        input_dir = "{0}/input/".format(os.path.dirname(os.path.realpath(__file__)))
        test_files = ["{0}{1}".format(input_dir, f) for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]

        for test_file in test_files:
            print("Working on input file: {0}".format(test_file))
            test_system = system.System("{0}".format(test_file))
