#!/usr/bin/env python

"""test_base.py: Some base functionality that can be used by all unit tests."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import os

def grab_input_files():
    """Search the input/ directory for use cases and return a list of all pertinent
    ones we can unit test with.
    """
    
    input_dir = "{0}/input/".format(os.path.dirname(os.path.realpath(__file__)))
    test_files = ["{0}{1}".format(input_dir, f) for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]

    return test_files
