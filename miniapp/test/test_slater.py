#!/usr/bin/env python

"""test_slater.py: Verify the Slater class functions as it's meant to."""

__author__      = "Salvatore Cardamone"
__email__       = "sav.cardamone@gmail.com"

import numpy as np
import unittest as ut
from scipy import integrate
from matplotlib import pyplot as plt

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

            test_walker = (
                np.random.uniform(-10.0, 10.0, size=(test_slater.num_mos[0],3)),
                np.random.uniform(-10.0, 10.0, size=(test_slater.num_mos[1],3))
            )

            (alpha_mat,beta_mat) = test_slater.matrix(test_walker)
            (alpha_lap,beta_lap) = test_slater.laplacian(test_walker)

    def test_2(self):
        """Evalaute the electron density for each of the test files. 

        Do numerical quadrature over all space and verify we recover the number of electrons.
        """
        test_files = test_base.grab_input_files()
        
        for test_file in test_files:
            
            print("Working on input file: {0}".format(test_file))
            
            test_system = system.System("{0}".format(test_file))
            test_slater = slater.Slater(test_system, "{0}".format(test_file))
            """
            ys = np.linspace(-10, 10, num=1000)
            zs = np.linspace(-10, 10, num=1000)
            grid = np.zeros(shape=(1000,1000), dtype=float)
            for iy in range(1000):
                for iz in range(1000):
                    grid[iy,iz] = test_slater.density(np.array([0.0, ys[iy], zs[iz]]))
            plt.contour(ys, zs, grid)
            plt.show()
            """

            # Full integral over R^3 to get the number of electrons a la Born (i.e. integral of square of
            # wavefunction over all space)
            density_sample = lambda z, y, x: test_slater.density(np.array([x, y, z]))
            density_integral = integrate.tplquad(
                density_sample, -10, 10, lambda x: -10, lambda x: 10, lambda x, y: -10, lambda x, y: 10
            )
            print("Expected Integral: {0}   Quadrature Result: {1} +/- {2}".format(test_slater.num_mos[0] + test_slater.num_mos[1], density_sample[0], density_sample[1]))
            # Verify we get the number of electrons from the integral, to within the quadrature error
            self.assertAlmostEqual(
                density_integral[0], test_slater.num_mos[0] + test_slater.num_mos[1], delta=density_integral[1]
            )
