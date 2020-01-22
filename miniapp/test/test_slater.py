import unittest
from numpy import array, transpose, zeros

# Plotting-based imports
from numpy import meshgrid, linspace
import matplotlib.pyplot as plt

from qmc.System import System
from qmc.Wavefunction import Slater

class TestSlater(unittest.TestCase):

    def test_1(self):
        """Build a STO-6G RHF wavefunction for H2.
        """
        print("\n")

        system_config = [('H', 0.0, 0.0, 0.0), ('H', 1.0, 0.0, 0.0)]
        system = System.System(system_config)

        mo_coeffs = (zeros((1,2)) + array([0.153, 0.153]), zeros((1,2)) + array([0.153, 0.153]))
        eval_pos = (zeros((1,3)) + array([1.0, 1.0, 1.0]), zeros((1,3)) + array([-1.0, -1.0, -1.0]))

        wfn = Slater.Slater(eval_pos, system, mo_coeffs)
        print(wfn)        

        # Want to scan between [-1, 2] along x-axis and [-1, 1] along y-axis
        nx = 200; ny = 200
        wfn_vals = zeros((nx,ny))
        xx, yy = meshgrid(linspace(-1, 2, nx), linspace(-1, 1, ny), sparse=False, indexing='ij')
        for ix in range(nx):
            for iy in range(ny):
                eval_pos = (zeros((1,3)) + array([-1 + ix*(3/nx), -1 + iy*(2/ny), 0]), zeros((1,3)) + array([0.0, 0.0, 0.0]))
                wfn.evaluate(eval_pos)
                wfn_vals[ix,iy] = wfn.value()
        
        fig, ax = plt.subplots()
        contours = ax.contour(xx, yy, wfn_vals)
        ax.clabel(contours, fontsize=9, inline=1)
        ax.set_xlabel("x / $\AA$")
        ax.set_ylabel("y / $\AA$")
        ax.set_title("H$_2$ STO-6G Wavefunction")
        ax.scatter(0.0, 0.0)
        ax.scatter(1.0, 0.0)
        ax.annotate("H", (0.0, 0.0))
        ax.annotate("H", (1.0, 0.0))
        plt.show()
        
if __name__ == "__main__":
    unittest.main()
