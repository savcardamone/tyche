import unittest as ut
import numpy as np

from miniapp.orbital import atomic_orbital as ao

class TestAtomicOrbital(ut.TestCase):

    def test_1(self):
        """Build a STO-6G RHF wavefunction for H2.
        """
        test_centre  = np.array([1.0, 1.0, 1.0])
        test_coeffs  = np.array([0.1, 0.2, 0.3])
        test_zetas   = np.array([8.0, 4.0, 2.0])
        test_ang_mom = np.array([1,1,0])
        
        test_ao = ao.AtomicOrbital(test_centre, test_coeffs, test_zetas, test_ang_mom)

