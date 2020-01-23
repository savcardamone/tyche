
import sys
import numpy as np
import untangle

from miniapp.orbital import lcao

class Slater():
    """Slater wavefunction class. Keeps alpha and beta spin components of the
    wavefunction separate.
    """
    
    def __init__(self, system, input_file):
        """Class constructor.
        """            
        # Keep hold of the atomic orbitals that form the LCAO
        self.aos = lcao.LCAO(system, input_file)

        # Pull out the level of the xml we want to be looping through
        xml = untangle.parse(input_file)
        molecular_orbitals = xml.Input.Wavefunction.MolecularOrbitals

        self.num_mos = (int(molecular_orbitals['num_alpha']), int(molecular_orbitals['num_beta']))
        self.mo_coeffs = (
            np.zeros((self.num_mos[0], self.aos.num_aos), dtype=float),
            np.zeros((self.num_mos[1], self.aos.num_aos), dtype=float)
        )

        imo_a = 0; imo_b = 0
        for molecular_orbital in molecular_orbitals.Coefficients:
            coeffs = np.fromstring(molecular_orbital.cdata, dtype=float, sep=',')
            if molecular_orbital['spin'] == "Alpha" or molecular_orbital['spin'] == "Both":
                self.mo_coeffs[0][imo_a,:] = coeffs; imo_a += 1
            if molecular_orbital['spin'] == "Beta"  or molecular_orbital['spin'] == "Both":
                self.mo_coeffs[1][imo_b,:] = coeffs; imo_b += 1
                
        if imo_a != self.num_mos[0] or imo_b != self.num_mos[1]:
            sys.exit("Discrepancy in number of molecular orbitals in input file.")

        if self.num_mos[0] == 0:
            self.alpha = np.zeros((1, 1), dtype=float)
        else:
            self.alpha = np.zeros((self.num_mos[0], self.num_mos[0]), dtype=float)

        if self.num_mos[1] == 0:
            self.beta = np.zeros((1, 1), dtype=float)
        else:
            self.beta = np.zeros((self.num_mos[1], self.num_mos[1]), dtype=float)
 
        self.alpha_det = 1.0; self.beta_det = 1.0        

        
    def __str__(self):
        """String representation of the Slater wavefunction.
        """

        np.set_printoptions(formatter={'float': '{:7.4f}'.format})
        return "Slater Wavefunction\n" \
            " * Alpha MO Coefficients\n   {0}\n * Beta MO Coefficients\n   {1}\n" \
            " * Alpha Slater Matrix\n   {2}\n * Beta Slater Matrix\n   {3}\n" \
            " * Alpha Determinant : {4:7.4f}\n * Beta Determinant  : {5:7.4f}\n * Slater Determinant: {6:7.4f}"\
            .format(self.mo_coeffs[0], self.mo_coeffs[1], \
                    self.alpha, self.beta, self.alpha_det, self.beta_det, self.value())

    def value(self):
        """Return the product of spin-Slater determinants.
        """
        return self.alpha_det * self.beta_det
        
    def evaluate(self, pos):
        """Evaluate the spin-Slater matrices and determinants for all electrons in
        the system. Note that this is brute force determinant calculating, and not
        Sherman-Morrison, so cubic scaling in number of electrons.
        """
        alpha_elec_aos = np.apply_along_axis(self.aos.evaluate, axis=1, arr=pos[0])
        beta_elec_aos  = np.apply_along_axis(self.aos.evaluate, axis=1, arr=pos[1])

        self.alpha = np.dot(self.mo_coeffs[0], np.transpose(alpha_elec_aos))
        self.beta  = np.dot(self.mo_coeffs[1], np.transpose(beta_elec_aos))
        
        self.alpha_det = np.linalg.det(self.alpha); self.beta_det = np.linalg.det(self.beta)

    def update(self, pos, iel):
        """Sherman-Morrison update of the appropriate spin-Slater matrix and
        corresponding determinant.
        """
        print("Sherman-Morrison updating unsupported.")
