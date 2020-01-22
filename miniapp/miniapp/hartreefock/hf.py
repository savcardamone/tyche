from math import pi
from numpy import array, ndarray, divide, sqrt, argsort, sort, diag, trace
from numpy.linalg import eig, norm

class HartreeFock():

    zeta = array([38.474970, 5.782948, 1.242567, 0.298073])
    num_aos = len(zeta)
    num_mos = 0
    energy_tolerance = 0.0001; density_tolerance = 0.001
    prev_energy = 0
    prev_density = []

    def __init__(self, num_elec):

        # Make sure we can pair electrons
        if num_elec % 2 != 0:
            raise Exception("Can't do a RHF with", num_elec, "electrons.")
        else:
            print("Restricted Hartree-Fock with", num_elec, "electron(s).")

        # We're RHF, so pair up spins in each molecular orbital
        self.num_mos = int(num_elec / 2)

        if self.num_mos > self.num_aos:
            raise Exception("Can't create", self.num_mos, "molecular orbital(s) from", self.num_aos, "atomic orbital(s).")
        else:
            print(self.num_aos, "atomic orbital(s) and", self.num_mos, "molecular orbital(s).")
            print("Zeta: ", self.zeta)

        self.prev_density = ndarray(shape=(self.num_aos,self.num_aos),dtype=float, order='C')

    def one_electron_integrals(self):

        def overlap_kernel(zeta_i, zeta_j):
            return pow(pi / (zeta_i + zeta_j), 1.5)

        def kinetic_kernel(zeta_i, zeta_j):
            return 3 * pow(pi, 1.5) * (zeta_i * zeta_j) / pow(zeta_i + zeta_j, 2.5)

        def nucattr_kernel(zeta_i, zeta_j):
            return (-4 * pi) / (zeta_i + zeta_j)

        # Initialise our matrices
        overlap = ndarray(shape=(self.num_aos,self.num_aos), dtype=float, order='C')
        kinetic = ndarray(shape=(self.num_aos,self.num_aos), dtype=float, order='C')
        nucattr = ndarray(shape=(self.num_aos,self.num_aos), dtype=float, order='C')

        for i_ao in range(self.num_aos):
            for j_ao in range(self.num_aos):
                overlap[i_ao,j_ao] = overlap_kernel(self.zeta[i_ao], self.zeta[j_ao])
                kinetic[i_ao,j_ao] = kinetic_kernel(self.zeta[i_ao], self.zeta[j_ao])
                nucattr[i_ao,j_ao] = nucattr_kernel(self.zeta[i_ao], self.zeta[j_ao])

        return overlap, kinetic, nucattr

    def two_electron_integrals(self):

        def tei_kernel(zeta_i, zeta_j, zeta_k, zeta_l):
            temp_1 = (zeta_i + zeta_j) * (zeta_k + zeta_l)
            temp_2 = sqrt(zeta_i + zeta_j + zeta_k + zeta_l)
            return 2 * pow(pi, 2.5) / (temp_1 * temp_2)

        teis = ndarray(shape=(self.num_aos,self.num_aos,self.num_aos,self.num_aos), dtype=float, order='C')

        for i_ao in range(self.num_aos):
            for j_ao in range(self.num_aos):
                for k_ao in range(self.num_aos):
                    for l_ao in range(self.num_aos):
                        teis[i_ao,j_ao,k_ao,l_ao] = tei_kernel(self.zeta[i_ao], self.zeta[j_ao], self.zeta[k_ao], self.zeta[l_ao])

        return teis

    def basis_transformation_matrix(self, overlap):

        # Get the eigenvalues and eigenvectors of the overlap matrix
        overlap_evals, overlap_evecs = eig(overlap)

        # Create diagonal matrix with entries given by inverse of eigenvalues of
        # overlap matrix
        try:
            inv_sqrt_evals = diag(divide(1., sqrt(overlap_evals)))
        except:
            raise Exception("Overlap matrix is not positive definite.")

        # Construct the basis transformation matrix and return it
        return overlap_evecs @ inv_sqrt_evals @ overlap_evecs.T

    def fock_matrix(self, core_hamiltonian, teis, density):

        fock = ndarray(shape=density.shape, dtype=float, order='C')

        for i_ao in range(self.num_aos):
            for j_ao in range(self.num_aos):

                fock[i_ao,j_ao] = core_hamiltonian[i_ao,j_ao]

                for k_ao in range(self.num_aos):
                    for l_ao in range(self.num_aos):

                        coulomb = teis[i_ao,k_ao,j_ao,l_ao]
                        exchange = teis[i_ao,k_ao,l_ao,j_ao]

                        fock[i_ao,j_ao] += density[k_ao,l_ao] * (coulomb - 0.5*exchange)

        return fock

    def density_matrix(self, overlap, basis_transform, fock):

        def ordered_eigensystem(matrix):

            # Generate the eigenvalues and eigenvectors of the matrix
            evals, evecs = eig(matrix)

            # Sort the eigenvalues in ascending order and keep a track of what index they
            # were originally assigned
            ordered_indices = argsort(evals)
            ordered_evals = sort(evals)

            # Order the eigenvectors in asceding order of their corresponding eigenvalues
            ordered_evecs = ndarray(shape=evecs.shape, dtype=float, order='C')
            ordered_transform = ndarray(shape=evecs.shape, dtype=float, order='C')
            for i_evec in range(len(ordered_evals)):
                ordered_evecs[:,i_evec] = evecs[:,ordered_indices[i_evec]]
                ordered_transform[i_evec,:] = basis_transform[ordered_indices[i_evec],:]

            # Return the ordered eigenvalues and corresponding eigenvectors
            return ordered_evals, ordered_evecs, ordered_transform

        # Transform Fock matrix to orthogonal basis
        fock = basis_transform.T @ fock @ basis_transform

        # Get the eigenvalues and eigenvectors of the input Fock matrix
        fock_evals, fock_evecs, new_transform = ordered_eigensystem(fock)
        # Transform the eigenvectors of the Fock matrix back to the original basis
        fock_evecs = new_transform @ fock_evecs

        # First of all we make sure the eigenvectors of the Fock matrix are normalised by the 
        # overlap matrix (these are molecular orbitals, afterall)
        for i_mo in range(self.num_aos):

            ao_coeffs = fock_evecs[:,i_mo]
            norm = ao_coeffs.T @ overlap @ ao_coeffs
            fock_evecs[:,i_mo] /= sqrt(norm)

        # Initialise the density matrix
        density = ndarray(shape=overlap.shape, dtype=float, order='C')

        # Loop over all elements in the density matrix and accumulate
        for i_ao in range(self.num_aos):
            for j_ao in range(self.num_aos):

                density[i_ao,j_ao] = 0.0

                # We accumulate only over occupied molecular orbitals! Note that we also have
                # access to the virtual orbitals at this point, but they're effectively discarded
                for i_mo in range(self.num_mos):
                    density[i_ao,j_ao] += 2 * fock_evecs[i_ao,i_mo] * fock_evecs[j_ao,i_mo]

        return fock_evecs, density

    def scf_energy(self, density, core_hamiltonian, fock):

        energy = 0.0

        for i_ao in range(self.num_aos):
            for j_ao in range(self.num_aos):
                energy += 0.5 * density[i_ao,j_ao] * (core_hamiltonian[i_ao,j_ao] + fock[i_ao,j_ao])

        return energy

    def check_convergence(self, energy, density):

        if abs(energy - self.prev_energy) < self.energy_tolerance:
            energy_converged = True
        else:
            energy_converged = False
        self.prev_energy = energy

        if norm(density - self.prev_density) < self.density_tolerance:
            density_converged = True
        else:
            density_converged = False
        self.prev_density = density

        return energy_converged, density_converged

    def mulliken(self, overlap, density):

        return trace(density @ overlap)

    def run(self, num_cycles):

        print("Hartree-Fock will run for a maximum of", num_cycles, "SCF iteration(s).")

        overlap, kinetic, nucattr = self.one_electron_integrals()
        core_hamiltonian = kinetic + nucattr
        teis = self.two_electron_integrals()

        basis_transform = self.basis_transformation_matrix(overlap)
        _, density = self.density_matrix(overlap, basis_transform, core_hamiltonian)

        energy = self.scf_energy(density, core_hamiltonian, core_hamiltonian)

        for i in range(num_cycles):

            fock = self.fock_matrix(core_hamiltonian, teis, density)
            fock_evecs, density = self.density_matrix(overlap, basis_transform, fock)

            energy = self.scf_energy(density, core_hamiltonian, fock)
            print("Iteration", i, "SCF Energy:", energy)

            energy_converged, density_converged = self.check_convergence(energy, density)

            if energy_converged and density_converged:
                print("SCF has converged!")
                for i_mo in range(self.num_mos):
                    print("Molecular Orbital", i_mo, "Coefficients :", fock_evecs[:,i_mo])
                print("Mulliken charge:", self.mulliken(overlap, density))
                break

            if i == num_cycles - 1:
                print("SCF failed to converge.")
                print("Energy Convergence Check:", energy_converged)
                print("Density Convergence Check:", density_converged)


        fock_mo_basis = ndarray(shape=(self.num_mos,self.num_mos), dtype=float, order='C')

        for i_mo in range(self.num_mos):
            for j_mo in range(self.num_mos):

                fock_mo_basis[i_mo,j_mo] = 0.0

                for i_ao in range(self.num_aos):
                    for j_ao in range(self.num_aos):

                        fock_mo_basis[i_mo,j_mo] += fock_evecs[i_ao,j_mo] * fock_evecs[j_ao,i_mo] * fock[i_ao,j_ao]

        print(fock_mo_basis)

if __name__ == "__main__":

    hf = HartreeFock(4)
    hf.run(2000)