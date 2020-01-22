from numpy.random import randint, normal, zeros
from numpy.linalg import norm

class Walker():

    def __init__(self):
        """Class constructor.
        """
        
        def place_electron(system):
            """Place an electron within the van der Waals radius of a random atom,
            weighted by its nuclear charge.
            """
            random_centre = randint(low=0, high=system.tot_z)
            acc_z = 0; selected_centre = 0
            for iatom in range(system.num_atoms):
                acc_z += system.atoms[iatom].Z
                if random_centre <= acc_z:
                    selected_centre = iatom
                    break

            # Return a Gaussian variate centred on the chosen atom
            return normal(system.atoms[selected_centre].pos, 2.0, size=3)

        # Initialise the array to store the alpha and beta electrons
        self.pos = (zeros((wfn.num_alpha,3)), zeros(wfn.num_beta,3))

        # Initialise alpha-spin electrons
        for ialpha in range(wfn.num_alpha):
            self.pos[0][ialpha,:] = place_electron(system)
        # Initialise beta-spin electrons
        for ibeta in range(wfn.num_beta):
            self.pos[0][ibeta,:] = place_electron(system)
        
    def __str__(self):



    def electron_nuclear_attraction(self, system):
        """Compute the nuclear-electron attractive potential energy contribution to
        the walker's local energy.
        """

        # Zero our accumulator for the electron-nuclear energy
        vne = 0.0
        
        for iatom in range(system.num_atoms):
            for ialpha in range(self.pos[0].shape[0]):
                Rij = self.pos[0][ialpha,:] - system.atoms[iatom].pos
                vne -= system.atoms[iatom].Z / sqrt(norm(Rij))
                
            for ibeta in range(self.pos[1].shape[0]):
                Rij = self.pos[0][ibeta,:] - system.atoms[iatom].pos
                vne -= system.atoms[iatom].Z / sqrt(norm(Rij))

        return vne
                
    def electron_electron_repulsion(self):

        vee = 0.0
        
         =
        for
