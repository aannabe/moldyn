#!/usr/bin/env python

import numpy as np

class Particle:
    """Class representing a particle in a simulation.

    Parameters
    ----------
    location : list of float, optional
        Initial position of the particle. Default is [None, None].
    velocity : list of float, optional
        Initial velocity of the particle. Default is [None, None].
    acceleration : list of float, optional
        Initial acceleration of the particle. Default is [None, None].

    Attributes
    ----------
    _MASS : float
        Mass of the particle.
    location : numpy.ndarray
        Current position of the particle.
    velocity : numpy.ndarray
        Current velocity of the particle.
    acceleration : numpy.ndarray
        Current acceleration of the particle.

    Methods
    -------
    get_kinetic_energy()
        Calculate the kinetic energy of the particle.
    get_momentum()
        Calculate the momentum of the particle.
    """

    def __init__(self, location = [None, None], velocity = [None, None], acceleration = [None, None]):
        self._MASS = 1.0
        self.location = np.array(location)
        self.velocity = np.array(velocity)
        self.acceleration = np.array(acceleration)
    #end __init__

    def get_kinetic_energy(self):
        """Calculate the kinetic energy of the particle.

        Returns
        -------
        float
            The kinetic energy of the particle.
        """
        kinetic_energy = self._MASS * (self.velocity**2).sum() * 0.5 # KE = (mv^2)/2
        return kinetic_energy
    # end get_kinetic_energy

    def get_momentum(self):
        """Calculate the momentum of the particle.

        Returns
        -------
        numpy.ndarray
            The momentum of the particle.
        """
        momentum = self._MASS * self.velocity # M = mv
        return momentum
    # end get_momentum

# end Particle

class Interaction:
    """Class representing particle interactions.

    Parameters
    ----------
    R_CUT : float
        Cutoff radius for interactions.

    Attributes
    ----------
    __R_CUT : float
        Cutoff radius for interactions.

    Methods
    -------
    LJ_potential(r)
        Calculate the Lennard-Jones potential for a given distance.
    LJ_force(r)
        Calculate the Lennard-Jones force for a given distance.
    """

    def __init__(self, R_CUT):
        self.__R_CUT = R_CUT
    # end __init__

    def LJ_potential(self, r):
        """Calculate the Lennard-Jones potential for a given distance.

        Parameters
        ----------
        r : float
            Distance between particles.

        Returns
        -------
        float
            Lennard-Jones potential.
        """
        r = float(r)
        if r < self.__R_CUT:
            V = 4.0 * ((1.0/r)**12.0 - (1.0/r)**6.0)
        else:
            V = 0
        return V
    # end LJ_potential
    
    def LJ_force(self, r):
        """Calculate the Lennard-Jones force for a given distance.

        Parameters
        ----------
        r : float
            Distance between particles.

        Returns
        -------
        float
            Lennard-Jones force.
        """
        r = float(r)
        if r < self.__R_CUT:
            F = 4.0 * ((-12.0)*r**(-13.0) + 6.0*r**(-7.0))
        else:
            F = 0
        return F
    # end LJ_force

# end Interaction

class Box:
    """Class representing a simulation box.

    Parameters
    ----------
    BOX_SIDE_LENGTH : float
        Side length of the simulation box.
    R_CUT : float
        Cutoff radius for interactions.

    Attributes
    ----------
    _BOX_SIDE_LENGTH : float
        Side length of the simulation box.
    _interaction : Interaction
        Interaction object for particle interactions.

    Methods
    -------
    set_periodic_boundary(r)
        Apply periodic boundary conditions to a distance.
    """

    def __init__(self, BOX_SIDE_LENGTH, R_CUT):
        self._BOX_SIDE_LENGTH = BOX_SIDE_LENGTH
        self._interaction = Interaction(R_CUT) # using Composition for more flexible code ("has-a" relationship)
    # end def

    def set_periodic_boundary(self, r):
        """Apply periodic boundary conditions to a distance.

        Parameters
        ----------
        r : float
            Distance to apply periodic boundary conditions.

        Returns
        -------
        float
            Distance with periodic boundary conditions applied.
        """
        # Using minimum-image convention: only interact with the nearest copy of the particle (real or image)
        if (abs(r) > 0.5 * self._BOX_SIDE_LENGTH):
            r *= 1 - self._BOX_SIDE_LENGTH/abs(r)
        return r
    # end def
# end Box    

class Ensemble(Box): # Ensemble inherits from Box
    """Class representing an ensemble of particles within a simulation box.

    Parameters
    ----------
    BOX_SIDE_LENGTH : float
        Side length of the simulation box.
    PARTICLE_NUMBER : int
        Number of particles in the ensemble.
    R_CUT : float
        Cutoff radius for interactions.

    Attributes
    ----------
    _PARTICLE_NUMBER : int
        Number of particles in the ensemble.
    _particles : list of Particle
        List containing Particle objects.
    kinetic_energy : float
        Total kinetic energy of the ensemble.
    potential_energy : float
        Total potential energy of the ensemble.
    total_energy : float
        Total energy of the ensemble.
    temperature : float
        Temperature of the ensemble.

    Methods
    -------
    generate_particles()
        Generate particles and set their initial positions.
    set_initial_velocities(VELOCITY_AMPLITUDE)
        Set initial velocities for all particles in the ensemble.
    calculate_ensemble_momentum()
        Calculate the total momentum of the ensemble.
    calculate_ensemble_energy()
        Calculate the total energy of the ensemble.
    get_ensemble_temperature()
        Calculate the temperature of the ensemble.
    normalize_momentum()
        Normalize the momentum of the ensemble.
    update_ensemble_acceleration()
        Update the acceleration of each particle in the ensemble.
    propagate_velocity_verlet(dt)
        Propagate particle positions and velocities using the Velocity-Verlet algorithm.
    change_temperature(eta)
        Change the temperature of the ensemble.
    get_particle_positions()
        Get the positions of all particles in the ensemble.
    plot_ensemble()
        Plot the ensemble configuration.
    """

    def __init__(self, BOX_SIDE_LENGTH, PARTICLE_NUMBER, R_CUT): # initialize an empty box
        super().__init__(BOX_SIDE_LENGTH, R_CUT)
        self._PARTICLE_NUMBER = PARTICLE_NUMBER
        self._particles = []
    # end __init__

    def generate_particles(self): # set the positions
        """Generate particles and set their initial positions."""
        particle_root_1 = int(np.ceil(np.sqrt(self._PARTICLE_NUMBER)))          # 1st internal value to evenly place the particles
        particle_root_2 = int(np.ceil(self._PARTICLE_NUMBER / particle_root_1)) # 2nd internal value to evenly place the particles
        added_particles = 0
        for i in range(particle_root_1):
            for j in range(particle_root_2):
                Lx = (i + 0.5) * self._BOX_SIDE_LENGTH / particle_root_1
                Ly = (j + 0.5) * self._BOX_SIDE_LENGTH / particle_root_2
                if added_particles < self._PARTICLE_NUMBER:
                    self._particles.append(Particle([Lx, Ly])) # using Composition for more flexible code ("has-a" relationship)
                    added_particles += 1
                else:
                    return
        # end for
    # end generate_particles

    def set_initial_velocities(self, VELOCITY_AMPLITUDE):
        """Set initial velocities for all particles in the ensemble.

        Parameters
        ----------
        VELOCITY_AMPLITUDE : float
            Amplitude for generating initial velocities.
        """
        for i in range(0, self._PARTICLE_NUMBER):
            Vx = VELOCITY_AMPLITUDE * (2 * np.random.rand() - 1)
            Vy = VELOCITY_AMPLITUDE * (2 * np.random.rand() - 1)
            self._particles[i].velocity = np.array([Vx, Vy])
        # end for
    # end set_initial_velocities

    def calculate_ensemble_momentum(self):
        """Calculate the total momentum of the ensemble.

        Returns
        -------
        float
            Total momentum of the ensemble.
        """
        total_momentum = 0
        for i in range(0, self._PARTICLE_NUMBER):
            total_momentum += self._particles[i].get_momentum()
        # end for
        return total_momentum
    # end calculate_ensemble_momentum

    def calculate_ensemble_energy(self):
        """Calculate the total energy of the ensemble.

        Returns
        -------
        float
            Total energy of the ensemble.
        """
        total_kinetic = 0
        total_potential = 0
        for i in range(0, self._PARTICLE_NUMBER):
            for j in range(0, self._PARTICLE_NUMBER):
                if i==j: # avoid self-interaction
                    continue
                pair_distance = self._particles[j].location - self._particles[i].location
                pair_distance[0] = self.set_periodic_boundary(pair_distance[0])
                pair_distance[1] = self.set_periodic_boundary(pair_distance[1])
                distance_norm = np.linalg.norm(pair_distance)
                total_potential = (total_potential + self._interaction.LJ_potential(distance_norm)) / 2.0 # divide by two because it counts twice
            total_kinetic = total_kinetic + self._particles[i].get_kinetic_energy()
        # end for
        total_energy = total_kinetic + total_potential
        self.kinetic_energy = total_kinetic
        self.potential_energy = total_potential
        self.total_energy = total_energy
        return self.total_energy
    # end calculate_ensemble_energy

    def get_ensemble_temperature(self):
        """Calculate the temperature of the ensemble.

        Returns
        -------
        float
            Temperature of the ensemble.
        """
        # KE_avg = 3/2 * k * T
        self.temperature = 2/3 * self.kinetic_energy / self._PARTICLE_NUMBER
        return self.temperature
    # end get_ensemble_temperature

    def normalize_momentum(self):
        """Normalize the momentum of the ensemble."""
        total_momentum = self.calculate_ensemble_momentum()
        for i in range(0, self._PARTICLE_NUMBER):
            delta_V = total_momentum / (self._particles[i]._MASS * self._PARTICLE_NUMBER)
            self._particles[i].velocity -= delta_V
        # end for
    # end normalize_momentum

    def update_ensemble_acceleration(self):
        """Update the acceleration of each particle in the ensemble."""
        for i in range(0, self._PARTICLE_NUMBER):
            particle_acceleration = 0
            for j in range(0, self._PARTICLE_NUMBER):
                if i==j: # avoid self-interaction
                    continue
                pair_distance = self._particles[j].location - self._particles[i].location
                pair_distance[0] = self.set_periodic_boundary(pair_distance[0])
                pair_distance[1] = self.set_periodic_boundary(pair_distance[1])
                distance_norm = np.linalg.norm(pair_distance)
                particle_acceleration += self._interaction.LJ_force(distance_norm) * (pair_distance / distance_norm)
            self._particles[i].acceleration = particle_acceleration
    # end update_ensemble_acceleration

    def propagate_velocity_verlet(self, dt): # this is the Velocity-Verlet algorithm
        """Propagate particle positions and velocities using the Velocity-Verlet algorithm.

        Parameters
        ----------
        dt : float
            Time step for the integration.
        """
        #self.update_ensemble_acceleration()
        # update the positions and partly velocity
        for i in range(0, self._PARTICLE_NUMBER):
            self._particles[i].location = (self._particles[i].location + self._particles[i].velocity * dt + 0.5 * self._particles[i].acceleration * dt**2) % self._BOX_SIDE_LENGTH
            self._particles[i].velocity = self._particles[i].velocity + 0.5 * self._particles[i].acceleration * dt
        # end for
        # complete the update of velocities
        self.update_ensemble_acceleration()
        for i in range(0, self._PARTICLE_NUMBER):
            self._particles[i].velocity = self._particles[i].velocity + 0.5 * self._particles[i].acceleration * dt
        # end for
    # end propagate_velocity_verlet

    def change_temperature(self, eta):
        """Change the temperature of the ensemble.

        Parameters
        ----------
        eta : float
            Scaling factor for velocities.
        """
        for i in range(0, self._PARTICLE_NUMBER):
            self._particles[i].velocity = self._particles[i].velocity * eta
        # end for
    # end change_temperature

    def get_particle_positions(self):
        """Get the positions of all particles in the ensemble.

        Returns
        -------
        list of numpy.ndarray
            List of particle positions.
        """
        positions = []
        for i in range(0, self._PARTICLE_NUMBER):
            location = self._particles[i].location
            positions.append(location)
        # end for
        return positions
    # end get_particle_positions

    def plot_ensemble(self):
        """Plot the ensemble configuration."""
        fig, ax = plt.subplots()
        ax.set_xlim(0, self._BOX_SIDE_LENGTH)
        ax.set_ylim(0, self._BOX_SIDE_LENGTH)
        for i in range(0, self._PARTICLE_NUMBER):
            location = self._particles[i].location
            ax.plot(location[0], location[1], 'ro')
        # end for
        plt.xticks([])
        plt.yticks([])
        plt.gca().set_aspect('equal')
        #plt.pause(0.001)
        #ax.cla()
        plt.show()
    # end plot_ensemble

# end Ensemble

