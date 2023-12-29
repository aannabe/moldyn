#!/usr/bin/env python

import numpy as np

class Particle:
    def __init__(self, location = [None, None], velocity = [None, None], acceleration = [None, None]):
        self._MASS = 1.0
        self.location = np.array(location)
        self.velocity = np.array(velocity)
        self.acceleration = np.array(acceleration)
    #end __init__

    def get_kinetic_energy(self):
        kinetic_energy = self._MASS * (self.velocity**2).sum() * 0.5 # KE = (mv^2)/2
        return kinetic_energy
    # end get_kinetic_energy

    def get_momentum(self):
        momentum = self._MASS * self.velocity # M = mv
        return momentum
    # end get_momentum

# end Particle

class Interaction:
    def __init__(self, R_CUT):
        self.__R_CUT = R_CUT
    # end __init__

    def LJ_potential(self, r):
        r = float(r)
        if r < self.__R_CUT:
            V = 4.0 * ((1.0/r)**12.0 - (1.0/r)**6.0)
        else:
            V = 0
        return V
    # end LJ_potential
    
    def LJ_force(self, r):
        r = float(r)
        if r < self.__R_CUT:
            F = 4.0 * ((-12.0)*r**(-13.0) + 6.0*r**(-7.0))
        else:
            F = 0
        return F
    # end LJ_force

# end Interaction

class Box:
    def __init__(self, BOX_SIDE_LENGTH, R_CUT):
        self._BOX_SIDE_LENGTH = BOX_SIDE_LENGTH
        self._interaction = Interaction(R_CUT) # using Composition for more flexible code ("has-a" relationship)
    # end def

    def set_periodic_boundary(self, r):
        # Using minimum-image convention: only interact with the nearest copy of the particle (real or image)
        if (abs(r) > 0.5 * self._BOX_SIDE_LENGTH):
            r *= 1 - self._BOX_SIDE_LENGTH/abs(r)
        return r
    # end def
# end Box    

class Ensemble(Box): # Ensemble inherits from Box
    def __init__(self, BOX_SIDE_LENGTH, PARTICLE_NUMBER, R_CUT): # initialize an empty box
        super().__init__(BOX_SIDE_LENGTH, R_CUT)
        self._PARTICLE_NUMBER = PARTICLE_NUMBER
        self._particles = []
    # end __init__

    def generate_particles(self): # set the positions
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
        for i in range(0, self._PARTICLE_NUMBER):
            Vx = VELOCITY_AMPLITUDE * (2 * np.random.rand() - 1)
            Vy = VELOCITY_AMPLITUDE * (2 * np.random.rand() - 1)
            self._particles[i].velocity = np.array([Vx, Vy])
        # end for
    # end set_initial_velocities

    def calculate_ensemble_momentum(self):
        total_momentum = 0
        for i in range(0, self._PARTICLE_NUMBER):
            total_momentum += self._particles[i].get_momentum()
        # end for
        return total_momentum
    # end calculate_ensemble_momentum

    def calculate_ensemble_energy(self):
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
        # KE_avg = 3/2 * k * T
        self.temperature = 2/3 * self.kinetic_energy / self._PARTICLE_NUMBER
        return self.temperature
    # end get_ensemble_temperature

    def normalize_momentum(self):
        total_momentum = self.calculate_ensemble_momentum()
        for i in range(0, self._PARTICLE_NUMBER):
            delta_V = total_momentum / (self._particles[i]._MASS * self._PARTICLE_NUMBER)
            self._particles[i].velocity -= delta_V
        # end for
    # end normalize_momentum

    def update_ensemble_acceleration(self):
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
        for i in range(0, self._PARTICLE_NUMBER):
            self._particles[i].velocity = self._particles[i].velocity * eta
        # end for
    # end change_temperature

    def get_particle_positions(self):
        positions = []
        for i in range(0, self._PARTICLE_NUMBER):
            location = self._particles[i].location
            positions.append(location)
        # end for
        return positions
    # end get_particle_positions

    def plot_ensemble(self):
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

