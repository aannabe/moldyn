#!/usr/bin/env python

import sys
sys.path.append('../moldyn/')
from base import *

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation, PillowWriter

###============ USER INPUT PARAMETERS ================

BOX_SIDE_LENGTH = 5.0   # Side legth of the simulation box
R_CUT = 4.0             # Cutoff radius beyond which interactions are ignored
PARTICLE_NUMBER = 25    # Total number of particles in the simulation
SPEED_AMPLITUDE = 12    # Initial speed amplitude of particles
TIMESTEP = 0.004        # Simulation timestep
TOTAL_TIME = 6.0        # Total time to propagate
COOL_FACTOR = 0.997     # Factor multiplying the velocities

###===================================================

ensemble = Ensemble(BOX_SIDE_LENGTH, PARTICLE_NUMBER, R_CUT)
ensemble.generate_particles()
ensemble.set_initial_velocities(SPEED_AMPLITUDE)
ensemble.normalize_momentum()
ensemble.update_ensemble_acceleration()

print('Initial Propagation ...')
for t in range(0, int(1.0/(TIMESTEP))):
    ensemble.propagate_velocity_verlet(TIMESTEP)
# end for
print('Initial Propagation Done.')

fig, ax = plt.subplots()

temperatures = []

def animate(t):
    ensemble.propagate_velocity_verlet(TIMESTEP)
    ensemble.change_temperature(COOL_FACTOR)
    #ensemble.normalize_momentum()
    ensemble.calculate_ensemble_energy()
    T = ensemble.get_ensemble_temperature()
    temperatures.append(T)
    ax.clear()
    ax.set_xlim(0, ensemble._BOX_SIDE_LENGTH)
    ax.set_ylim(0, ensemble._BOX_SIDE_LENGTH)
    plt.xticks([])
    plt.yticks([])
    plt.gca().set_aspect('equal')
    images = []
    for i in range(0, ensemble._PARTICLE_NUMBER):
        location = ensemble._particles[i].location
        line, = ax.plot(location[0], location[1], 'ro', markersize=8)
        images.append(line)
    running_avg = np.array(temperatures)[-200:].mean()
    text = plt.text(0.6, 5.2, 'Temperature = {:.2f}'.format(running_avg), fontsize = 20)
    images.append(text)
    # end for
    return images
# end animate
 
ani = FuncAnimation(fig, animate, interval=1, blit=True, repeat=True, frames=int(TOTAL_TIME/TIMESTEP))
ani.save("animation.gif", dpi=150, writer=PillowWriter(fps = 25))

