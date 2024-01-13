Introduction
============

The key user input parameters and usage are shown in the below code snippet:

.. code-block:: python

    BOX_SIDE_LENGTH = 5.0   # Side legth of the simulation box
    R_CUT = 3.0             # Cutoff radius beyond which interactions are ignored
    PARTICLE_NUMBER = 25    # Total number of particles in the simulation
    SPEED_AMPLITUDE = 11    # Initial speed amplitude of particles
    TIMESTEP = 0.005        # Simulation timestep
    TOTAL_TIME = 7.0        # Total time to propagate
    COOL_FACTOR = 0.997     # Factor multiplying the velocities
    
    ensemble = Ensemble(BOX_SIDE_LENGTH, PARTICLE_NUMBER, R_CUT)
    ensemble.generate_particles()
    ensemble.set_initial_velocities(SPEED_AMPLITUDE)
    ensemble.normalize_momentum()
    ensemble.update_ensemble_acceleration()

See an example `here <https://github.com/aannabe/moldyn/blob/main/examples/cool_temperature.py>`_
