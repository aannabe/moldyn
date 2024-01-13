#!/usr/bin/env python

from context import *

def test_normalize_momentum():
    ensemble = Ensemble(5.0, 25, 3.0)
    ensemble.generate_particles()
    ensemble.set_initial_velocities(10.0)
    ensemble.normalize_momentum()
    assert ensemble.calculate_ensemble_momentum()[0] < 1e-6
    assert ensemble.calculate_ensemble_momentum()[1] < 1e-6
