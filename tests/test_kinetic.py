#!/usr/bin/env python

from context import *

def test_kinetic():
    particle = Particle(velocity=[1.0, 1.0])
    kinetic_energy = particle.get_kinetic_energy()
    assert abs(kinetic_energy - 1.0) < 1e-6
