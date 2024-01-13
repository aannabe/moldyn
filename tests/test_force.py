#!/usr/bin/env python

from context import *

def test_force():
    interaction = Interaction(3.0)
    long_force = interaction.LJ_force(4.0)
    short_force = interaction.LJ_force(1.0)
    assert abs(long_force) < 1e-6
    assert abs(short_force - -24.0) < 1e-6
