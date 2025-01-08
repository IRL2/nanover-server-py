"""
Tests for :mod:`nanover.openmm.imd`.
"""

from queue import Queue
import numpy as np

import pytest
import openmm as mm
from openmm.unit import nanometer
from nanover.openmm import imd
from nanover.app import NanoverImdApplication
from nanover.openmm.serializer import deserialize_simulation
from nanover.trajectory import FrameData
from nanover.imd.particle_interaction import ParticleInteraction

from .simulation_utils import (
    basic_system,
    basic_simulation,
    basic_simulation_with_imd_force,
    BASIC_SIMULATION_POSITIONS,
    empty_imd_force,
    assert_basic_simulation_topology,
    single_atom_system,
    single_atom_simulation,
    single_atom_simulation_with_imd_force,
    ARGON_SIMULATION_POSITION,
    assert_single_atom_simulation_topology,
)


def test_create_imd_force(empty_imd_force):
    """
    The force created has the expected parameters per particle.
    """
    num_per_particle_parameters = empty_imd_force.getNumPerParticleParameters()
    parameter_names = [
        empty_imd_force.getPerParticleParameterName(i)
        for i in range(num_per_particle_parameters)
    ]
    assert parameter_names == ["fx", "fy", "fz"]


def assert_fresh_force_particle_parameters(
    force: mm.CustomExternalForce,
    system: mm.System,
):
    """
    Assert that a freshly populated imd force has the expected per-particle
    parameters for a given system.
    """
    # The first int is a reference to the particle the force applies to,
    # the following tuple is the parameters in x, y, and z.
    num_particles = system.getNumParticles()
    expectations = [[i, (0.0, 0.0, 0.0)] for i in range(num_particles)]
    particle_parameters = [force.getParticleParameters(i) for i in range(num_particles)]
    assert particle_parameters == expectations


def test_populate_imd_force(empty_imd_force, basic_system):
    """
    When populating the imd force, there is the right number of particles,
    the parameters are set to 0, and they refer to the expected particles.
    """
    force = empty_imd_force
    imd.populate_imd_force(force, basic_system)
    assert_fresh_force_particle_parameters(force, basic_system)


def test_add_imd_force_to_system_parameters(basic_system):
    """
    The force returned by :func:`imd.add_imd_force_to_system` has the expected
    per particle parameters.
    """
    force = imd.add_imd_force_to_system(basic_system)
    assert_fresh_force_particle_parameters(force, basic_system)


def test_add_imd_force_to_system_force_is_in_system(basic_system):
    """
    When using :func:`imd.add_imd_force_to_system`, the force is indeed added to
    the system.
    """
    force_added = imd.add_imd_force_to_system(basic_system)
    force_obtained = basic_system.getForce(0)
    # The forces are the same if by modifying one we also modify the other.
    force_added.setParticleParameters(0, 0, (1.0, 2.0, 3.0))
    parameters = force_obtained.getParticleParameters(0)
    assert parameters == [0, (1.0, 2.0, 3.0)]
