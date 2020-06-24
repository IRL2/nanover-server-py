# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Tests for :mod:`narupa.openmm.imd`.
"""

import pytest
from simtk import openmm as mm
from narupa.openmm import imd

from .simulation_utils import basic_system


@pytest.fixture
def empty_imd_force():
    return imd.create_imd_force()


def test_create_imd_force(empty_imd_force):
    """
    The force created has the expected parameters per particle.
    """
    num_per_particle_parameters = empty_imd_force.getNumPerParticleParameters()
    parameter_names = [
        empty_imd_force.getPerParticleParameterName(i)
        for i in range(num_per_particle_parameters)
    ]
    assert parameter_names == ['fx', 'fy', 'fz']


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
    particle_parameters = [
        force.getParticleParameters(i)
        for i in range(num_particles)
    ]
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
    The force returned by :fun:`imd.add_imd_force_to_system` has the expected
    per particle parameters.
    """
    force = imd.add_imd_force_to_system(basic_system)
    assert_fresh_force_particle_parameters(force, basic_system)


def test_add_imd_force_to_system_force_is_in_system(basic_system):
    """
    When using :fun:`imd.add_imd_force_to_system`, the force is indeed added to
    the system.
    """
    force_added = imd.add_imd_force_to_system(basic_system)
    force_obtained = basic_system.getForce(0)
    # The forces are the same if by modifying one we also modify the other.
    force_added.setParticleParameters(0, 0, (1.0, 2.0, 3.0))
    parameters = force_obtained.getParticleParameters(0)
    assert parameters == [0, (1.0, 2.0, 3.0)]
