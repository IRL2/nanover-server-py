"""
Tests for :mod:`nanover.openmm.imd`.
"""
import numpy as np
import openmm as mm
from openmm import unit

from nanover.imd import ParticleInteraction
from nanover.openmm import imd, OpenMMSimulation
from openmm.unit import nanometer

from nanover.testing.servers import make_app_server
from simulation_utils import (
    basic_simulation,
    basic_system,
    empty_imd_force,
)
from simulation_utils import build_basic_simulation


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


# TODO: remove if we do different velocity reset later
def test_velocity_reset_linear_motion():
    """Tests that velocity reset removes mean linear motion."""
    with make_app_server() as app_server:
        simulation = OpenMMSimulation.from_simulation(build_basic_simulation())
        simulation.reset(app_server)
        simulation.include_velocities = True

        particles = [0]

        app_server.imd.insert_interaction(
            "interaction.test",
            ParticleInteraction(
                position=[1000, 0, 0],
                particles=particles,
                type="constant",
                scale=100000,
                reset_velocities=True,
            )
        )
        simulation.advance_by_one_step()

        prev_velocities = simulation.make_regular_frame().particle_velocities[particles]
        prev_magnitude = np.linalg.norm(np.average(prev_velocities, axis=0))

        app_server.imd.remove_interaction("interaction.test")
        simulation.advance_by_one_step()

        next_velocities = simulation.make_regular_frame().particle_velocities[particles]
        next_magnitude = np.linalg.norm(np.average(next_velocities, axis=0))

        assert prev_magnitude > 0.05
        assert np.isclose(next_magnitude, 0)
