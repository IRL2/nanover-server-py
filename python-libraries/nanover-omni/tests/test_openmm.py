import numpy as np
import pytest

from nanover.openmm import imd
from nanover.imd import ParticleInteraction
from nanover.omni.openmm import OpenMMSimulation
from nanover.app import NanoverImdApplication

from common import app_server

from openmm_simulation_utils import (
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


@pytest.fixture
def example_openmm(app_server, single_atom_simulation):
    sim = OpenMMSimulation.from_simulation(single_atom_simulation)
    sim.load()
    sim.reset(app_server)
    yield sim


@pytest.fixture
def single_atom_app_and_simulation_with_constant_force(
    app_server, single_atom_simulation
):
    sim = OpenMMSimulation.from_simulation(single_atom_simulation)
    sim.load()
    sim.reset(app_server)

    # Add a constant force interaction with the force positioned far
    # from the origin (the initial position of the atom)
    interaction = ParticleInteraction(
        interaction_type="constant",
        position=(0.0, 0.0, 1000.0),
        particles=[0],
        scale=1,
    )

    app_server.imd.insert_interaction("interaction.test", interaction)

    return app_server, sim


def test_auto_force(app_server, single_atom_simulation):
    """
    Test that interactions work if the imd force isn't added manually.
    """
    omni_sim = OpenMMSimulation.from_simulation(single_atom_simulation)
    omni_sim.load()
    omni_sim.reset(app_server)

    def get_position():
        positions = omni_sim.simulation.context.getState(
            getPositions=True
        ).getPositions(asNumpy=True)
        return np.asarray(positions[0])

    # add an interaction far to the right
    prev_pos = get_position()
    next_pos = list(prev_pos)
    next_pos[0] += 100

    interaction = ParticleInteraction(
        interaction_type="constant",
        position=next_pos,
        particles=[0],
        scale=10,
    )
    app_server.imd.insert_interaction("interaction.test", interaction)

    # run some simulation steps
    for _ in range(50):
        omni_sim.advance_by_one_step()

    # check the atom moved some way to the right
    curr_pos = get_position()
    assert curr_pos[0] - prev_pos[0] >= 1


def test_step_interval(example_openmm):
    """
    Test that advancing by one step increments the dynamics steps by frame_interval.
    """
    for i in range(5):
        assert (
            example_openmm.simulation.currentStep == i * example_openmm.frame_interval
        )
        example_openmm.advance_by_one_step()


def test_work_done(single_atom_app_and_simulation_with_constant_force):
    """
    Test that the calculated user work done on a single atom system gives the
    expected result.
    """

    app, sim = single_atom_app_and_simulation_with_constant_force

    def get_velocity():
        velocities = sim.simulation.context.getState(
            getVelocities=True
        ).getVelocities(asNumpy=True)
        return np.asarray(velocities[0])

    initial_vel = get_velocity()

    for _ in range(11):
        sim.advance_to_next_report()
        step_count = sim.simulation.currentStep
        print("\n%s" % step_count)
        velocity = get_velocity()
        print(velocity)

    final_vel = get_velocity()
    print(sim.work_done)
    assert (final_vel[2] - initial_vel[2]) > 0.
