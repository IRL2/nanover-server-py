import numpy as np
import pytest

from nanover.openmm import imd
from nanover.imd import ParticleInteraction
from nanover.omni.openmm import OpenMMSimulation
from nanover.app import NanoverImdApplication, NanoverImdClient

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

    yield app_server, sim


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


def test_work_done_server(single_atom_app_and_simulation_with_constant_force):
    """
    Test that the calculated user work done on a single atom system gives the
    expected numerical result within the code.
    """

    # For a simulation with a frame interval of 5 simulation steps and a simulation
    # step size of 2 fs, using a constant force to accelerate the atom at 1 nm ps-1
    # (for Ar this is a force of 40 kJ mol-1 nm-1), after 100 steps the work done on
    # an Argon atom should be 20.0 kJ mol-1 analytically. Allowing for numerical error,
    # the expected value should be slightly greater than 20.0 kJ mol-1

    app, sim = single_atom_app_and_simulation_with_constant_force

    # Add step to account for zeroth (topology) frame where force is not applied
    for _ in range(101):
        sim.advance_to_next_report()

    # Check that 1.01 ps have passed (and hence that force has been applied for 1 ps)
    assert sim.simulation.context.getTime()._value == pytest.approx(1.01, abs=10e-12)
    assert sim.work_done == pytest.approx(20.0, abs=0.05)


def test_work_done_frame(single_atom_app_and_simulation_with_constant_force):
    """
    Test that the calculated user work done on a single atom system that appears
    in the frame is equal to the user work done as calculated in the OpenMMSimulation.
    """
    app, sim = single_atom_app_and_simulation_with_constant_force

    for _ in range(11):
        sim.advance_to_next_report()

    with NanoverImdClient.connect_to_single_server(
        port=app.port, address="localhost"
    ) as client:
        client.subscribe_to_frames()
        client.wait_until_first_frame()
        assert client.current_frame.user_work_done == sim.work_done
