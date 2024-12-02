import sys
from io import StringIO

from pathlib import Path

import numpy as np
import pytest
from openmm.app import StateDataReporter

from nanover.openmm import imd, serializer
from nanover.imd import ParticleInteraction
from nanover.omni.openmm import OpenMMSimulation
from nanover.app import NanoverImdApplication, NanoverImdClient

from common import app_server
from nanover.trajectory import FrameData

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


@pytest.fixture
def basic_system_app_and_simulation(app_server, basic_simulation):
    sim = OpenMMSimulation.from_simulation(basic_simulation)
    sim.load()
    sim.reset(app_server)

    yield app_server, sim


@pytest.fixture
def basic_system_app_and_simulation_with_constant_force(app_server, basic_simulation):
    sim = OpenMMSimulation.from_simulation(basic_simulation)
    sim.load()
    sim.reset(app_server)

    # Add a constant force interaction with the force positioned far
    # from the origin and a harmonic force sqrt(2) from the origin,
    # acting on different atoms
    interaction_1 = ParticleInteraction(
        interaction_type="constant",
        position=(0.0, 0.0, 1000.0),
        particles=[0],
        scale=1,
    )
    interaction_2 = ParticleInteraction(
        interaction_type="spring",
        position=(1.0, 0.0, 1.0),
        particles=[6],
        scale=1,
    )

    app_server.imd.insert_interaction("interaction.test1", interaction_1)
    app_server.imd.insert_interaction("interaction.test2", interaction_2)

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


def test_save_state_basic_system(basic_system_app_and_simulation_with_constant_force):
    """
    Test that the state of the system can be serialized and deserialized correctly
    by testing that the velocities of the simulation being serialized approximately
    equal those of the simulation loaded after serialization/deserialization.
    """
    app, sim = basic_system_app_and_simulation_with_constant_force

    # Run the simulation for a few steps
    for _ in range(11):
        sim.advance_to_next_report()

    # Get velocities of state before serialization/deserialization
    velocities = sim.simulation.context.getState(getVelocities=True).getVelocities()

    # Serialize/deserialize simulation
    serialized_simulation = serializer.serialize_simulation(
        sim.simulation, save_state=True
    )
    sim_2 = serializer.deserialize_simulation(serialized_simulation)

    # Get velocities of state after serialization/deserialization
    loaded_velocities = sim_2.context.getState(getVelocities=True).getVelocities()

    # Assert that all components of velocities are approximately equal before and after
    # serialization/deserialization procedure.
    for i in range(len(velocities)):
        assert velocities[i].x == pytest.approx(loaded_velocities[i].x, abs=2.0e-7)
        assert velocities[i].y == pytest.approx(loaded_velocities[i].y, abs=2.0e-7)
        assert velocities[i].z == pytest.approx(loaded_velocities[i].z, abs=2.0e-7)


def test_instantaneous_temperature_no_interaction(basic_system_app_and_simulation):
    """
    Test that the instantaneous temperature calculated by NanoVer is equal to the
    instantaneous temperature calculated by the StateDataReporter of OpenMM for a
    simulation without iMD interactions.
    """
    app, sim = basic_system_app_and_simulation

    # Save the output of the StateDataReporter to a variable
    old_stdout = sys.stdout
    sys.stdout = state_data_output = StringIO()

    # Attach StateDataReporter to the simulation
    sim.simulation.reporters.append(
        StateDataReporter(sys.stdout, 1, step=True, temperature=True, append=True)
    )

    # Advance the simulation
    for _ in range(11):
        sim.advance_to_next_report()

    # Retrieve the final temperature from StateDataReporter output
    sys.stdout = old_stdout
    state_data_temperature = float(state_data_output.getvalue().split(",")[-1])

    with NanoverImdClient.connect_to_single_server(
        port=app.port, address="localhost"
    ) as client:
        client.subscribe_to_frames()
        client.wait_until_first_frame()
        assert client.current_frame.system_temperature == state_data_temperature


def test_instantaneous_temperature_imd_interaction(
    basic_system_app_and_simulation_with_constant_force,
):
    """
    Test that the instantaneous temperature calculated by NanoVer is equal to the
    instantaneous temperature calculated by the StateDataReporter of OpenMM to within
    a tolerance (see `issue #324 <https://github.com/IRL2/nanover-server-py/issues/324>`__).
    """
    app, sim = basic_system_app_and_simulation_with_constant_force

    # Save the output of the StateDataReporter to a variable
    old_stdout = sys.stdout
    sys.stdout = state_data_output = StringIO()

    # Attach StateDataReporter to the simulation
    sim.simulation.reporters.append(
        StateDataReporter(sys.stdout, 1, step=True, temperature=True, append=True)
    )

    # Advance the simulation
    for _ in range(101):
        sim.advance_to_next_report()

    # Retrieve the final temperature from StateDataReporter output
    sys.stdout = old_stdout
    state_data_temperature = float(state_data_output.getvalue().split(",")[-1])

    with NanoverImdClient.connect_to_single_server(
        port=app.port, address="localhost"
    ) as client:
        client.subscribe_to_frames()
        client.wait_until_first_frame()
        frame_temperature = client.current_frame.system_temperature
        # Check that the temperature during the iMD interaction is within
        # 1% of the temperature calculated by the StateDataReporter (including
        # the effect of the iMD interaction on the temperature)
        assert frame_temperature == pytest.approx(state_data_temperature, rel=1.0e-2)


def test_reset_gives_equal_frames(app_server):
    """
    Test that resetting the simulation gives frames with equal positions, velocities, and forces etc.
    """

    def fetch_data(frame_data: FrameData):
        return {
            "positions": np.array(frame_data.particle_positions).flatten(),
            "velocities": np.array(frame_data.particle_velocities).flatten(),
            "forces": np.array(frame_data.particle_forces_system).flatten(),
        }

    sim = OpenMMSimulation.from_xml_path(Path(__file__).parent / "hiv1_complex.xml")
    sim.load()

    sim.include_forces = True
    sim.include_velocities = True

    sim.reset(app_server)
    prev_data = fetch_data(sim.make_regular_frame())

    sim.reset(app_server)
    next_data = fetch_data(sim.make_regular_frame())

    assert all(np.all(prev_data[key] - next_data[key]) == 0 for key in prev_data)
