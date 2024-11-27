from contextlib import contextmanager
from typing import Set

import numpy as np
import pytest
from joblib.testing import param
from openmm import CustomExternalForce

from nanover.openmm import serializer
from nanover.imd import ParticleInteraction
from nanover.omni.openmm import OpenMMSimulation
from nanover.app import NanoverImdClient

from common import (
    make_app_server,
    make_connected_client_from_runner,
    make_connected_client_from_app_server,
    connect_and_retrieve_first_frame_from_app_server,
)

from openmm_simulation_utils import (
    build_single_atom_simulation,
    build_basic_simulation,
)


@pytest.fixture
def example_openmm():
    with make_app_server() as app_server:
        sim = make_example_openmm()
        sim.load()
        sim.reset(app_server)
        yield sim


def make_example_openmm():
    return OpenMMSimulation.from_simulation(build_single_atom_simulation())


@pytest.fixture
def single_atom_app_and_simulation_with_constant_force():
    with make_single_atom_app_and_simulation_with_constant_force() as (app_server, sim):
        yield app_server, sim


@contextmanager
def make_single_atom_app_and_simulation_with_constant_force():
    with make_app_server() as app_server:
        omm_sim = build_single_atom_simulation()
        sim = OpenMMSimulation.from_simulation(omm_sim)
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
def basic_system_app_and_simulation_with_constant_force():
    with make_app_server() as app_server:
        sim = OpenMMSimulation.from_simulation(build_basic_simulation())
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


def test_auto_force():
    """
    Test that interactions work if the imd force isn't added manually.
    """
    with make_app_server() as app_server:
        omni_sim = OpenMMSimulation.from_simulation(build_single_atom_simulation())
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


def test_force_manager_masses(example_openmm):
    """
    Test that the force manager has the correct masses for the simulated system.
    """
    for _ in range(10):
        example_openmm.advance_by_one_step()
        assert example_openmm.imd_force_manager.masses == pytest.approx([40])


# TODO: could generalise for all three MDs
def test_report_frame_forces(basic_system_app_and_simulation_with_constant_force):
    """
    Test that user forces are reported within the frame.
    """
    app, sim = basic_system_app_and_simulation_with_constant_force
    sim.advance_by_one_step()
    frame = connect_and_retrieve_first_frame_from_app_server(app)

    assert frame.user_forces_index == [0, 6]


def test_sparse_user_forces(basic_system_app_and_simulation_with_constant_force):
    """
    Test that the sparse user forces exist in the frame data when a user applies an iMD force,
    check that the size of the array of forces is equal to the size of the array of corresponding
    indices, and check that none of the elements of the sparse forces array are zero.
    """
    app, sim = basic_system_app_and_simulation_with_constant_force
    sim.advance_by_one_step()
    frame = connect_and_retrieve_first_frame_from_app_server(app)

    assert frame.user_forces_sparse
    assert frame.user_forces_index
    assert len(frame.user_forces_sparse) >= 1
    assert len(frame.user_forces_sparse) == len(frame.user_forces_index)
    assert np.all(frame.user_forces_sparse) != 0.0


# TODO: update for actual system or use system from original test
def test_sparse_user_forces_elements(
    basic_system_app_and_simulation_with_constant_force,
):
    """
    Test that the values of the sparse user forces are approximately as expected from the initial
    positions and the position from which the user force is applied, for a constant force acting
    on the C atoms.
    """
    app, sim = basic_system_app_and_simulation_with_constant_force
    sim.advance_by_one_step()
    frame = connect_and_retrieve_first_frame_from_app_server(app)

    # For a mass-weighted constant force applied at [0.0, 0.0, 1.0] to the COM of the C atoms
    mass_weighted_user_forces_t0 = [[0.0, 0.0, -6.0], [0.0, 0.0, -6.0]]
    assert set(frame.user_forces_index) == {0, 6}
    for i in range(len(frame.user_forces_index)):
        assert frame.user_forces_sparse[i] == pytest.approx(
            mass_weighted_user_forces_t0[i], abs=3e-3
        )


def test_apply_interactions(basic_system_app_and_simulation_with_constant_force):
    """
    Interactions are applied and the computed forces are passed to the imd
    force object.
    """
    app, sim = basic_system_app_and_simulation_with_constant_force
    sim.advance_by_one_step()

    assert_imd_force_affected_particles(
        sim.imd_force_manager.imd_force,
        expected_affected_indices={0, 6},
    )


def assert_imd_force_affected_particles(imd_force: CustomExternalForce, expected_affected_indices: Set[int]):
    """
    Assert that the given imd force is applying force only to the expected particle indices.
    """
    num_particles = imd_force.getNumParticles()

    def particle_is_affected(index: int):
        index, force = imd_force.getParticleParameters(index)
        return any(component != 0 for component in force)

    actual_affected_indices = {index for index in range(num_particles) if particle_is_affected(index)}
    assert actual_affected_indices == expected_affected_indices
