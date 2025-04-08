import sys
from io import StringIO
from contextlib import redirect_stdout, contextmanager

from pathlib import Path
from typing import Set

import numpy as np
import pytest
from openmm import CustomExternalForce
from openmm.app import StateDataReporter

from nanover.openmm import serializer
from nanover.imd import ParticleInteraction
from nanover.omni.openmm import OpenMMSimulation

from nanover.trajectory import FrameData

from common import (
    make_app_server,
    connect_and_retrieve_first_frame_from_app_server,
    make_loaded_sim,
    make_loaded_sim_with_interactions,
)

from openmm_simulation_utils import (
    build_single_atom_simulation,
    build_basic_simulation,
)

UNIT_SIMULATION_BOX_VECTORS = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]


@pytest.fixture
def example_openmm():
    with make_loaded_sim(make_example_openmm()) as sim:
        yield sim


def make_example_openmm():
    return OpenMMSimulation.from_simulation(build_basic_simulation())


@pytest.fixture
def single_atom_app_and_simulation_with_constant_force():
    with make_single_atom_app_and_simulation_with_constant_force() as (app_server, sim):
        yield app_server, sim


@contextmanager
def make_single_atom_app_and_simulation_with_constant_force():
    with make_loaded_sim_with_interactions(
        OpenMMSimulation.from_simulation(build_single_atom_simulation()),
        # Add a constant force interaction with the force positioned far
        # from the origin (the initial position of the atom)
        ParticleInteraction(
            interaction_type="constant",
            position=(0.0, 0.0, 1000.0),
            particles=[0],
            scale=1,
        ),
    ) as (app_server, sim):
        yield app_server, sim


@pytest.fixture
def basic_system_app_and_simulation():
    with make_loaded_sim_with_interactions(
        OpenMMSimulation.from_simulation(build_basic_simulation())
    ) as (app_server, sim):
        yield app_server, sim


@pytest.fixture
def basic_system_app_and_simulation_with_constant_force():
    with make_loaded_sim_with_interactions(
        OpenMMSimulation.from_simulation(build_basic_simulation()),
        # Add a constant force interaction with the force positioned far
        # from the origin and a harmonic force sqrt(2) from the origin,
        # acting on different atoms
        ParticleInteraction(
            interaction_type="constant",
            position=(0.0, 0.0, 1000.0),
            particles=[0],
            scale=1,
        ),
        ParticleInteraction(
            interaction_type="spring",
            position=(1.0, 0.0, 1.0),
            particles=[6],
            scale=1,
        ),
    ) as (app_server, sim):
        yield app_server, sim


@pytest.fixture
def basic_system_app_and_simulation_with_constant_force_old():
    with make_loaded_sim_with_interactions(
        OpenMMSimulation.from_simulation(build_basic_simulation()),
        # Add a constant force interaction with the force positioned at
        # 1 nm along the positive x axis
        ParticleInteraction(
            interaction_type="constant",
            position=(0.0, 0.0, 1.0),
            particles=[0, 4],
            scale=1,
        ),
    ) as (app_server, sim):
        yield app_server, sim


@pytest.fixture
def basic_system_app_and_simulation_with_complex_interactions():
    with make_loaded_sim_with_interactions(
        OpenMMSimulation.from_simulation(build_basic_simulation()),
        ParticleInteraction(
            position=(2.0, 3.0, 1.0),
            particles=(0, 1, 4),
            interaction_type="spring",
        ),
        ParticleInteraction(
            position=(10.0, 20.0, 0.0),
            particles=(4, 5),
            interaction_type="spring",
        ),
    ) as (app_server, sim):
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
    expected numerical result within the code, even across resets.
    """

    # For a simulation with a frame interval of 5 simulation steps and a simulation
    # step size of 2 fs, using a constant force to accelerate the atom at 1 nm ps-1
    # (for Ar this is a force of 40 kJ mol-1 nm-1), after 100 steps the work done on
    # an Argon atom should be 20.0 kJ mol-1 analytically. Allowing for numerical error,
    # the expected value should be slightly greater than 20.0 kJ mol-1

    app, sim = single_atom_app_and_simulation_with_constant_force

    for _ in range(3):
        # Add step to account for zeroth (topology) frame where force is not applied
        for _ in range(101):
            sim.advance_to_next_report()

        # Check that 1.01 ps have passed (and hence that force has been applied for 1 ps)
        assert sim.simulation.context.getTime()._value == pytest.approx(
            1.01, abs=10e-12
        )
        assert sim.work_done == pytest.approx(20.0, abs=0.05)

        # Reset simulation and check work done is now zero
        sim.reset(app)
        assert sim.work_done == 0


def test_work_done_frame(basic_system_app_and_simulation_with_complex_interactions):
    """
    Test that the calculated user work done on a system that appears in the frame is equal
    to the user work done as calculated in the OpenMMSimulation.
    """
    app, sim = basic_system_app_and_simulation_with_complex_interactions

    for _ in range(3):
        sim.reset(app)

        for _ in range(11):
            sim.advance_to_next_report()

        frame = connect_and_retrieve_first_frame_from_app_server(app)
        assert frame.user_work_done == sim.work_done


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
        assert velocities[i].x == pytest.approx(loaded_velocities[i].x, abs=2.0e-6)
        assert velocities[i].y == pytest.approx(loaded_velocities[i].y, abs=2.0e-6)
        assert velocities[i].z == pytest.approx(loaded_velocities[i].z, abs=2.0e-6)


def test_instantaneous_temperature_no_interaction(basic_system_app_and_simulation):
    """
    Test that the instantaneous temperature calculated by NanoVer is equal to the
    instantaneous temperature calculated by the StateDataReporter of OpenMM for a
    simulation without iMD interactions.
    """
    app, sim = basic_system_app_and_simulation

    # Save the output of the StateDataReporter to a variable
    with redirect_stdout(StringIO()) as state_data_output:
        # Attach StateDataReporter to the simulation
        sim.simulation.reporters.append(
            StateDataReporter(sys.stdout, 1, step=True, temperature=True, append=True)
        )

        # Advance the simulation
        for _ in range(11):
            sim.advance_to_next_report()

        state_data_temperature = float(state_data_output.getvalue().split(",")[-1])

    frame = connect_and_retrieve_first_frame_from_app_server(app)
    assert frame.system_temperature == pytest.approx(state_data_temperature, abs=1e-12)


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
    with redirect_stdout(StringIO()) as state_data_output:

        # Attach StateDataReporter to the simulation
        sim.simulation.reporters.append(
            StateDataReporter(sys.stdout, 1, step=True, temperature=True, append=True)
        )

        # Advance the simulation
        for _ in range(101):
            sim.advance_to_next_report()

        state_data_temperature = float(state_data_output.getvalue().split(",")[-1])

    # Check that the temperature during the iMD interaction is within
    # 1% of the temperature calculated by the StateDataReporter (including
    # the effect of the iMD interaction on the temperature)
    frame = connect_and_retrieve_first_frame_from_app_server(app)
    assert frame.system_temperature == pytest.approx(state_data_temperature, rel=1.0e-2)


def test_reset_gives_equal_frames():
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
    sim.include_forces = True
    sim.include_velocities = True

    with make_loaded_sim(sim):
        prev_data = fetch_data(sim.make_regular_frame())
        sim.reset(sim.app_server)
        next_data = fetch_data(sim.make_regular_frame())

    assert all(np.all(prev_data[key] - next_data[key]) == 0 for key in prev_data)


def test_force_manager_masses(basic_system_app_and_simulation):
    """
    Test that the force manager has the correct masses for the simulated system.
    """
    _, sim = basic_system_app_and_simulation

    for _ in range(10):
        sim.advance_by_one_step()
        assert sim.imd_force_manager.masses == pytest.approx([12, 1, 1, 1, 12, 1, 1, 1])


# TODO: could generalise for both OMM and ASE
def test_report_frame_forces(basic_system_app_and_simulation_with_complex_interactions):
    """
    Test that user forces are reported within the frame.
    """
    app, sim = basic_system_app_and_simulation_with_complex_interactions
    sim.advance_by_one_step()
    frame = connect_and_retrieve_first_frame_from_app_server(app)

    assert frame.user_forces_index == [0, 1, 4, 5]


# TODO: could generalise for both OMM and ASE
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


def test_apply_interactions(basic_system_app_and_simulation_with_complex_interactions):
    """
    Interactions are applied and the computed forces are passed to the imd
    force object.
    """
    app, sim = basic_system_app_and_simulation_with_complex_interactions
    sim.advance_by_one_step()

    assert_imd_force_affected_particles(
        sim.imd_force_manager.imd_force,
        expected_affected_indices={0, 1, 4, 5},
    )


def test_remove_interaction_partial(
    basic_system_app_and_simulation_with_complex_interactions,
):
    """
    When an interaction is removed, the corresponding forces are reset.
    """
    app, sim = basic_system_app_and_simulation_with_complex_interactions

    sim.advance_by_one_step()
    app.imd.remove_interaction("interaction.0")
    sim.advance_by_one_step()

    assert_imd_force_affected_particles(
        sim.imd_force_manager.imd_force,
        expected_affected_indices={4, 5},
    )


def test_remove_interaction_complete(
    basic_system_app_and_simulation_with_complex_interactions,
):
    """
    When all interactions are removed, all the corresponding forces are reset.
    """
    app, sim = basic_system_app_and_simulation_with_complex_interactions

    sim.advance_by_one_step()
    app.imd.remove_interaction("interaction.0")
    app.imd.remove_interaction("interaction.1")
    sim.advance_by_one_step()

    assert_imd_force_affected_particles(
        sim.imd_force_manager.imd_force,
        expected_affected_indices=set(),
    )


def test_velocities_and_forces(basic_system_app_and_simulation_with_constant_force):
    """
    Test the particle velocities and particle forces that can be optionally included
    when running OpenMM simulations. Assert that these arrays exist, have the same
    length as the particle positions array and are non-zero.
    """
    app, sim = basic_system_app_and_simulation_with_constant_force

    sim.include_forces = True
    sim.include_velocities = True

    sim.advance_by_one_step()
    frame = connect_and_retrieve_first_frame_from_app_server(app)

    assert frame.particle_velocities
    assert frame.particle_forces_system
    assert len(frame.particle_velocities) == len(frame.particle_positions)
    assert len(frame.particle_forces_system) == len(frame.particle_positions)
    assert np.all(frame.particle_velocities) != 0.0
    assert np.all(frame.particle_forces_system) != 0.0


def assert_imd_force_affected_particles(
    imd_force: CustomExternalForce, expected_affected_indices: Set[int]
):
    """
    Assert that the given imd force is applying force only to the expected particle indices.
    """
    num_particles = imd_force.getNumParticles()

    def particle_is_affected(index: int):
        index, force = imd_force.getParticleParameters(index)
        return any(component != 0 for component in force)

    actual_affected_indices = {
        index for index in range(num_particles) if particle_is_affected(index)
    }
    assert actual_affected_indices == expected_affected_indices


# TODO: could generalise for both MDs
def test_sparse_user_forces_elements(
    basic_system_app_and_simulation_with_constant_force_old,
):
    """
    Test that the values of the sparse user forces are approximately as expected from the initial
    positions and the position from which the user force is applied, for a constant force acting
    on the C atoms.
    """
    app, sim = basic_system_app_and_simulation_with_constant_force_old
    sim.advance_by_one_step()
    frame = connect_and_retrieve_first_frame_from_app_server(app)

    # For a mass-weighted constant force applied at [0.0, 0.0, 1.0] to the COM of the C atoms
    mass_weighted_user_forces_t0 = [[0.0, 0.0, -6.0], [0.0, 0.0, -6.0]]
    assert set(frame.user_forces_index) == {0, 4}
    for i in range(len(frame.user_forces_index)):
        assert frame.user_forces_sparse[i] == pytest.approx(
            mass_weighted_user_forces_t0[i], abs=3e-3
        )


def test_velocities_and_forces_single_atom():
    """
    Numerically test the optionally included velocities and forces being passed
    from OpenMM. This test checks that the velocities and forces arrays have the
    same length as the particle positions array, that the forces array is the
    same as the user forces array (which should be true for the second frame of a
    single atom system using the Verlet integrator with a constant force), and
    then numerically checks that the values of the velocities and forces arrays
    are as expected.
    """

    with make_app_server() as app_server:
        sim = OpenMMSimulation.from_simulation(build_single_atom_simulation())

        sim.include_forces = True
        sim.include_velocities = True

        sim.load()
        sim.reset(app_server)

        app_server.imd.insert_interaction(
            "interaction.0",
            ParticleInteraction(
                position=(0.0, 0.0, 1.0),
                particles=[0],
                interaction_type="constant",
            ),
        )

        sim.advance_by_one_step()
        sim.advance_by_one_step()
        frame = connect_and_retrieve_first_frame_from_app_server(app_server)

    # The force is a constant force which should cause the particle to accelerate
    # at 1 nm ps^-1. Thus the expected force (along a single axis) for an argon
    # atom with a mass of 40 amu is 40 kJ mol^-1 nm^-1. The force is applied from
    # the frame with index 1, with a simulation step size of 2 fs and a frame and
    # force interval of 5 simulation steps. Therefore, the expected velocity after
    # 5 simulation steps (0.01 ps) is 0.01 nm ps^-1.
    expected_forces = [0.0, 0.0, 40.0]
    expected_velocities = [0.0, 0.0, 0.01]

    assert frame.particle_velocities
    assert frame.particle_forces_system
    assert frame.user_forces_sparse
    assert len(frame.particle_velocities) == len(frame.particle_positions)
    assert len(frame.particle_forces_system) == len(frame.particle_positions)

    assert len(frame.user_forces_sparse) == 1
    assert frame.user_forces_sparse[0] == pytest.approx(expected_forces)

    for i in range(len(frame.particle_forces_system)):
        assert frame.particle_velocities[i] == pytest.approx(
            expected_velocities, abs=1e-7
        )


def test_pbc_enforcement():
    """
    Test that PBC wrapping of positions defaults to off, doesn't wrap positions when off, and does wrap positions when on.
    """
    # use simulation with PBC changed to 1x1x1 with some atom positions falling outside the PBC
    omm_sim = build_basic_simulation()
    omm_sim.context.setPeriodicBoxVectors(*UNIT_SIMULATION_BOX_VECTORS)
    sim = OpenMMSimulation.from_simulation(omm_sim)
    sim.load()

    with make_app_server() as app_server:
        sim.reset(app_server)
        sim.advance_by_one_step()

    def out_of_bounds(coord):
        return coord < 0 or coord > 1

    def get_sim_position_coords(sim):
        for position in sim.make_regular_frame().particle_positions:
            for coord in position:
                yield coord

    # should default to not PBC wrapping coords
    assert not sim.enforce_pbc

    # without PBC wrapping, some coords should fall outside the box
    sim.enforce_pbc = False
    assert any(out_of_bounds(coord) for coord in get_sim_position_coords(sim))

    # with PBC wrapping, no coords should fall outside the box
    sim.enforce_pbc = True
    assert not any(out_of_bounds(coord) for coord in get_sim_position_coords(sim))
