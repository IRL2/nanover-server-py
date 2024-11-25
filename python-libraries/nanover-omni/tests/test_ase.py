import pytest
import numpy as np
from ase import units, Atoms
import ase.units as ase_units
from ase.calculators.lj import LennardJones
from ase.md import VelocityVerlet

from nanover.imd import ParticleInteraction
from nanover.app import NanoverImdClient
from nanover.omni import OmniRunner
from nanover.omni.ase import ASESimulation
from nanover.ase.converter import ASE_TIME_UNIT_TO_PS, ANG_TO_NM, EV_TO_KJMOL

from common import app_server


@pytest.fixture
def example_ase(app_server, example_dynamics):
    sim = ASESimulation.from_ase_dynamics(example_dynamics)
    sim.load()
    sim.reset(app_server)
    yield sim


@pytest.fixture
def example_ase_app_sim(app_server, example_dynamics):
    sim = ASESimulation.from_ase_dynamics(example_dynamics)
    sim.include_forces = True
    sim.include_velocities = True
    sim.load()
    sim.reset(app_server)
    yield app_server, sim


@pytest.fixture
def example_ase_app_sim_constant_force_interaction(example_ase_app_sim):
    app, sim = example_ase_app_sim

    # Add iMD interaction
    sim.app_server.imd.insert_interaction(
        "interaction.0",
        ParticleInteraction(
            position=(1.0, 2.0, 10.0),
            particles=[0],
            interaction_type="constant",
        ),
    )
    yield app, sim


@pytest.fixture
def example_dynamics():
    atoms = Atoms("Ar", positions=[(0, 0, 0)], masses=[40.0], cell=[2, 2, 2])
    atoms.calc = LennardJones()
    dynamics = VelocityVerlet(atoms, timestep=0.5 * ase_units.fs)
    yield dynamics


@pytest.fixture
def multiple_atom_dynamics():
    # Generate system of multiple Argon atoms
    atoms = Atoms(
        ["Ar" for _ in range(10)],
        positions=[(0, 0, i) for i in range(10)],
        cell=[2, 2, 2],
    )
    atoms.calc = LennardJones()
    dynamics = VelocityVerlet(atoms, timestep=0.5 * ase_units.fs)
    yield dynamics


@pytest.fixture
def multiple_atom_ase_app_sim(app_server, multiple_atom_dynamics):
    sim = ASESimulation.from_ase_dynamics(multiple_atom_dynamics)
    sim.include_forces = True
    sim.include_velocities = True
    sim.load()
    sim.reset(app_server)
    yield app_server, sim


@pytest.fixture
def multiple_atom_ase_app_sim_multiple_interactions(multiple_atom_ase_app_sim):
    app, sim = multiple_atom_ase_app_sim

    # Add iMD interactions
    sim.app_server.imd.insert_interaction(
        "interaction.0",
        ParticleInteraction(
            position=(1.0, 2.0, 10.0),
            particles=[0, 4],
            interaction_type="constant",
        ),
    )
    sim.app_server.imd.insert_interaction(
        "interaction.1",
        ParticleInteraction(
            position=(-1.0, 10.0, 2.0),
            particles=[2, 5, 8],
            interaction_type="spring",
        ),
    )

    yield app, sim


def test_step_interval(example_ase):
    """
    Test that advancing by one step increments the dynamics steps by frame_interval.
    """
    for i in range(5):
        assert (
            example_ase.dynamics.get_number_of_steps() == i * example_ase.frame_interval
        )
        example_ase.advance_by_one_step()


def test_dynamics_interaction(example_ase):
    """
    Test that example dynamics responds to interactions.
    """
    example_ase.app_server.imd.insert_interaction(
        "interaction.0",
        ParticleInteraction(
            position=(0.0, 0.0, 10.0),
            particles=[0],
            interaction_type="constant",
        ),
    )
    for _ in range(31):
        example_ase.advance_by_one_step()

    positions = example_ase.atoms.get_positions()
    (x, y, z) = positions[0]

    # Applying a force of 1 kJ mol-1 nm-1 for
    # t = (0.5 fs * 5 simulation steps per advance * 30 advances) = 75 fs
    # using the velocity verlet algorithm should move the atom by 0.028125 nm
    # using s = u*t + 0.5*a*(t^2). 31 advances performed because there are no
    # iMD forces on the first frame. Allow for numerical error with pytest.approx:
    assert z == pytest.approx(0.028125, abs=1e-8)


def test_simulation_time(example_ase_app_sim):
    """
    Test that the simulation time delivered in the frame matches the elapsed
    simulation time (in ps).
    """
    # Check consistency of unit conversion:
    assert (1.0 / (1e3 * ase_units.fs)) == ASE_TIME_UNIT_TO_PS

    app, sim = example_ase_app_sim

    # Advance simulation by 75 fs (0.075 ps)
    for _ in range(30):
        sim.advance_by_one_step()
    time_elapsed_ps = sim.dynamics.get_time() * ASE_TIME_UNIT_TO_PS
    assert time_elapsed_ps == 75.0 * 1e-3

    # Check that time delivered to client is the same as the time elapsed
    # in the simulation
    with NanoverImdClient.connect_to_single_server(
        port=app.port, address="localhost"
    ) as client:
        client.subscribe_to_frames()
        client.wait_until_first_frame()
        assert client.current_frame.simulation_time == time_elapsed_ps


def test_sparse_user_forces(multiple_atom_ase_app_sim_multiple_interactions):
    """
    Test that the sparse user forces exist in the frame data when a user applies an iMD force,
    check that the size of the array of forces is equal to the size of the array of corresponding
    indices, and check that none of the elements of the sparse forces array are zero.
    """
    app, sim = multiple_atom_ase_app_sim_multiple_interactions

    # Advance simulation by 30 steps with interaction applied
    for _ in range(30):
        sim.advance_by_one_step()

    with NanoverImdClient.connect_to_single_server(
        port=app.port, address="localhost"
    ) as client:
        client.subscribe_to_frames()
        client.wait_until_first_frame()
        frame = client.current_frame
        assert frame.user_forces_sparse
        assert frame.user_forces_index
        assert len(frame.user_forces_index) == 5
        assert len(frame.user_forces_sparse) == len(frame.user_forces_index)
        for f_i in np.array(frame.user_forces_sparse).flatten():
            assert f_i != 0.0


def test_user_energy(example_ase_app_sim_constant_force_interaction):
    """
    Test that the user energy exists in the frame when a user applies a force in iMD,
    and check that that energy equals the expected energy for the constant force
    applied to the single atom (the user energy should be equal to the magnitude
    of the constant force).
    """
    app, sim = example_ase_app_sim_constant_force_interaction
    # Advance simulation by 30 steps with interaction applied
    for _ in range(30):
        sim.advance_by_one_step()

    with NanoverImdClient.connect_to_single_server(
        port=app.port, address="localhost"
    ) as client:
        client.subscribe_to_frames()
        client.wait_until_first_frame()
        frame = client.current_frame
        assert frame.user_energy
        assert frame.user_energy == pytest.approx(
            np.sqrt(np.sum(np.square(frame.user_forces_sparse))), abs=1e-6
        )


def test_particle_forces_system_single_atom(
    example_ase_app_sim_constant_force_interaction,
):
    """
    Test that checks that the system particle forces for a single atom system are zero
    when an iMD force is applied to the atom.
    """
    app, sim = example_ase_app_sim_constant_force_interaction
    # Advance simulation by 30 steps with interaction applied
    for _ in range(30):
        sim.advance_by_one_step()

    with NanoverImdClient.connect_to_single_server(
        port=app.port, address="localhost"
    ) as client:
        client.subscribe_to_frames()
        client.wait_until_first_frame()
        frame = client.current_frame
        assert frame.particle_forces_system
        assert frame.user_forces_sparse
        # Assert elements are approximately zero (allow for small
        # numerical error)
        for f_i in np.squeeze(frame.particle_forces_system):
            assert f_i == pytest.approx(0.0, abs=5e-6)


def test_velocity_unit_conversion(example_ase_app_sim_constant_force_interaction):
    """
    Test that the units of velocity are correctly converted to NanoVer units
    when delivered in the frame data.
    """
    app, sim = example_ase_app_sim_constant_force_interaction
    # Advance simulation by 30 steps with interaction applied
    for _ in range(30):
        sim.advance_by_one_step()

    particle_velocities = np.array(
        sim.atoms.get_velocities() * (ANG_TO_NM / ASE_TIME_UNIT_TO_PS)
    ).flatten()

    with NanoverImdClient.connect_to_single_server(
        port=app.port, address="localhost"
    ) as client:
        client.subscribe_to_frames()
        client.wait_until_first_frame()
        frame = client.current_frame
        assert frame.particle_velocities
        particle_velocities_frame = np.array(frame.particle_velocities).flatten()
        for v_i in range(len(particle_velocities_frame)):
            assert particle_velocities[v_i] == pytest.approx(
                particle_velocities_frame[v_i], rel=1e-7
            )


def test_system_force_unit_conversion(multiple_atom_ase_app_sim):
    """
    Test that the units of force are correctly converted to NanoVer units
    when delivered in the frame data for the system forces.
    """
    app, sim = multiple_atom_ase_app_sim
    # Advance simulation by 30 steps with interaction applied
    for _ in range(30):
        sim.advance_by_one_step()

    # Total force equal to system force in this case, ASE atoms object stores the
    # total forces (sum of the system forces and iMD forces)
    particle_forces = np.array(
        sim.atoms.get_forces() * (EV_TO_KJMOL / ANG_TO_NM)
    ).flatten()

    with NanoverImdClient.connect_to_single_server(
        port=app.port, address="localhost"
    ) as client:
        client.subscribe_to_frames()
        client.wait_until_first_frame()
        frame = client.current_frame
        assert frame.particle_forces_system
        particle_forces_frame = np.array(frame.particle_forces_system).flatten()
        for f_i in range(len(particle_forces_frame)):
            assert particle_forces[f_i] == pytest.approx(
                particle_forces_frame[f_i], rel=1e-7
            )


def test_imd_force_unit_conversion(example_ase_app_sim_constant_force_interaction):
    """
    Test that the units of force are correctly converted to NanoVer units
    when delivered in the frame data for the iMD forces.
    """
    app, sim = example_ase_app_sim_constant_force_interaction
    # Advance simulation by 30 steps with interaction applied
    for _ in range(30):
        sim.advance_by_one_step()

    # Total force equal to iMD force in this case, ASE atoms object stores the
    # total forces (sum of the system forces and iMD forces)
    particle_forces = np.array(
        sim.atoms.get_forces() * (EV_TO_KJMOL / ANG_TO_NM)
    ).flatten()

    with NanoverImdClient.connect_to_single_server(
        port=app.port, address="localhost"
    ) as client:
        client.subscribe_to_frames()
        client.wait_until_first_frame()
        frame = client.current_frame
        assert frame.user_forces_sparse
        particle_forces_frame = np.array(frame.user_forces_sparse).flatten()
        for f_i in range(len(particle_forces_frame)):
            assert particle_forces[f_i] == pytest.approx(
                particle_forces_frame[f_i], rel=1e-7
            )


def test_work_done_server(example_ase_app_sim_constant_force_interaction):
    """
    Test that the calculated user work done on a single atom system gives the
    expected numerical result within the code.
    """

    # For a simulation with a frame interval of 5 simulation steps and a simulation
    # step size of 0.5 fs, using a constant force to accelerate the atom at 1 nm ps-1
    # (for Ar this is a force of 40 kJ mol-1 nm-1), after 401 steps the work done on
    # an Argon atom should be 20.0 kJ mol-1 analytically. Allowing for numerical error,
    # the expected value should be slightly greater than 20.0 kJ mol-1

    app, sim = example_ase_app_sim_constant_force_interaction

    print(sim.atoms.get_masses())

    # Add step to account for zeroth (topology) frame where force is not applied
    for _ in range(401):
        sim.advance_to_next_report()

    # Check that 1.0025 ps have passed (and hence that force has been applied for 1 ps)
    assert sim.dynamics.get_time() * ASE_TIME_UNIT_TO_PS == pytest.approx(1.0025, abs=10e-12)
    assert sim.work_done == pytest.approx(20.0, abs=1e-6)


def test_work_done_frame(example_ase_app_sim_constant_force_interaction):
    """
    Test that the calculated user work done on a single atom system that appears
    in the frame is equal to the user work done as calculated in the ASESimulation.
    """
    app, sim = example_ase_app_sim_constant_force_interaction

    for _ in range(11):
        sim.advance_to_next_report()

    with NanoverImdClient.connect_to_single_server(
        port=app.port, address="localhost"
    ) as client:
        client.subscribe_to_frames()
        client.wait_until_first_frame()
        assert client.current_frame.user_work_done == sim.work_done
