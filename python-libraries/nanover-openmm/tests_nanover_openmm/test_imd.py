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


@pytest.fixture
def app_simulation_and_reporter(basic_simulation_with_imd_force):
    simulation, imd_force = basic_simulation_with_imd_force
    with NanoverImdApplication.basic_server(port=0) as app:
        reporter = imd.NanoverImdReporter(
            frame_interval=3,
            force_interval=4,
            include_velocities=False,
            include_forces=False,
            imd_force=imd_force,
            imd_state=app.imd,
            frame_publisher=app.frame_publisher,
        )
        simulation.reporters.append(reporter)
        yield app, simulation, reporter


@pytest.fixture
def app_simulation_and_reporter_with_interactions(app_simulation_and_reporter):
    app, simulation, reporter = app_simulation_and_reporter
    app.imd.insert_interaction(
        "interaction.0",
        ParticleInteraction(
            position=(2.0, 3.0, 1.0),
            particles=(0, 1, 4),
            interaction_type="spring",
        ),
    )
    app.imd.insert_interaction(
        "interaction.1",
        ParticleInteraction(
            position=(10.0, 20.0, 0.0),
            particles=(4, 5),
            interaction_type="spring",
        ),
    )
    return app, simulation, reporter


@pytest.fixture
def app_simulation_and_reporter_with_velocities_and_forces(
    basic_simulation_with_imd_force,
):
    simulation, imd_force = basic_simulation_with_imd_force
    with NanoverImdApplication.basic_server(port=0) as app:
        reporter = imd.NanoverImdReporter(
            frame_interval=3,
            force_interval=4,
            include_velocities=True,
            include_forces=True,
            imd_force=imd_force,
            imd_state=app.imd,
            frame_publisher=app.frame_publisher,
        )
        simulation.reporters.append(reporter)
        yield app, simulation, reporter


@pytest.fixture
def app_simulation_and_reporter_with_constant_force_interactions(
    app_simulation_and_reporter,
):
    app, simulation, reporter = app_simulation_and_reporter
    app.imd.insert_interaction(
        "interaction.0",
        ParticleInteraction(
            position=(0.0, 0.0, 1.0),
            particles=(0, 4),
            interaction_type="constant",
        ),
    )
    return app, simulation, reporter


@pytest.fixture
def single_atom_app_simulation_and_reporter(single_atom_simulation_with_imd_force):
    simulation, imd_force = single_atom_simulation_with_imd_force
    with NanoverImdApplication.basic_server(port=0) as app:
        reporter = imd.NanoverImdReporter(
            frame_interval=5,
            force_interval=5,
            include_velocities=True,
            include_forces=True,
            imd_force=imd_force,
            imd_state=app.imd,
            frame_publisher=app.frame_publisher,
        )
        simulation.reporters.append(reporter)
        yield app, simulation, reporter


@pytest.fixture
def single_atom_app_simulation_and_reporter_with_constant_force_interaction(
    single_atom_app_simulation_and_reporter,
):
    app, simulation, reporter = single_atom_app_simulation_and_reporter
    app.imd.insert_interaction(
        "interaction.0",
        ParticleInteraction(
            position=(0.0, 0.0, 1.0),
            particles=[0],
            interaction_type="constant",
        ),
    )
    return app, simulation, reporter


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


@pytest.mark.parametrize("number_of_forces", (0, 1, 2))
def test_get_imd_forces_from_system(basic_system, number_of_forces):
    """
    :func:`imd.get_imd_forces_from_system` returns all the compatible forces.
    """
    for _ in range(number_of_forces):
        imd.add_imd_force_to_system(basic_system)
    compatible_forces = imd.get_imd_forces_from_system(basic_system)
    assert len(compatible_forces) == number_of_forces


class TestNanoverImdReporter:
    # The name of the method is part of the OpenMM API. It cannot be made to
    # conform PEP8.
    # noinspection PyPep8Naming
    @pytest.mark.parametrize(
        "simulation_step, expected_step", zip(range(7), (3, 2, 1, 1, 2, 1))
    )
    def test_describeNextReport(
        self, app_simulation_and_reporter, simulation_step, expected_step
    ):
        """
        :meth:`NanoverImdReporter.describeNextReport` returns the expected value
        for step.
        """
        # We use a frame interval of 3 and a force interval of 4
        # current step:   0 1 2 3 4 5 6
        # step frame:     3 2 1 3 2 1 3
        # step forces:    4 3 2 1 4 3 2
        # step to return: 3 2 1 2 2 1 2
        expectation = (expected_step, True, False, False, True, False)
        _, simulation, reporter = app_simulation_and_reporter
        reporter.frame_interval = 3
        reporter.force_interval = 4
        simulation.currentStep = simulation_step
        next_report = reporter.describeNextReport(simulation)
        assert next_report == expectation

    def test_report_first_frame_attributes(self, app_simulation_and_reporter):
        """
        When reporting the first frame, the expected attributes are assigned.
        """
        _, simulation, reporter = app_simulation_and_reporter
        simulation.step(1)
        assert reporter.masses == pytest.approx([12, 1, 1, 1, 12, 1, 1, 1])
        assert reporter.n_particles == 8

    def test_report_send_first_frame(self, app_simulation_and_reporter):
        """
        When reporting the first frame, the reporter sends the topology
        and the positions.
        """
        app, simulation, reporter = app_simulation_and_reporter
        request_id = app.frame_publisher._get_new_request_id()
        frame_queues = app.frame_publisher.frame_queues
        with frame_queues.one_queue(request_id, Queue) as publisher_queue:
            simulation.step(1)
            frames = list(publisher_queue.queue)
        assert len(frames) == 1
        assert frames[0].frame_index == 0
        frame = FrameData(frames[0].frame)
        assert_basic_simulation_topology(frame)
        assert np.allclose(frame.particle_positions, BASIC_SIMULATION_POSITIONS)

    def test_report_frame_forces(self, app_simulation_and_reporter_with_interactions):
        """
        Test that user forces are reported within the frame.
        """
        app, simulation, reporter = app_simulation_and_reporter_with_interactions
        request_id = app.frame_publisher._get_new_request_id()
        frame_queues = app.frame_publisher.frame_queues
        with frame_queues.one_queue(request_id, Queue) as publisher_queue:
            simulation.step(10)
            frames = list(publisher_queue.queue)

        # frame 0: topology special case
        # frame 1: no previous frame so no forces?
        # frame 2: forces from frame 1?
        frame = FrameData(frames[2].frame)

        assert set(frame.user_forces_index) == {0, 1, 4, 5}

    def test_sparse_user_forces(self, app_simulation_and_reporter_with_interactions):
        """
        Test that the sparse user forces exist in the frame data when a user applies an iMD force,
        check that the size of the array of forces is equal to the size of the array of corresponding
        indices, and check that none of the elements of the sparse forces array are zero.
        """
        app, simulation, reporter = app_simulation_and_reporter_with_interactions
        request_id = app.frame_publisher._get_new_request_id()
        frame_queues = app.frame_publisher.frame_queues
        with frame_queues.one_queue(request_id, Queue) as publisher_queue:
            simulation.step(10)
            frames = list(publisher_queue.queue)
        frame = FrameData(frames[2].frame)
        assert frame.user_forces_sparse
        assert frame.user_forces_index
        assert len(frame.user_forces_sparse) >= 1
        assert len(frame.user_forces_sparse) == len(frame.user_forces_index)
        assert np.all(frame.user_forces_sparse) != 0.0

    def test_sparse_user_forces_elements(
        self, app_simulation_and_reporter_with_constant_force_interactions
    ):
        """
        Test that the values of the sparse user forces are approximately as expected from the initial
        positions and the position from which the user force is applied, for a constant force acting
        on the C atoms.
        """
        app, simulation, reporter = (
            app_simulation_and_reporter_with_constant_force_interactions
        )
        request_id = app.frame_publisher._get_new_request_id()
        frame_queues = app.frame_publisher.frame_queues
        with frame_queues.one_queue(request_id, Queue) as publisher_queue:
            simulation.step(10)
            frames = list(publisher_queue.queue)
        frame = FrameData(frames[2].frame)
        # For a mass-weighted constant force applied at [0.0, 0.0, 1.0] to the COM of the C atoms
        mass_weighted_user_forces_t0 = [[0.0, 0.0, -6.0], [0.0, 0.0, -6.0]]
        assert set(frame.user_forces_index) == {0, 4}
        for i in range(len(frame.user_forces_index)):
            assert frame.user_forces_sparse[i] == pytest.approx(
                mass_weighted_user_forces_t0[i], abs=2e-3
            )

    @pytest.mark.parametrize("interval", (1, 2, 3, 4))
    def test_send_frame_frequency(self, app_simulation_and_reporter, interval):
        """
        The expected number of frames is sent.
        """
        app, simulation, reporter = app_simulation_and_reporter
        reporter.frame_interval = interval
        request_id = app.frame_publisher._get_new_request_id()
        frame_queues = app.frame_publisher.frame_queues
        n_steps = 13
        with frame_queues.one_queue(request_id, Queue) as publisher_queue:
            simulation.step(n_steps)
            frames = list(publisher_queue.queue)
        assert len(frames) == (n_steps // reporter.frame_interval) + 1

    @staticmethod
    def assert_forces(reporter, affected_indices, unaffected_indices):
        num_particles = reporter.imd_force.getNumParticles()
        parameters = [
            reporter.imd_force.getParticleParameters(i) for i in range(num_particles)
        ]

        forces_affected = np.array(
            [parameters[particle][1] for particle in affected_indices]
        )
        forces_unaffected = np.array(
            [parameters[particle][1] for particle in unaffected_indices]
        )
        assert np.all(forces_affected != 0)
        assert np.all(forces_unaffected == 0)

    def test_apply_interactions(self, app_simulation_and_reporter_with_interactions):
        """
        Interactions are applied and the computed forces are passed to the imd
        force object.
        """
        app, simulation, reporter = app_simulation_and_reporter_with_interactions

        reporter.force_interval = 1
        simulation.step(1)

        self.assert_forces(
            reporter,
            affected_indices=(0, 1, 4, 5),
            unaffected_indices=(2, 3, 6, 7),
        )

    def test_remove_interaction_partial(
        self, app_simulation_and_reporter_with_interactions
    ):
        """
        When an interaction is removed, the corresponding forces are reset.
        """
        app, simulation, reporter = app_simulation_and_reporter_with_interactions

        reporter.force_interval = 1
        simulation.step(1)
        app.imd.remove_interaction("interaction.0")
        simulation.step(1)

        self.assert_forces(
            reporter,
            affected_indices=(4, 5),
            unaffected_indices=(0, 1, 2, 3, 6, 7),
        )

    def test_remove_interaction_complete(
        self, app_simulation_and_reporter_with_interactions
    ):
        """
        When all interactions are removed, all the corresponding forces are
        reset.
        """
        app, simulation, reporter = app_simulation_and_reporter_with_interactions

        reporter.force_interval = 1
        simulation.step(1)
        app.imd.remove_interaction("interaction.0")
        app.imd.remove_interaction("interaction.1")
        simulation.step(1)

        self.assert_forces(
            reporter,
            affected_indices=[],
            unaffected_indices=(0, 1, 2, 3, 4, 5, 6, 7),
        )

    def test_interactions_interval(self, app_simulation_and_reporter_with_interactions):
        """
        Interactions are updated when expected.
        """
        app, simulation, reporter = app_simulation_and_reporter_with_interactions

        reporter.force_interval = 3
        # Interactions are ignored until the first update that will happen at
        # step 3.
        simulation.step(3)
        # The interactions have been picked up. At step 3. The next update will
        # happen at step 6.
        app.imd.remove_interaction("interaction.0")
        simulation.step(1)
        # Since we are not at step 6 yet, the forces must account for the
        # interaction we removed.
        self.assert_forces(
            reporter,
            affected_indices=(0, 1, 4, 5),
            unaffected_indices=(2, 3, 6, 7),
        )

    @pytest.mark.parametrize(
        "simulation_step, expected_step", zip(range(7), (3, 2, 1, 1, 2, 1))
    )
    def test_describeNextReportIncludingVelocitiesAndForces(
        self,
        app_simulation_and_reporter_with_velocities_and_forces,
        simulation_step,
        expected_step,
    ):
        """
        :meth:`NanoverImdReporter.describeNextReport` returns the expected value
        for step.
        """
        # We use a frame interval of 3 and a force interval of 4
        # current step:   0 1 2 3 4 5 6
        # step frame:     3 2 1 3 2 1 3
        # step forces:    4 3 2 1 4 3 2
        # step to return: 3 2 1 2 2 1 2
        expectation = (expected_step, True, True, True, True, False)
        _, simulation, reporter = app_simulation_and_reporter_with_velocities_and_forces
        reporter.frame_interval = 3
        reporter.force_interval = 4
        simulation.currentStep = simulation_step
        next_report = reporter.describeNextReport(simulation)
        assert next_report == expectation

    def test_velocities_and_forces(
        self, app_simulation_and_reporter_with_velocities_and_forces
    ):
        """
        Test the particle velocities and particle forces that can be optionally included
        when running OpenMM simulations. Assert that these arrays exist, have the same
        length as the particle positions array and are non-zero.
        """
        app, simulation, reporter = (
            app_simulation_and_reporter_with_velocities_and_forces
        )
        request_id = app.frame_publisher._get_new_request_id()
        frame_queues = app.frame_publisher.frame_queues
        with frame_queues.one_queue(request_id, Queue) as publisher_queue:
            simulation.step(10)
            frames = list(publisher_queue.queue)
        frame = FrameData(frames[2].frame)
        assert frame.particle_velocities
        assert frame.particle_forces
        assert len(frame.particle_velocities) == len(frame.particle_positions)
        assert len(frame.particle_forces) == len(frame.particle_positions)
        assert np.all(frame.particle_velocities) != 0.0
        assert np.all(frame.particle_forces) != 0.0

    def test_velocities_and_forces_single_atom(
        self, single_atom_app_simulation_and_reporter_with_constant_force_interaction
    ):
        app, simulation, reporter = (
            single_atom_app_simulation_and_reporter_with_constant_force_interaction
        )
        """
        Numerically test the optionally included velocities and forces being passed
        from OpenMM. This test checks that the velocities and forces arrays have the
        same length as the particle positions array, that the forces array is the 
        same as the user forces array (which should be true for the second frame of a 
        single atom system using the Verlet integrator with a constant force), and 
        then numerically checks that the values of the velocities and forces arrays 
        are as expected.
        """
        # The force is a constant force which should cause the particle to accelerate
        # at 1 nm ps^-1. Thus the expected force (along a single axis) for an argon
        # atom with a mass of 40 amu is 40 kJ mol^-1 nm^-1. The force is applied from
        # the frame with index 1, with a simulation step size of 2 fs and a frame and
        # force interval of 5 simulation steps. Therefore, the expected velocity after
        # 5 simulation steps (0.01 ps) is 0.01 nm ps^-1.
        expected_forces = [0.0, 0.0, 40.0]
        expected_velocities = [0.0, 0.0, 0.01]
        request_id = app.frame_publisher._get_new_request_id()
        frame_queues = app.frame_publisher.frame_queues
        with frame_queues.one_queue(request_id, Queue) as publisher_queue:
            simulation.step(10)
            frames = list(publisher_queue.queue)
        frame = FrameData(frames[2].frame)
        assert frame.particle_velocities
        assert frame.particle_forces
        assert frame.user_forces_sparse
        assert len(frame.particle_velocities) == len(frame.particle_positions)
        assert len(frame.particle_forces) == len(frame.particle_positions)
        for i in range(len(frame.particle_forces)):
            assert frame.particle_forces[i] == pytest.approx(
                frame.user_forces_sparse[i], abs=1e-10
            )
            assert frame.particle_forces[i] == pytest.approx(expected_forces, abs=1e-10)
            assert frame.particle_velocities[i] == pytest.approx(
                expected_velocities, abs=1e-7
            )
            print(frame.particle_velocities[i])
