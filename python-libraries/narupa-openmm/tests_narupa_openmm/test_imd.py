# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Tests for :mod:`narupa.openmm.imd`.
"""

from queue import Queue
import numpy as np

import pytest
from simtk import openmm as mm
from simtk.unit import nanometer
from narupa.openmm import imd
from narupa.app import NarupaImdApplication
from narupa.openmm.serializer import deserialize_simulation
from narupa.trajectory import FrameData
from narupa.imd.particle_interaction import ParticleInteraction

from .simulation_utils import (
    basic_system,
    basic_simulation,
    basic_simulation_with_imd_force,
    BASIC_SIMULATION_POSITIONS,
    empty_imd_force,
    assert_basic_simulation_topology,
)


@pytest.fixture
def app_simulation_and_reporter(basic_simulation_with_imd_force):
    simulation, imd_force = basic_simulation_with_imd_force
    with NarupaImdApplication.basic_server(port=0) as app:
        reporter = imd.NarupaImdReporter(
            frame_interval=3,
            force_interval=4,
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
        'interaction.0',
        ParticleInteraction(
            position=(2.0, 3.0, 1.0),
            particles=(0, 1, 4),
            interaction_type='spring',
        )
    )
    app.imd.insert_interaction(
        'interaction.1',
        ParticleInteraction(
            position=(10.0, 20.0, 0.0),
            particles=(4, 5),
            interaction_type='spring',
        )
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


@pytest.mark.parametrize('number_of_forces', (0, 1, 2))
def test_get_imd_forces_from_system(basic_system, number_of_forces):
    """
    :func:`imd.get_imd_forces_from_system` returns all the compatible forces.
    """
    for _ in range(number_of_forces):
        imd.add_imd_force_to_system(basic_system)
    compatible_forces = imd.get_imd_forces_from_system(basic_system)
    assert len(compatible_forces) == number_of_forces


class TestNarupaImdReporter:
    # The name of the method is part of the OpenMM API. It cannot be made to
    # conform PEP8.
    # noinspection PyPep8Naming
    @pytest.mark.parametrize('simulation_step, expected_step',
                             zip(range(7), (3, 2, 1, 1, 2, 1)))
    def test_describeNextReport(
            self, app_simulation_and_reporter, simulation_step, expected_step):
        """
        :meth:`NarupaImdReporter.describeNextReport` returns the expected value
        for step.
        """
        # We use a frame interval of 3 and a force interval of 4
        # current step:   0 1 2 3 4 5 6
        # step frame:     3 2 1 3 2 1 3
        # step forces:    4 3 2 1 4 3 2
        # step to return: 3 2 1 2 2 1 2
        expectation = (expected_step, True, False, False, False, True)
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

    @pytest.mark.parametrize('interval', (1, 2, 3, 4))
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
            reporter.imd_force.getParticleParameters(i)
            for i in range(num_particles)
        ]

        forces_affected = np.array([
            parameters[particle][1]
            for particle in affected_indices
        ])
        forces_unaffected = np.array([
            parameters[particle][1]
            for particle in unaffected_indices
        ])
        assert np.all(forces_affected != 0)
        assert np.all(forces_unaffected == 0)

    def test_apply_interactions(
            self, app_simulation_and_reporter_with_interactions):
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
            self, app_simulation_and_reporter_with_interactions):
        """
        When an interaction is removed, the corresponding forces are reset.
        """
        app, simulation, reporter = app_simulation_and_reporter_with_interactions

        reporter.force_interval = 1
        simulation.step(1)
        app.imd.remove_interaction('interaction.0')
        simulation.step(1)

        self.assert_forces(
            reporter,
            affected_indices=(4, 5),
            unaffected_indices=(0, 1, 2, 3, 6, 7),
        )

    def test_remove_interaction_complete(
            self, app_simulation_and_reporter_with_interactions):
        """
        When all interactions are removed, all the corresponding forces are
        reset.
        """
        app, simulation, reporter = app_simulation_and_reporter_with_interactions

        reporter.force_interval = 1
        simulation.step(1)
        app.imd.remove_interaction('interaction.0')
        app.imd.remove_interaction('interaction.1')
        simulation.step(1)

        self.assert_forces(
            reporter,
            affected_indices=[],
            unaffected_indices=(0, 1, 2, 3, 4, 5, 6, 7),
        )

    def test_interactions_interval(
            self, app_simulation_and_reporter_with_interactions):
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
        app.imd.remove_interaction('interaction.0')
        simulation.step(1)
        # Since we are not at step 6 yet, the forces must account for the
        # interaction we removed.
        self.assert_forces(
            reporter,
            affected_indices=(0, 1, 4, 5),
            unaffected_indices=(2, 3, 6, 7),
        )
