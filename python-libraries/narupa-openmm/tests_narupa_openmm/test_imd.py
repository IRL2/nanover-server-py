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
    basic_simulation_xml,
    BASIC_SIMULATION_POSITIONS,
    empty_imd_force,
)


@pytest.fixture
def app_simulation_and_reporter(basic_simulation_xml):
    with NarupaImdApplication.basic_server(port=0) as app:
        imd_force = imd.create_imd_force()
        simulation = deserialize_simulation(
            basic_simulation_xml, imd_force=imd_force)

        reporter = imd.NarupaImdReporter(
            frame_interval=3,
            force_interval=4,
            imd_force=imd_force,
            imd_service=app.imd,
            frame_publisher=app.frame_publisher,
        )
        simulation.reporters.append(reporter)
        yield app, simulation, reporter


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
        # The residue names are defined as METH1 and METH2 in the original
        # simulation. However, we built the simulation from an XML file and
        # therefore we used a PDB for the topology. Because the residue name
        # is stored on 3 columns in PDB files, the residue names are truncated.
        assert frame.residue_names == ['MET', 'MET']
        assert frame.residue_chains == [0, 1]
        assert frame.particle_names == ['C1', 'H2', 'H3', 'H4'] * 2
        assert frame.particle_elements == [6, 1, 1, 1] * 2
        assert frame.particle_residues == [0] * 4 + [1] * 4
        assert frame.bond_pairs == [
            [0, 1], [0, 2], [0, 3],  # First residue
            [4, 5], [4, 6], [4, 7],  # Second residue
        ]
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

    def test_apply_forces(self, app_simulation_and_reporter):
        app, simulation, reporter = app_simulation_and_reporter
        app.imd.insert_interaction(
            ParticleInteraction(
                interaction_id='0',
                position=(2.0, 3.0, 1.0),
                particles=(0, 1, 4),
                interaction_type='spring',
            )
        )
        app.imd.insert_interaction(
            ParticleInteraction(
                interaction_id='1',
                position=(10.0, 20.0, 0.0),
                particles=(4, 5),
                interaction_type='spring',
            )
        )

        reporter.force_interval = 1
        simulation.step(1)

        num_particles = reporter.imd_force.getNumParticles()
        parameters = [
            reporter.imd_force.getParticleParameters(i)
            for i in range(num_particles)
        ]

        forces_affected = np.array([
            parameters[particle][1]
            for particle in (0, 1, 4, 5)
        ])
        forces_unaffected = np.array([
            parameters[particle][1]
            for particle in (2, 3, 6, 7)
        ])
        assert np.all(forces_affected != 0)
        assert np.all(forces_unaffected == 0)

