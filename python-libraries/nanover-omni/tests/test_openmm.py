import numpy as np
import pytest

from nanover.openmm import imd
from nanover.imd import ParticleInteraction
from nanover.omni.openmm import OpenMMSimulation
from nanover.app import NanoverImdApplication

from common import app_server, ARGON_XML_PATH
from nanover.openmm import serializer

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


def test_auto_force(app_server):
    """
    Test that interactions work if the imd force isn't added manually.
    """
    with open(ARGON_XML_PATH) as infile:
        omm_sim = serializer.deserialize_simulation(infile)
    omni_sim = OpenMMSimulation.from_simulation(omm_sim)
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
