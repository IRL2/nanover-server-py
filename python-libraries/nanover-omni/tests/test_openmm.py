import pytest

from nanover.omni.openmm import OpenMMSimulation

from common import app_server, ARGON_XML_PATH
from nanover.openmm.imd import add_imd_force_to_system


@pytest.fixture
def example_openmm(app_server):
    sim = OpenMMSimulation.from_xml_path(ARGON_XML_PATH)
    sim.load()
    sim.reset(app_server)
    yield sim


def test_step_interval(example_openmm):
    """
    Test that advancing by one step increments the dynamics steps by frame_interval.
    """
    for i in range(5):
        assert (
            example_openmm.simulation.currentStep == i * example_openmm.frame_interval
        )
        example_openmm.advance_by_one_step()
