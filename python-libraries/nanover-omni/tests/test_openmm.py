import pytest

from nanover.omni.openmm import OpenMMSimulation

from common import app_server, ARGON_XML_PATH


@pytest.fixture
def example_openmm(app_server):
    sim = OpenMMSimulation.from_xml_path(ARGON_XML_PATH)
    sim.load()
    sim.reset(app_server)
    yield sim
