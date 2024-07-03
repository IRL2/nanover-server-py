import pytest

from nanover.omni.omni import Simulation
from test_openmm import example_openmm
from test_playback import example_playback
from common import app_server

SIMULATION_FIXTURES = (
    "example_openmm",
    "example_playback",
)


@pytest.mark.parametrize("sim_fixture", SIMULATION_FIXTURES)
def test_reset(sim_fixture, app_server, request):
    sim: Simulation = request.getfixturevalue(sim_fixture)
    sim.reset(app_server)


@pytest.mark.parametrize("sim_fixture", SIMULATION_FIXTURES)
def test_load(sim_fixture, request):
    sim: Simulation = request.getfixturevalue(sim_fixture)
    sim.load()
