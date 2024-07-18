from pathlib import Path

import pytest

from nanover.app import NanoverImdApplication

EXAMPLES_PATH = Path(__file__).parent
RECORDING_PATH_TRAJ = EXAMPLES_PATH / "nanotube-example-recording.traj"
RECORDING_PATH_STATE = EXAMPLES_PATH / "nanotube-example-recording.state"
ARGON_XML_PATH = EXAMPLES_PATH / "argon_simulation.xml"


@pytest.fixture
def app_server():
    with NanoverImdApplication.basic_server(port=0) as app_server:
        yield app_server
