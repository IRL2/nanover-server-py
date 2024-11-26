from contextlib import contextmanager
from pathlib import Path
from nanover.app import NanoverImdApplication, NanoverImdClient
from nanover.omni import OmniRunner

EXAMPLES_PATH = Path(__file__).parent
RECORDING_PATH_TRAJ = EXAMPLES_PATH / "nanotube-example-recording.traj"
RECORDING_PATH_STATE = EXAMPLES_PATH / "nanotube-example-recording.state"
ARGON_XML_PATH = EXAMPLES_PATH / "argon_simulation.xml"


@contextmanager
def make_loaded_sim(sim):
    with make_app_server() as app_server:
        sim.load()
        sim.reset(app_server)
        yield sim


@contextmanager
def make_runner(simulations):
    with OmniRunner.with_basic_server(*simulations, port=0) as runner:
        yield runner


@contextmanager
def make_connected_client_from_runner(runner):
    with NanoverImdClient.connect_to_single_server(
        port=runner.app_server.port
    ) as client:
        yield client


@contextmanager
def make_app_server():
    with NanoverImdApplication.basic_server(port=0) as app_server:
        yield app_server



