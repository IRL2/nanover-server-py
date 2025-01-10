from typing import Generator, Tuple

import pytest
from nanover.imd import ImdClient
from nanover.imd import ImdServer
from nanover.imd.particle_interaction import ParticleInteraction


@pytest.fixture
def imd_server() -> Generator[ImdServer, None, None]:
    with ImdServer(address="localhost", port=0) as server:
        yield server


@pytest.fixture
def imd_server_client(imd_server) -> Generator[Tuple[ImdServer, ImdClient], None, None]:
    with ImdClient.insecure_channel(
        address="localhost", port=imd_server.port
    ) as client:
        yield imd_server, client


@pytest.fixture
def interaction():
    return ParticleInteraction()
