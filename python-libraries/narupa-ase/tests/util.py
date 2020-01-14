# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import pytest
from ase import Atoms
from narupa.imd import ImdServer, ImdClient


@pytest.fixture
def imd_server_client():
    with ImdServer(address='localhost', port=0) as server:
        with ImdClient(address='localhost', port=server.port) as client:
            yield server, client


def co_atoms():
    d = 1.1
    co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)],
               cell=[2, 2, 2],
               pbc=[1, 1, 1])
    return co

