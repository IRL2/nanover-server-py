# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import pytest
from ase import Atoms
from narupa.imd import ImdServer


@pytest.fixture
def imd_server():
    with ImdServer(address='localhost', port=0) as server:
        yield server


def co_atoms():
    d = 1.1
    co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)],
               cell=[2, 2, 2],
               pbc=[1, 1, 1])
    return co

