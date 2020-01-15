# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import pytest
from ase import Atoms

from narupa.imd.imd_client import ImdClient


def co_atoms():
    d = 1.1
    co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)],
               cell=[2, 2, 2],
               pbc=[1, 1, 1])
    return co


@pytest.fixture
def imd_client():
    with ImdClient.insecure_channel(address='localhost', port=54322) as client:
        yield client

