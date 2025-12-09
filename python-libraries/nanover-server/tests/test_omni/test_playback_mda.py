from unittest.mock import patch

import pytest
from MDAnalysis import Universe

from nanover.mdanalysis import UniverseSimulation

from .common import make_loaded_sim, PDB_PATH, DCD_PATH


@pytest.fixture
def example_playback_mda():
    with make_loaded_sim(make_example_playback_mda()) as sim:
        yield sim


def make_example_playback_mda():
    universe = Universe(topology=PDB_PATH, trajectory=DCD_PATH)
    return UniverseSimulation.from_universe(universe)
