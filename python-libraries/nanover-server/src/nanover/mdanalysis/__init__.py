"""
Module providing conversion and utility methods for working with NanoVer and MDAnalysis.
"""

from .simulation import UniverseSimulation
from .converter import mdanalysis_to_frame_data, frame_data_to_mdanalysis
from .universe import (
    NanoverParser,
    NanoverReader,
    explosion_mask,
    universe_from_recording,
    universes_from_recording,
)
