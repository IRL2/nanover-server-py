"""
Module providing conversion and utility methods for working with NanoVer and MDAnalysis.
"""

from .converter import (
    frame_data_to_mdanalysis as frame_data_to_mdanalysis,
    mdanalysis_to_frame_data as mdanalysis_to_frame_data,
)
from .simulation import UniverseSimulation as UniverseSimulation
from .universe import (
    universe_from_recording as universe_from_recording,
    universes_from_recording as universes_from_recording,
)
