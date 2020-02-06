import os
from contextlib import contextmanager
from io import StringIO
from tempfile import mkstemp

import nglview
import numpy as np

from narupa.app import NarupaImdClient
from narupa.ase.converter import frame_data_to_ase
from narupa.mdanalysis import frame_data_to_mdanalysis, mdanalysis_to_frame_data
from narupa.trajectory import FrameData
from nglview import NGLWidget

import MDAnalysis as mda


class NGLClient(NarupaImdClient):
    def __init__(self, dynamic_bonds=False, *args, update_callback=None,
                 **kwargs):
        self._view = None
        super().__init__(*args, **kwargs)
        self.update_callback = update_callback
        self.dynamic_bonds = dynamic_bonds

    @property
    def view(self):
        if self._view is None or self.dynamic_bonds:
            self._view = show_framedata(self.latest_frame)
        return self._view

    def _on_frame_received(self, frame_index: int, frame):
        super()._on_frame_received(frame_index, frame)
        self.view.set_coordinates(
            {0: np.array(self.latest_frame.particle_positions) * 10}
        )
        if self.update_callback is not None:
            self.update_callback(self.universe)


class FrameDataStructure(nglview.Structure):
    def __init__(self, frame, ext='pdb', params={}):
        super().__init__()
        self.path = ''
        self.ext = ext
        self.params = params
        self._frame = frame

    def get_structure_string(self):
        return frame_data_to_pdb(self._frame)


def show_framedata(frame, **kwargs):
    structure = FrameDataStructure(frame)
    return NGLWidget(structure, **kwargs)


def mda_to_pdb_str(universe: mda.Universe):
    with StringIO() as str_io, mda.coordinates.PDB.PDBWriter(str_io) as writer:
        writer.write(universe.atoms)
        pdb = str_io.getvalue()
    return pdb


def frame_data_to_pdb(frame: FrameData) -> str:
    universe = frame_data_to_mdanalysis(frame)
    pdb = mda_to_pdb_str(universe)
    return pdb
