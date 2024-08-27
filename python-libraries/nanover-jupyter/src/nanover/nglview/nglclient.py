"""
Provides an NGLView python client to connect to a NanoVer server in order to visualise
the molecular system from within a Jupyter notebook (or iPython interface).
"""


from io import StringIO

import nglview
import numpy as np

from nanover.app import NanoverImdClient
from nanover.mdanalysis import frame_data_to_mdanalysis
from nanover.trajectory import FrameData
from nglview import NGLWidget

import MDAnalysis as mda


class NGLClient(NanoverImdClient):
    """
    A python client that enables visualisation of the molecular system via
    an NGLView widget.

    Example
    =======

    .. code-block:: python

        from nanover.nglview import NGLClient
        client = NGLClient.autoconnect()
        client.view

    :param dynamic_bonds: A boolean flag that dictates whether bonds should be dynamically updated
        during the simulation.
    :param *args: Additional arguments passed to the parent class (NanoverImdClient) constructor.
    :param update_callback: An optional callback function executed each time a new frame is
        received.
    :param **kwargs: Additional arguments passed to the parent class (NanoverImdClient) constructor.
    """
    def __init__(self, dynamic_bonds=False, *args, update_callback=None, **kwargs):
        self._view = None
        super().__init__(*args, **kwargs)
        self.subscribe_to_frames()
        self.update_callback = update_callback
        self.dynamic_bonds = dynamic_bonds

    @property
    def view(self):
        """
        Returns an NGLView widget to visualise the molecular system.
        """
        if self._view is None or self.dynamic_bonds:
            self._view = frame_data_to_nglwidget(self.latest_frame)
        return self._view

    def _on_frame_received(self, frame_index: int, frame):
        """
        On receiving the latest frame, defines the new coordinates of the atoms
        in the molecular system in Angstrom for visualisation using NGLView.
        """
        super()._on_frame_received(frame_index, frame)
        self.view.set_coordinates(
            {0: np.array(self.latest_frame.particle_positions) * 10}
        )
        if self.update_callback is not None:
            self.update_callback(self.universe)


class FrameDataStructure(nglview.Structure):
    def __init__(self, frame, ext="pdb", params={}):
        super().__init__()
        self.path = ""
        self.ext = ext
        self.params = params
        self._frame = frame

    def get_structure_string(self):
        return frame_data_to_pdb(self._frame)


def frame_data_to_nglwidget(frame, **kwargs):
    structure = FrameDataStructure(frame)
    return NGLWidget(structure, **kwargs)


def fill_empty_fields(universe: mda.Universe):
    """
    Set the PDB-specific fields with their default values.

    Some topology fields are specific to PDB files and are often missing
    from Universes. This function set these fields to their default values if
    they are not present already.
    """
    defaults_per_atom = (
        {
            "altLocs": " ",
            "occupancies": 1.0,
            "tempfactors": 0.0,
        },
        len(universe.atoms),
    )
    defaults_per_residue = (
        {
            "icodes": " ",
        },
        len(universe.residues),
    )
    all_defaults = (defaults_per_atom, defaults_per_residue)
    for source_of_defaults, n_elements in all_defaults:
        for key, default_value in source_of_defaults.items():
            if not hasattr(universe.atoms, key):
                universe.add_TopologyAttr(key, [default_value] * n_elements)


def mda_to_pdb_str(universe: mda.Universe):
    fill_empty_fields(universe)
    with StringIO() as str_io, mda.coordinates.PDB.PDBWriter(str_io) as writer:
        writer.filename = ""  # See https://github.com/MDAnalysis/mdanalysis/issues/2512
        writer.write(universe.atoms)
        pdb = str_io.getvalue()
    return pdb


def frame_data_to_pdb(frame: FrameData) -> str:
    universe = frame_data_to_mdanalysis(frame)
    pdb = mda_to_pdb_str(universe)
    return pdb
