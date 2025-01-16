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
    :param args: Additional arguments passed to the parent class (NanoverImdClient) constructor.
    :param update_callback: An optional callback function executed each time a new frame is
        received (CURRENTLY UNUSED).
    :param kwargs: Additional arguments passed to the parent class (NanoverImdClient) constructor.
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
        # TODO: Add functionality to update callback functions to allow widget customisation


class FrameDataStructure(nglview.Structure):
    """
    Subclass of the nglview.Structure class that converts FrameData
    objects to formatted strings that can be read by NGLView to
    visualise the molecular system.

    :param frame: The FrameData object containing the data from the
        molecular simulation.
    :param ext: The file extension for the structure representation
        that is passed to NGLView, which defaults to PDB.
    :param params: A dictionary of loading parameters that are passed
        to NGLView (see parent class).
    """

    def __init__(self, frame, ext="pdb", params={}):
        super().__init__()
        self.path = ""
        self.ext = ext
        self.params = params
        self._frame = frame

    def get_structure_string(self):
        """
        A function that overrides the get_structure_string function of
        the parent class to convert the frame to a PDB formatted
        string to be read by NGLView.

        :return: A PDB string of the molecular structure defined in the
            frame.
        """
        return frame_data_to_pdb(self._frame)


def frame_data_to_nglwidget(frame, **kwargs):
    """
    Function that takes a FrameData object and outputs an NGLView
    widget for visualisation of the molecular system described by
    the frame.

    :param frame: The FrameData object containing the data from the
        molecular simulation.
    :param kwargs: Additional keyword arguments passed to the
        NGLWidget constructor.
    :return: An NGLView widget to visualise the molecular system
        described by the frame.
    """
    structure = FrameDataStructure(frame)
    return NGLWidget(structure, **kwargs)


def fill_empty_fields(universe: mda.Universe):
    """
    Set the PDB-specific fields with their default values.

    Some topology fields are specific to PDB files and are often missing
    from Universes. This function set these fields to their default values if
    they are not present already.

    :param universe: The MDAnalysis universe describing the current frame.
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
    """
    Converts an MDAnalysis Universe to a PDB string.

    :param universe: The MDAnalysis universe describing the current frame.
    :return: A PDB string of the molecular structure defined in the
        MDAnalysis universe.
    """
    fill_empty_fields(universe)
    with StringIO() as str_io, mda.coordinates.PDB.PDBWriter(str_io) as writer:
        writer.filename = ""  # See https://github.com/MDAnalysis/mdanalysis/issues/2512
        writer.write(universe.atoms)
        pdb = str_io.getvalue()
    return pdb


def frame_data_to_pdb(frame: FrameData) -> str:
    """
    Converts a FrameData object to a PDB string, by first reading the frame as
    an MDAnalysis universe and then converting it to a PDB string.

    :param frame: The FrameData object containing the data from the
        molecular simulation.
    :return: A PDB string of the molecular structure defined in the
        frame.
    """
    universe = frame_data_to_mdanalysis(frame)
    pdb = mda_to_pdb_str(universe)
    return pdb
