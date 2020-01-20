import os
from contextlib import contextmanager
from tempfile import mkstemp

import nglview
import numpy as np

from narupa.app import NarupaImdClient
from narupa.ase.converter import frame_data_to_ase
from nglview import NGLWidget


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
            atoms = frame_data_to_ase(
                self.first_frame,
                topology=True,
                positions=False,
            )
            atoms.set_positions(
                np.array(self.latest_frame.particle_positions) * 10
            )
            self._view = show_ase(atoms)
        return self._view

    def _on_frame_received(self, frame_index: int, frame):
        super()._on_frame_received(frame_index, frame)
        self.view.set_coordinates(
            {0: np.array(self.latest_frame.particle_positions) * 10}
        )
        if self.update_callback is not None:
            self.update_callback(self.universe)


@contextmanager
def mkstemp_wrapper(*args, **kwargs):
    # NamedTemporaryFile cannot be used here because it makes it impossible
    # on Windows to the file for writing. mkstemp is a bit less restrictive
    # in this regard.
    fd, fname = mkstemp(*args, **kwargs)
    yield fname
    # On windows, the file must be closed before it can be removed.
    os.close(fd)
    os.remove(fname)


def _get_structure_string(write_method, suffix='.pdb'):
    with mkstemp_wrapper(suffix=suffix) as fname:
        write_method(fname)
        with open(fname) as fh:
            return fh.read()


class ASEStructure(nglview.Structure):
    def __init__(self, ase_atoms, ext='pdb', params={}):
        super().__init__()
        self.path = ''
        self.ext = ext
        self.params = params
        self._ase_atoms = ase_atoms

    def get_structure_string(self):
        return _get_structure_string(self._ase_atoms.write)


def show_ase(ase_atoms, **kwargs):
    """
    Examples
    --------
    >>> import nglview as nv
    >>> from ase import Atom, Atoms
    >>> dimer = Atoms([Atom('X', (0, 0, 0)),
    ...                Atom('X', (0, 0, 1))])
    >>> dimer.set_positions([(1, 2, 3), (4, 5, 6.2)])
    >>> w = nv.show_ase(dimer)
    >>> w # doctest: +SKIP
    """
    structure = ASEStructure(ase_atoms)
    return NGLWidget(structure, **kwargs)


nglview.adaptor._get_structure_string = _get_structure_string