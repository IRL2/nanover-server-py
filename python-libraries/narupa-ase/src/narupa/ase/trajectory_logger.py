# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module containing a trajectory logging class that can be used to output
portable trajectory files from an ASE molecular dynamics simulation.
"""

import datetime
import os
from typing import Optional

import ase.io
from ase import Atoms, units  # type: ignore
from ase.io.formats import filetype, ioformats
from ase.md import Langevin


class UnsupportedFormatError(Exception):
    pass


def validate_ase_can_append_filename(filename: str):
    """
    :raises UnsupportedFormatError: if ase is unable to append the format
        implied by the file extension
    """
    format = _get_format(filename)
    validate_ase_can_append_format(format)


def validate_ase_can_append_format(format: str):
    """
    :raises UnsupportedFormatError: if ase is unable to append the desired
        format.
    """
    try:
        # TODO: when ASE fixes their typo we can do the correct check
        if not ioformats[format].can_write:  # or not ioformats[format].can_append:
            raise UnsupportedFormatError
    except KeyError:
        raise UnsupportedFormatError


def _get_format(filename):
    return filetype(filename, read=False, guess=False)


def _ase_supports_file_descriptor(format: str) -> bool:
    try:
        return ioformats[format].acceptsfd
    except KeyError:
        raise UnsupportedFormatError


class TrajectoryLogger:
    """
    Trajectory logging class for use with ASE simulations.

    Can be attached to an ASE simulation, resulting in frames being written automatically:

    >>> from ase.calculators.emt import EMT
    >>> from ase.lattice.cubic import FaceCenteredCubic
    >>> atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], symbol="Cu", size=(2, 2, 2), pbc=True)
    >>> atoms.calc = EMT()
    >>> dynamics = Langevin(atoms, timestep=0.5, temperature_K=300, friction=1.0)
    >>> with TrajectoryLogger(atoms, 'example.xyz') as logger:
    ...     dynamics.attach(TrajectoryLogger(atoms, 'example.xyz'), interval=10) # attach an XYZ logger.

    :param atoms: ASE :class:`Atoms` from which to write data.
    :param filename: Path to filename to write to.
    :param format: Format to use, as supported by ASE. If not specified, derived from filename.
    :param timestamp: Whether to append a timestamp to the file name. Use to avoid overwriting the same file if
    dynamics is reset.
    :param parallel:  Default is to write on master process only.  Set to `False` to write from all processes.
    :param kwargs: Keyword arguments to be passed to the underlying :fun:`ase.io.write` method.

    Note that even if an ASE format supports appending, if that file requires additional data between frames (e.g.
    headers), then the resulting output may not be valid.

    If ASE supports writing directly a file descriptor for the given format, then that will be used for performance,
    otherwise, the file will be reopened and appended to each frame, which will negatively impact performance.

    """
    format: str

    def __init__(self, atoms: Atoms, filename: str, format: Optional[str] = None, timestamp=True, parallel=True,
                 **kwargs):

        self.frame_index = 0
        self.atoms = atoms
        self.base_path = filename
        self.parallel = parallel
        self._kwargs = kwargs
        self._timestamp = timestamp
        self.current_path = _generate_filename(self.base_path, self.timestamping)
        self._file_descriptor = None

        if format is not None:
            validate_ase_can_append_format(format)
            self.format = format
        else:
            validate_ase_can_append_filename(filename)
            self.format = _get_format(filename)
        self._format_supports_file_descriptor = _ase_supports_file_descriptor(self.format)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @property
    def timestamping(self) -> bool:
        """
        Indicates whether this logger is appending timestamps to the names of any files it produces.

        :return: `True`, if appending timestamps, `False` otherwise.
        """
        return self._timestamp

    def write(self):
        """
        Writes the current state of the atoms to file.
        """

        if self._format_supports_file_descriptor:
            if not self._file_descriptor:
                self._file_descriptor = open(self.current_path, 'w')
            file = self._file_descriptor
        else:
            file = self.current_path

        ase.io.write(file,
                     self.atoms,
                     format=self.format,
                     parallel=self.parallel,
                     append=False,
                     **self._kwargs)
        self.frame_index += 1

    def reset(self):
        """
        Resets the logger, restarting logging with a new file.

        If the logger is set to use timestamps, a new file will be generated with the current time.
        Otherwise, the file will be overwritten.

        ..note

        This method is used in Narupa to produce new trajectory files whenever the
        simulation is reset by the user.

        """
        self.frame_index = 0
        self.close()
        self.current_path = _generate_filename(self.base_path, self.timestamping)

    def __call__(self):
        """
        Method to allow the logger to be called by ASE molecular dynamics logging utility.
        """
        self.write()

    def close(self):
        """
        Closes the current file being written to, if it exists.
        """
        if self._file_descriptor:
            self._file_descriptor.close()
            self._file_descriptor = None


def _get_timestamp():
    now = datetime.datetime.now()
    timestamp = f'{now:%Y_%m_%d__%H_%M_%S}_{int(now.microsecond / 10000):02d}'
    return timestamp


def _generate_filename(path: str, add_timestamp: bool):
    if not add_timestamp:
        return path
    split_path = os.path.splitext(path)
    timestamp = _get_timestamp()
    new_path = f'{split_path[0]}_{timestamp}{split_path[1]}'
    return new_path
