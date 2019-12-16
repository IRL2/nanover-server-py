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
from ase import Atoms


class TrajectoryLogger:
    """
    Trajectory logging class for use with ASE simulations.

    :param atoms: ASE :class:`Atoms` from which to write data.
    :param filename: Path to filename to write to.
    :param format: Format to use, as supported by ASE. If not specified, derived from filename.
    :param timestamp: Whether to append a timestamp to the file name. Use to avoid overwriting the same file if
    dynamics is reset.
    :param parallel:  Default is to write on master process only.  Set to `False` to write from all processes.

    :param kwargs: Keyword arguments to be passed to the underlying :fun:`ase.io.write` method.
    """

    def __init__(self, atoms: Atoms, filename: str, format: Optional[str] = None, timestamp=True, parallel=True,
                 **kwargs):
        self.frame_index = 0
        self.atoms = atoms
        self.base_path = filename
        self.format = format
        self.parallel = parallel
        self._kwargs = kwargs
        self._timestamp = timestamp
        self.current_path = _generate_filename(self.base_path, self.timestamping)

    @property
    def timestamping(self) -> bool:
        """
        Indicates whether this logger is appending timestamps to the names of any files it produces.

        :return: `True`, if appending timestamps, `False` otherwise.
        """
        return self._timestamp

    def write(self):
        """
        Writes the current state of the atoms, overwriting if this is the first time the method has been
        called, appending otherwise.
        """
        should_append = self.frame_index != 0
        ase.io.write(self.current_path,
                     self.atoms,
                     format=self.format,
                     parallel=self.parallel,
                     append=should_append,
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
        self.current_path = _generate_filename(self.base_path, self.timestamping)

    def __call__(self):
        """
        Method to allow the logger to be called by ASE molecular dynamics logging utility.
        """
        self.write()


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
