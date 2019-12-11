import datetime
import os
from typing import Optional

import ase.io
from ase import Atoms


def _get_timestamp():
    now = datetime.datetime.now()
    timestamp = f'{now:%Y:%m:%d_%H_%M_%S}_{int(now.microsecond / 10000):02d}'
    return timestamp


def _generate_filename(path: str, add_timestamp: bool):
    if not add_timestamp:
        return path
    split_path = os.path.splitext(path)
    timestamp = _get_timestamp()
    new_path = f'{split_path[0]}_{timestamp}{split_path[1]}'
    return new_path


class TrajectoryLogger:
    """
    Trajectory logging class for use with ASE simulations.

    :param atoms: ASE :class:`Atoms` from which to write data.
    :param filename: Path to filename to write to.
    :param format: Format to use, as supported by ASE. If not specified, derived from filename.
    :param parallel:  Default is to write on master process only.  Set to `False` to write from all processes.
    :param timestep: Whether to append a timestamp to the file name. Use to avoid overwriting the same file if
    dynamics is reset.
    :param kwargs: Keyword arguments to be passed to the underlying :fun:`ase.io.write` method.
    """

    def __init__(self, atoms: Atoms, filename: str, format: Optional[str] = None, timestamp=True, parallel=True,
                 **kwargs):
        self.frame_index = 0
        self.atoms = atoms
        self._original_filename = filename
        self.format = format
        self.parallel = parallel
        self._kwargs = kwargs
        self._timestamp = timestamp
        self.filename = _generate_filename(self._original_filename, self.timestamping)

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
        ase.io.write(self.filename,
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
        """
        self.frame_index = 0
        self.filename = _generate_filename(self._original_filename, self.timestamping)

    def __call__(self):
        self.write()
