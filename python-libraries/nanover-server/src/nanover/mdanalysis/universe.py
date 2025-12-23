"""
Facilities to read a NanoVer trajectory recording into an MDAnalysis Universe.

.. code:: python

    import MDAnalysis as mda
    from nanover.mdanalysis import universe_from_recording

    u = universe_from_recording('input.traj')

    # or if there are multiple sessions in the recording:
    universes = universes_from_recording('input.traj')

.. note::
    A NanoVer trajectory recording can have its topology change over time. It
    can even contain trajectories for unrelated simulations. The topology in an
    MDAnalysis Universe is constant. Only the frames corresponding to the first
    topology are read in the singular `universe_from_recording` function.


"""

import warnings
from contextlib import suppress
from itertools import islice
from os import PathLike
from typing import NamedTuple, Type, Callable

from MDAnalysis import Universe
from MDAnalysis.coordinates.base import ProtoReader
from MDAnalysis.coordinates.timestep import Timestep
from MDAnalysis.lib.util import openany
from MDAnalysis.core.topologyattrs import (
    Atomnames,
    Atomtypes,
    Bonds,
    Elements,
    Resids,
    Resnames,
    Segids,
    ChainIDs,
    TopologyAttr,
)
from MDAnalysis.core.topology import Topology
from MDAnalysis.topology.base import TopologyReaderBase
import numpy as np

from nanover.trajectory import FrameData as FrameData, MissingDataError
import nanover.trajectory.keys as keys

from .converter import _to_chemical_symbol, frame_data_to_mdanalysis
from nanover.recording.reading import (
    MessageZipReader,
    RecordingIndexEntry,
    NanoverRecordingReader,
)


class KeyConversion(NamedTuple):
    attribute: Type[TopologyAttr]
    conversion: Callable


def _as_is(value):
    return value


def _trimmed(value):
    return value.strip()


KEY_TO_ATTRIBUTE = {
    keys.PARTICLE_NAMES: KeyConversion(Atomnames, _trimmed),
    keys.RESIDUE_NAMES: KeyConversion(Resnames, _as_is),
    keys.RESIDUE_IDS: KeyConversion(Resids, _as_is),
    keys.CHAIN_NAMES: KeyConversion(Segids, _as_is),
}


FIRST_FRAME_REQUIRED = {
    keys.PARTICLE_POSITIONS,
    keys.PARTICLE_COUNT,
}


def universe_from_recording(path: PathLike[str], *, convert_units=True):
    """
    Read and convert a NanoVer recording into an mdanalysis Universe, ignore all frames after a frame_index
    reset.
    """
    return Universe(
        path,
        format=NanoverReader,
        topology_format=NanoverParser,
        convert_units=convert_units,
    )


def universes_from_recording(path: PathLike[str], *, convert_units=True):
    """
    Decompose a NanoVer trajectory into an mdanalysis Universe for each session of simulation (determined
    by frame_index resets).
    """
    index_entries: list[RecordingIndexEntry] = []
    universes: list[Universe] = []
    first_particle_frame = FrameData()
    first_frame = last_frame = None

    def frame_begins_next_universe(frame: FrameData):
        return frame.frame_dict.get(keys.FRAME_INDEX, None) == 0

    def finalise_prev_universe():
        nonlocal first_particle_frame, first_frame, last_frame

        reader = MessageZipReader.from_path(path)
        reader.index = list(index_entries)

        try:
            universe = frame_data_to_mdanalysis(first_particle_frame)
            universe.trajectory = NanoverReaderBase(
                reader, filename=path, convert_units=convert_units
            )
            universes.append(universe)
        except Exception as e:
            warnings.warn(
                f"Failed to extract universe in frames #{first_frame}-{last_frame}: {e}"
            )

        index_entries.clear()
        first_particle_frame = FrameData()
        first_frame = last_frame = None

    with NanoverRecordingReader.from_path(path) as reader:
        for i, entry in enumerate(reader):
            frame = reader.get_frame_from_entry(entry)

            if frame is None:
                continue

            if first_frame is None:
                first_frame = i
            last_frame = i

            if frame_begins_next_universe(frame) and index_entries:
                finalise_prev_universe()
            # aggregate initial frames until there is position and topology information
            if not is_valid_first_frame(first_particle_frame):
                first_particle_frame.update(frame)
            index_entries.append(entry)

    if is_valid_first_frame(first_particle_frame) and index_entries:
        finalise_prev_universe()

    return universes


class NanoverParser(TopologyReaderBase):
    def parse(self, **kwargs):
        with openany(self.filename, mode="rb") as infile:
            # We assume the full topology is in the first frame with a particle
            # count greater than 0. This will be true only most of the time.
            # TODO: implement a more reliable way to get the full topology
            try:
                reader = MessageZipReader.from_io(infile)
                first_frame = _trim_start_frame_reader(reader)
            except StopIteration:
                raise IOError("The file does not contain any frame.")

            attrs = []
            for frame_key, (attribute, converter) in KEY_TO_ATTRIBUTE.items():
                with suppress(MissingDataError):
                    values = first_frame[frame_key]
                    attrs.append(attribute([converter(value) for value in values]))

            with suppress(MissingDataError):
                converted_elements = _to_chemical_symbol(first_frame.particle_elements)
                attrs.append(Atomtypes(converted_elements))
                attrs.append(Elements(converted_elements))

            # TODO: generate these values if they are not part of the FrameData
            residx = first_frame.particle_residues
            segidx = first_frame.residue_chains
            n_atoms = first_frame.particle_count
            n_residues = first_frame.residue_count
            n_chains = first_frame.chain_count

            with suppress(MissingDataError):
                chain_ids_per_chain = first_frame.chain_names
                chain_ids_per_particle = [
                    chain_ids_per_chain[segidx[residx[atom]]] for atom in range(n_atoms)
                ]
                attrs.append(ChainIDs(chain_ids_per_particle))

            with suppress(MissingDataError):
                order = None
                with suppress(MissingDataError):
                    order = first_frame.bond_orders

                attrs.append(
                    Bonds(
                        first_frame.bond_pairs,
                        guessed=False,
                        order=order,
                    )
                )

            return Topology(
                n_atoms,
                n_residues,
                n_chains,
                attrs=attrs,
                atom_resindex=residx,
                residue_segindex=segidx,
            )


class NanoverReaderBase(ProtoReader):
    units = {
        "time": "ps",
        "length": "nm",
        "velocity": "nm/ps",
        "force": "kJ/(mol*nm)",
    }

    def __init__(
        self,
        reader: MessageZipReader,
        *,
        filename=None,
        convert_units=True,
        **kwargs,
    ):
        super().__init__()

        self._current_frame_index: int | None = None
        self.convert_units = convert_units
        self.filename = filename
        self.reader = reader

        _strip_reader(self.reader)
        first_frame = _trim_start_frame_reader(self.reader)
        remainder = _trim_end_frame_reader(self.reader)
        self.n_atoms = first_frame.particle_count

        if remainder > 0:
            warnings.warn(
                f"The simulation contains changes to the topology after the "
                f"first frame. Only the frames with the initial topology are "
                f"accessible in this Universe. There are {remainder} "
                f"unread frames with a potentially different topology."
            )

        self._read_frame(0)

    def close(self):
        super().close()
        self.reader.close()

    @property
    def n_frames(self):
        return len(self.reader)

    def _frame_to_timestep(self, frame: int, frame_at_index: FrameData):
        ts = Timestep(self.n_atoms)
        ts.frame = frame

        try:
            ts.positions = frame_at_index.particle_positions
        except MissingDataError as e:
            raise Exception(f"No particle positions in trajectory frame {frame}") from e

        with suppress(MissingDataError):
            ts.time = frame_at_index.simulation_time
        with suppress(MissingDataError):
            ts.triclinic_dimensions = frame_at_index.box_vectors
        with suppress(MissingDataError):
            ts.velocities = frame_at_index.particle_velocities
        try:
            ts.forces = frame_at_index.particle_forces
        except MissingDataError:
            with suppress(MissingDataError):
                ts.forces = frame_at_index.particle_forces_system

        ts.data.update(frame_at_index.frame_dict)
        self._add_user_forces_to_ts(frame_at_index, ts)

        if self.convert_units:
            self.convert_pos_from_native(ts._pos)  # in-place !
            if ts.dimensions is not None:
                self.convert_pos_from_native(ts.dimensions[:3])  # in-place!
            if ts.has_velocities:
                # converts nm/ps to A/ps units
                self.convert_velocities_from_native(ts._velocities)
            if ts.has_forces:
                # converts kJ/(mol*nm) to kJ/(mol*A)
                self.convert_forces_from_native(ts._forces)

        self.ts = ts
        return ts

    def _read_next_timestep(self, ts=None):
        # unsupported
        assert ts is None

        if self._current_frame_index is None:
            frame = 0
        else:
            frame = self._current_frame_index + 1
        return self._read_frame(frame)

    def _reopen(self):
        # no current frame, next frame is 0
        self._current_frame_index = None

    def _add_user_forces_to_ts(self, frame: FrameData, ts):
        """
        Read the user forces from the frame if they are available and
        converts them to the units used by MDAnalysis if required.
        """
        try:
            indices = frame.user_forces_index
            sparse = frame.user_forces_sparse
        except MissingDataError:
            return
        forces = np.zeros((self.n_atoms, 3), dtype=np.float32)
        for index, force in zip(indices, sparse):
            forces[index, :] = force
        if self.convert_units:
            self.convert_forces_from_native(forces)
        ts.data["user_forces"] = forces

    def _read_frame(self, frame: int):
        self._current_frame_index = frame

        try:
            entry = self.reader[frame]
        except IndexError as err:
            raise EOFError(err) from None

        message = self.reader.get_message_from_entry(entry)
        frame_at_index = FrameData.unpack_from_dict(message["frame"])
        frame_at_index["elapsed"] = entry.metadata["timestamp"]

        return self._frame_to_timestep(frame, frame_at_index)


def _strip_reader(reader: MessageZipReader) -> None:
    reader.index = [
        entry for entry in reader.index if "frame" in entry.metadata["types"]
    ]


def _trim_start_frame_reader(reader: MessageZipReader) -> FrameData:
    """
    Change reader index to ignore initial frames that don't contain positions
    and return an aggregate of all skipped frames.
    """
    first_frame = FrameData()
    for i, entry in enumerate(reader):
        message = reader.get_message_from_entry(entry)
        if "frame" in message:
            first_frame.update(FrameData.unpack_from_dict(message["frame"]))
            if is_valid_first_frame(first_frame):
                reader.index = reader.index[i:]
                return first_frame
    raise ValueError("No valid first frame.")


def _trim_end_frame_reader(reader: MessageZipReader):
    """
    Change reader index to ignore all frames after and including the next
    frame_index reset.
    """
    for i, entry in enumerate(reader):
        message = reader.get_message_from_entry(entry)
        if "frame" in message:
            if message["frame"].get(keys.FRAME_INDEX, None) == 0 and i > 0:
                remainder = len(reader.index) - i
                reader.index = reader.index[:i]
                return remainder
    return 0


def has_topology(frame: FrameData) -> bool:
    topology_keys = set(list(KEY_TO_ATTRIBUTE.keys()) + [keys.PARTICLE_ELEMENTS])
    return bool(topology_keys.intersection(frame.frame_dict.keys()))


# Copied from the documentation of itertools. See
# <https://docs.python.org/3/library/itertools.html#itertools.batched>. We will
# be able to use itertools.batch when python 3.12+ will be our minimum
# supported version.
def batched(iterable, n):
    # batched('ABCDEFG', 3) --> ABC DEF G
    if n < 1:
        raise ValueError("n must be at least one")
    it = iter(iterable)
    while batch := tuple(islice(it, n)):
        yield batch


def explosion_mask(trajectory, max_displacement):
    """
    A mask to select the frames that are NOT explosions.

    :param trajectory: The trajectory to mask.
    :param max_displacement: The maximum displacement in Angstrom along
        a given axis before the frame is considered as exploding.
    :return: A list with one boolean per frame in the trajectory.
        The boolean is True is the frame is NOT exploding.
    :raises KeyError: if the frames in the trajectory do not have
        a reset counter. This can be the case for nanover trajectories
        recorded from a server that does not keep track of resets in
        the frames, or if the universe has not been built from a nanover
        trajectory recording.

    Here is an example of how to write a trajectory that excludes the
    exploding frames:

    .. code:: python

        import MDAnalysis as mda
        from nanover.mdanalysis import NanoverParser, NanoverReader, explosion_mask

        u = mda.Universe(
            'hello.traj',
            format=NanoverReader,
            topology_format=NanoverParser,
        )
        mask = explosion_mask(u.trajectory, 100)
        u.atoms.write('hello.pdb')
        u.atoms.write('hello.xtc', frames=mask)

    """
    mask = []
    first = trajectory[0]
    previous = first.positions
    prev_reset = first.data["system.reset.counter"]
    for i, ts in enumerate(trajectory):
        reset = ts.data["system.reset.counter"]
        diff = np.abs(ts.positions - previous)
        has_reset = reset != prev_reset
        is_explosion = (diff.max() > 100 and not has_reset) or np.any(
            ~np.isfinite(ts.positions)
        )
        mask.append(not is_explosion)
        previous = ts.positions
        prev_reset = reset
    return mask


def is_valid_first_frame(frame: FrameData):
    return all(key in frame for key in FIRST_FRAME_REQUIRED)


class NanoverReader(NanoverReaderBase):
    def __init__(self, filename, convert_units=True, **kwargs):
        super().__init__(
            MessageZipReader.from_path(filename),
            filename=filename,
            convert_units=convert_units,
            **kwargs,
        )
