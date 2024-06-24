"""
Facilities to read a NanoVer trajectory recording into an MDAnalysis Universe.

.. code:: python

    import MDAnalysis as mda
    from nanover.mdanalysis import NanoverReader, NanoverParser

    u = mda.Universe(
        'input.traj',
        format=NanoverReader,
        topology_format=NanoverParser,
    )

.. note::
    A NanoVer trajectory recording can have its topology change over time. It
    can even contain trajectories for unrelated simulations. The topology in an
    MDAnalysis Universe is constant. Only the frames corresponding to the first
    topology are read in a Universe.

"""

import warnings
from itertools import takewhile, chain, islice
from typing import NamedTuple, Type, Callable

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

from nanover.trajectory import FrameData
from nanover.trajectory.frame_data import (
    PARTICLE_COUNT,
    RESIDUE_COUNT,
    CHAIN_COUNT,
    PARTICLE_ELEMENTS,
    PARTICLE_NAMES,
    PARTICLE_RESIDUES,
    RESIDUE_NAMES,
    RESIDUE_CHAINS,
    RESIDUE_IDS,
    CHAIN_NAMES,
    MissingDataError,
)

from .recordings import (
    Unpacker,
    iter_trajectory_recording,
    iter_trajectory_with_elapsed_integrated,
    advance_to_first_particle_frame,
    advance_to_first_coordinate_frame,
)
from .converter import _to_chemical_symbol


class KeyConversion(NamedTuple):
    attribute: Type[TopologyAttr]
    conversion: Callable


def _as_is(value):
    return value


def _trimmed(value):
    return value.strip()


KEY_TO_ATTRIBUTE = {
    PARTICLE_NAMES: KeyConversion(Atomnames, _trimmed),
    RESIDUE_NAMES: KeyConversion(Resnames, _as_is),
    RESIDUE_IDS: KeyConversion(Resids, _as_is),
    CHAIN_NAMES: KeyConversion(Segids, _as_is),
}


class NanoverParser(TopologyReaderBase):
    def parse(self, **kwargs):
        with openany(self.filename, mode="rb") as infile:
            data = infile.read()
        unpacker = Unpacker(data)
        # We assume the full topology is in the first frame with a particle
        # count greater than 0. This will be true only most of the time.
        # TODO: implement a more reliable way to get the full topology
        try:
            _, _, frame = next(
                advance_to_first_particle_frame(iter_trajectory_recording(unpacker))
            )
        except StopIteration:
            raise IOError("The file does not contain any frame.")

        attrs = []
        for frame_key, (attribute, converter) in KEY_TO_ATTRIBUTE.items():
            try:
                values = frame.arrays[frame_key]
            except MissingDataError:
                pass
            else:
                attrs.append(attribute([converter(value) for value in values]))

        try:
            elements = frame.arrays[PARTICLE_ELEMENTS]
        except MissingDataError:
            pass
        else:
            converted_elements = _to_chemical_symbol(elements)
            attrs.append(Atomtypes(converted_elements))
            attrs.append(Elements(converted_elements))

        # TODO: generate these values if they are not part of the FrameData
        residx = frame.arrays[PARTICLE_RESIDUES]
        segidx = frame.arrays[RESIDUE_CHAINS]
        n_atoms = int(frame.values[PARTICLE_COUNT])
        n_residues = int(frame.values[RESIDUE_COUNT])
        n_chains = int(frame.values[CHAIN_COUNT])

        try:
            chain_ids_per_chain = frame.arrays[CHAIN_NAMES]
        except MissingDataError:
            pass
        else:
            chain_ids_per_particle = [
                chain_ids_per_chain[segidx[residx[atom]]] for atom in range(n_atoms)
            ]
            attrs.append(ChainIDs(chain_ids_per_particle))

        try:
            attrs.append(
                Bonds(frame.bond_pairs, guessed=False, order=frame.bond_orders)
            )
        except MissingDataError:
            pass

        return Topology(
            n_atoms,
            n_residues,
            n_chains,
            attrs=attrs,
            atom_resindex=residx,
            residue_segindex=segidx,
        )


class NanoverReader(ProtoReader):
    units = {"time": "ps", "length": "nm", "velocity": "nm/ps"}

    def __init__(self, filename, convert_units=True, **kwargs):
        super().__init__()
        self.filename = filename
        self.convert_units = convert_units
        with openany(filename, mode="rb") as infile:
            data = infile.read()
        unpacker = Unpacker(data)
        recording = advance_to_first_coordinate_frame(
            iter_trajectory_with_elapsed_integrated(iter_trajectory_recording(unpacker))
        )
        try:
            _, _, first_frame = next(recording)
        except StopIteration:
            raise IOError("Empty trajectory.")
        self.n_atoms = first_frame.particle_count

        non_topology_frames = takewhile(
            lambda frame: not has_topology(frame),
            map(lambda record: record[2], recording),
        )
        self._frames = list(chain([first_frame], non_topology_frames))
        self.n_frames = len(self._frames)
        self._read_frame(0)
        reminder = list(recording)
        if reminder:
            warnings.warn(
                f"The simulation contains changes to the topology after the "
                f"first frame. Only the frames with the initial topology are "
                f"accessible in this Universe. There are {len(reminder)} "
                f"unread frames with a different topology."
            )

    def _read_frame(self, frame):
        self._current_frame_index = frame
        try:
            frame_at_index = self._frames[frame]
        except IndexError as err:
            raise EOFError(err) from None

        ts = Timestep(self.n_atoms)
        ts.frame = frame
        ts.positions = frame_at_index.particle_positions
        try:
            ts.time = frame_at_index.simulation_time
        except MissingDataError:
            pass
        try:
            ts.triclinic_dimensions = frame_at_index.box_vectors
        except MissingDataError:
            pass
        try:
            ts.velocities = frame_at_index.particle_velocities
        except MissingDataError:
            pass
        try:
            ts.forces = frame_at_index.particle_forces
        except MissingDataError:
            pass

        ts.data.update(frame_at_index.values)
        self._add_user_forces_to_ts(frame_at_index, ts)

        if self.convert_units:
            self.convert_pos_from_native(ts._pos)  # in-place !
            if ts.dimensions is not None:
                self.convert_pos_from_native(ts.dimensions[:3])  # in-place!
            if ts.has_velocities:
                # converts nm/ps to A/ps units
                self.convert_velocities_from_native(ts._velocities)

        self.ts = ts
        return ts

    def _read_next_timestep(self):
        if self._current_frame_index is None:
            frame = 0
        else:
            frame = self._current_frame_index + 1
        return self._read_frame(frame)

    def _reopen(self):
        self._current_frame_index = None

    def _add_user_forces_to_ts(self, frame, ts):
        """
        Read the user forces from the frame if they are available.
        """
        try:
            indices = frame.arrays["forces.user.index"]
            sparse = frame.arrays["forces.user.sparse"]
        except KeyError:
            return
        forces = np.zeros((self.n_atoms, 3), dtype=np.float32)
        for index, force in zip(indices, batched(sparse, n=3)):
            forces[index, :] = force
        ts.data["user_forces"] = forces


def has_topology(frame: FrameData) -> bool:
    topology_keys = set(list(KEY_TO_ATTRIBUTE.keys()) + [PARTICLE_ELEMENTS])
    return bool(topology_keys.intersection(frame.array_keys))


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
