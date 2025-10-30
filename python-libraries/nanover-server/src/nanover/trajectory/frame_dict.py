from typing import Any

from nanover.utilities.packing import (
    PackingPair,
    pack_identity,
    pack_vec3,
    pack_uint8,
    pack_uint32,
    pack_bond,
)
from . import keys

FrameDict = dict[str, Any]


FRAME_PACKERS: dict[str, PackingPair] = {
    keys.PARTICLE_POSITIONS: pack_vec3,
    keys.PARTICLE_VELOCITIES: pack_vec3,
    keys.PARTICLE_FORCES: pack_vec3,
    keys.PARTICLE_FORCES_SYSTEM: pack_vec3,
    keys.PARTICLE_ELEMENTS: pack_uint8,
    keys.PARTICLE_RESIDUES: pack_uint32,
    keys.BOND_PAIRS: pack_bond,
    keys.BOND_ORDERS: pack_uint8,
    keys.RESIDUE_CHAINS: pack_uint32,
    keys.BOX_VECTORS: pack_vec3,
    keys.USER_FORCES_INDEX: pack_uint32,
    keys.USER_FORCES_SPARSE: pack_vec3,
}
"""Mapping of framedata keys to the packer used to pack/unpack between rich data and data ready for MessagePack"""


def merge_frame_dicts(a: FrameDict, b: FrameDict):
    """
    Merge two frame dictionaries, taking all the data of the two frames, preferring data from the second where they
    overlap. Key "frame.index" == 0 is treated as a frame reset whose presence in either frame is propagated to the
    merged frame, and whose presence in the second frame means the data from first frame should be discarded.
    """

    a_reset = a.get(keys.FRAME_INDEX, None) == 0
    b_reset = b.get(keys.FRAME_INDEX, None) == 0

    merged = {}

    if not b_reset:
        merged.update(a)

    merged.update(b)

    if a_reset or b_reset:
        merged[keys.FRAME_INDEX] = 0

    return merged


def _pack_frame_dict(frame: FrameDict) -> FrameDict:
    packed = {}

    for key in frame:
        packer = FRAME_PACKERS.get(key, pack_identity)
        packed[key] = packer.pack(frame[key])

    return packed


def _unpack_frame_dict(frame: FrameDict) -> FrameDict:
    unpacked = {}

    for key in frame:
        try:
            packer = FRAME_PACKERS.get(key, pack_identity)
            unpacked[key] = packer.unpack(frame[key])
        except Exception as e:
            raise RuntimeError(f"unpack {key} {type(frame[key])} failed {e}") from e

    return unpacked


frame_dict_packer = PackingPair(
    pack=_pack_frame_dict,
    unpack=_unpack_frame_dict,
)
"""Pair of functions for packing/unpacking frame data dictionaries of complex types into dictionaries of simpler MessagePack types."""
