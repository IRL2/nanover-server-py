from dataclasses import dataclass
from functools import partial
from typing import Iterable, Callable, TypeVar, Generic

import numpy as np

from nanover.trajectory.frame_data import (
    PARTICLE_COUNT,
    CHAIN_COUNT,
    RESIDUE_COUNT,
    SIMULATION_COUNTER,
    PARTICLE_POSITIONS,
    PARTICLE_ELEMENTS,
    PARTICLE_RESIDUES,
    BOND_PAIRS,
    RESIDUE_CHAINS,
    BOX_VECTORS,
    FrameData,
)


P = TypeVar("P")
U = TypeVar("U")


@dataclass(kw_only=True)
class PackingPair(Generic[U, P]):
    pack: Callable[[U], P]
    unpack: Callable[[P], U]


def pack_array(values: Iterable, *, dtype: str):
    return np.fromiter(values, dtype=dtype).tobytes()


def unpack_array(buffer: bytes, *, dtype: str):
    return np.frombuffer(buffer, dtype=dtype)


def make_bytes_packer(dtype: str):
    return PackingPair(
        pack=partial(pack_array, dtype=dtype),
        unpack=partial(unpack_array, dtype=dtype),
    )


pack_float32 = make_bytes_packer(np.float32)
pack_uint32 = make_bytes_packer(np.uint32)
pack_uint8 = make_bytes_packer(np.uint8)

pack_identity = PackingPair(pack=lambda value: value, unpack=lambda value: value)
pack_force_list = PackingPair(pack=list, unpack=list)
pack_force_int = PackingPair(pack=int, unpack=int)


converters: dict[str, PackingPair] = {
    PARTICLE_COUNT: pack_force_int,
    CHAIN_COUNT: pack_force_int,
    RESIDUE_COUNT: pack_force_int,
    SIMULATION_COUNTER: pack_force_int,
    PARTICLE_POSITIONS: pack_float32,
    PARTICLE_ELEMENTS: pack_uint8,
    PARTICLE_RESIDUES: pack_uint32,
    BOND_PAIRS: pack_uint32,
    RESIDUE_CHAINS: pack_uint32,
    BOX_VECTORS: pack_float32,
}


def pack_grpc_frame(frame: FrameData):
    data = {}

    for key in frame.value_keys:
        converter = converters.get(key, pack_identity)
        data[key] = converter.pack(frame.values[key])

    for key in frame.array_keys:
        converter = converters.get(key, pack_force_list)
        data[key] = converter.pack(frame.arrays[key])

    return data


def pack_dict_frame(frame: dict):
    for key in frame:
        converter = converters.get(key, pack_identity)
        frame[key] = converter.pack(frame[key])

    return frame


def unpack_dict_frame(frame: dict):
    for key in frame:
        converter = converters.get(key, pack_identity)
        frame[key] = converter.unpack(frame[key])

    return frame
