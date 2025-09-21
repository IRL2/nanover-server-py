import collections
from dataclasses import dataclass
from functools import partial
from typing import Iterable, Callable, TypeVar, Generic

import numpy as np
import numpy.typing as npt

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
    FRAME_INDEX, PARTICLE_VELOCITIES, PARTICLE_FORCES, PARTICLE_FORCES_SYSTEM,
)


P = TypeVar("P")
U = TypeVar("U")


@dataclass(kw_only=True)
class PackingPair(Generic[U, P]):
    pack: Callable[[U], P]
    unpack: Callable[[P], U]


def pack_array(values: Iterable, *, dtype: npt.DTypeLike):
    return np.fromiter(values, dtype=dtype).tobytes()


def unpack_array(buffer: bytes, *, dtype: npt.DTypeLike):
    return list(np.frombuffer(buffer, dtype=dtype))


def make_bytes_packer(dtype: npt.DTypeLike):
    return PackingPair(
        pack=partial(pack_array, dtype=dtype),
        unpack=partial(unpack_array, dtype=dtype),
    )


pack_float32 = make_bytes_packer(np.float32)
pack_uint32 = make_bytes_packer(np.uint32)
pack_uint8 = make_bytes_packer(np.uint8)

pack_identity = PackingPair(pack=lambda value: value, unpack=lambda value: value)  # type: ignore
pack_force_list = PackingPair(pack=list, unpack=list)  # type: ignore
pack_force_int = PackingPair(pack=int, unpack=int)


converters: dict[str, PackingPair] = {
    FRAME_INDEX: pack_force_int,
    PARTICLE_COUNT: pack_force_int,
    CHAIN_COUNT: pack_force_int,
    RESIDUE_COUNT: pack_force_int,
    SIMULATION_COUNTER: pack_force_int,
    PARTICLE_POSITIONS: pack_float32,
    PARTICLE_VELOCITIES: pack_float32,
    PARTICLE_FORCES: pack_float32,
    PARTICLE_FORCES_SYSTEM: pack_float32,
    PARTICLE_ELEMENTS: pack_uint8,
    PARTICLE_RESIDUES: pack_uint32,
    BOND_PAIRS: pack_uint32,
    RESIDUE_CHAINS: pack_uint32,
    BOX_VECTORS: pack_float32,
}


def convert_dict_frame_to_grpc_frame(dict_frame):
    grpc_frame = FrameData()

    for key, value in dict_frame.items():
        if isinstance(value, collections.abc.Sequence) and not isinstance(value, str):
            grpc_frame.arrays[key] = value
        else:
            grpc_frame.values[key] = value

    return grpc_frame


def convert_grpc_frame_to_dict_frame(grpc_frame):
    return unpack_dict_frame(pack_grpc_frame(grpc_frame))


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
    packed = {}

    for key in frame:
        converter = converters.get(key, pack_identity)
        packed[key] = converter.pack(frame[key])

    return packed


def unpack_dict_frame(frame: dict):
    unpacked = {}

    for key in frame:
        converter = converters.get(key, pack_identity)
        unpacked[key] = converter.unpack(frame[key])

    return unpacked
