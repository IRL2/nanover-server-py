from dataclasses import dataclass
from typing import Iterable, Callable, TypeVar, Generic, Any

import numpy as np
import numpy.typing as npt

import nanover.trajectory.keys as keys

P = TypeVar("P")
U = TypeVar("U")


@dataclass(kw_only=True)
class PackingPair(Generic[U, P]):
    """
    Pair of functions for packing and unpacking rich data into and from simpler types understood by MessagePack.
    """

    pack: Callable[[U], P]
    unpack: Callable[[P], U]


def pack_array(values: Iterable, *, dtype: npt.DTypeLike) -> bytes:
    """
    Serialize an array of numbers into bytes of the given dtype.
    """
    if isinstance(values, np.ndarray):
        return values.astype(dtype=dtype, copy=False).flatten().tobytes()
    else:
        try:
            return np.fromiter(values, dtype=dtype).tobytes()
        except ValueError as e:
            raise ValueError(
                f"Failed to create array of {dtype} from {type(values)} ({values})"
            ) from e


def unpack_array(buffer: bytes, *, dtype: npt.DTypeLike) -> npt.NDArray:
    """
    Deserialize an array of numeric data of the given dtype from bytes.
    """
    return np.frombuffer(buffer, dtype=dtype)


def make_bytes_packer(dtype: npt.DTypeLike, shape: tuple[int, ...] = (-1,)):
    return PackingPair(
        pack=lambda array: pack_array(array, dtype=dtype),
        unpack=lambda data: unpack_array(data, dtype=dtype).reshape(shape),
    )


pack_identity = PackingPair(pack=lambda value: value, unpack=lambda value: value)  # type: ignore

pack_uint32 = make_bytes_packer(np.uint32)
pack_uint8 = make_bytes_packer(np.uint8)

pack_vec3 = make_bytes_packer(np.float32, shape=(-1, 3))
pack_bond = make_bytes_packer(np.uint32, shape=(-1, 2))

converters: dict[str, PackingPair] = {
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


def _pack_dict_frame(frame: dict) -> dict[str, Any]:
    packed = {}

    for key in frame:
        converter = converters.get(key, pack_identity)
        packed[key] = converter.pack(frame[key])

    return packed


def _unpack_dict_frame(frame: dict) -> dict[str, Any]:
    unpacked = {}

    for key in frame:
        try:
            converter = converters.get(key, pack_identity)
            unpacked[key] = converter.unpack(frame[key])
        except Exception as e:
            raise RuntimeError(f"unpack {key} {type(frame[key])} failed {e}") from e

    return unpacked


frame_dict_packer = PackingPair(
    pack=_pack_dict_frame,
    unpack=_unpack_dict_frame,
)
"""Pair of functions for packing/unpacking frame data dictionaries of complex types into dictionaries of simpler MessagePack types."""
