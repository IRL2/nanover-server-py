import collections
from dataclasses import dataclass
from typing import Iterable, Callable, TypeVar, Generic, Any

import numpy as np
import numpy.typing as npt
from nanover.testing.utilities import simplify_numpy

from nanover.protocol.trajectory import GetFrameResponse
from .frame_data2 import FrameData as FrameData2
from .frame_data import (
    FrameData,
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
    FRAME_INDEX,
    PARTICLE_VELOCITIES,
    PARTICLE_FORCES,
    PARTICLE_FORCES_SYSTEM,
    BOND_ORDERS,
    USER_FORCES_INDEX,
    USER_FORCES_SPARSE,
)
from nanover.utilities.change_buffers import DictionaryChange

P = TypeVar("P")
U = TypeVar("U")


@dataclass(kw_only=True)
class PackingPair(Generic[U, P]):
    pack: Callable[[U], P]
    unpack: Callable[[P], U]


def pack_array(values: Iterable, *, dtype: npt.DTypeLike) -> bytes:
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
    return np.frombuffer(buffer, dtype=dtype)


def make_bytes_packer(dtype: npt.DTypeLike, shape: tuple[int, ...] = (-1,)):
    return PackingPair(
        pack=lambda array: pack_array(array, dtype=dtype),
        unpack=lambda data: unpack_array(data, dtype=dtype).reshape(shape),
    )


pack_uint32 = make_bytes_packer(np.uint32)
pack_uint8 = make_bytes_packer(np.uint8)

pack_vec3 = make_bytes_packer(np.float32, shape=(-1, 3))
pack_bond = make_bytes_packer(np.uint32, shape=(-1, 2))

pack_identity = PackingPair(pack=lambda value: value, unpack=lambda value: value)  # type: ignore
pack_force_list = PackingPair(pack=list, unpack=list)  # type: ignore
pack_force_int = PackingPair(pack=int, unpack=int)


converters: dict[str, PackingPair] = {
    FRAME_INDEX: pack_force_int,
    PARTICLE_COUNT: pack_force_int,
    CHAIN_COUNT: pack_force_int,
    RESIDUE_COUNT: pack_force_int,
    SIMULATION_COUNTER: pack_force_int,
    PARTICLE_POSITIONS: pack_vec3,
    PARTICLE_VELOCITIES: pack_vec3,
    PARTICLE_FORCES: pack_vec3,
    PARTICLE_FORCES_SYSTEM: pack_vec3,
    PARTICLE_ELEMENTS: pack_uint8,
    PARTICLE_RESIDUES: pack_uint32,
    BOND_PAIRS: pack_bond,
    BOND_ORDERS: pack_uint8,
    RESIDUE_CHAINS: pack_uint32,
    BOX_VECTORS: pack_vec3,
    USER_FORCES_INDEX: pack_uint32,
    USER_FORCES_SPARSE: pack_vec3,
}


def convert_dict_frame_to_grpc_frame(dict_frame) -> FrameData:
    grpc_frame = FrameData()

    for key, value in simplify_numpy(dict_frame).items():
        if isinstance(value, collections.abc.Sequence) and not isinstance(value, str):
            grpc_frame.arrays[key] = value
        else:
            grpc_frame.values[key] = value

    return grpc_frame


def convert_dict_state_to_dictionary_change(dict_state) -> DictionaryChange:
    return DictionaryChange(
        updates=dict_state["updates"],
        removals=dict_state["removals"],
    )


def convert_GetFrameResponse_to_framedata2(response: GetFrameResponse) -> FrameData2:
    frame = FrameData2(convert_grpc_frame_to_dict_frame(FrameData(response.frame)))
    frame.frame_index = response.frame_index
    return frame


def convert_framedata2_to_GetFrameResponse(frame: FrameData2) -> GetFrameResponse:
    return GetFrameResponse(
        frame_index=frame.frame_index,
        frame=convert_dict_frame_to_grpc_frame(frame).raw,
    )


def convert_grpc_frame_to_dict_frame(grpc_frame) -> dict[str, Any]:
    return unpack_dict_frame(pack_grpc_frame(grpc_frame))


def pack_grpc_frame(frame: FrameData) -> dict[str, Any]:
    data = {}

    for key in frame.value_keys:
        converter = converters.get(key, pack_identity)
        data[key] = converter.pack(frame.values[key])

    for key in frame.array_keys:
        converter = converters.get(key, pack_force_list)
        data[key] = converter.pack(frame.arrays[key])

    return data


def pack_dict_frame(frame: dict) -> dict[str, Any]:
    packed = {}

    for key in frame:
        converter = converters.get(key, pack_identity)
        packed[key] = converter.pack(frame[key])

    return packed


def unpack_dict_frame(frame: dict) -> dict[str, Any]:
    unpacked = {}

    for key in frame:
        converter = converters.get(key, pack_identity)
        unpacked[key] = converter.unpack(frame[key])

    return unpacked
