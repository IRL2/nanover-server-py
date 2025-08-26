from functools import partial
from typing import Iterable, Callable, Any

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


def pack_array(dtype: str, values: Iterable):
    return np.fromiter(values, dtype=dtype).tobytes()


pack_float32 = partial(pack_array, "f")
pack_uint32 = partial(pack_array, "I")
pack_uint8 = partial(pack_array, "B")

converters: dict[str, Callable[[Any], Any]] = {
    PARTICLE_COUNT: int,
    CHAIN_COUNT: int,
    RESIDUE_COUNT: int,
    SIMULATION_COUNTER: int,
    PARTICLE_POSITIONS: pack_float32,
    PARTICLE_ELEMENTS: pack_uint8,
    PARTICLE_RESIDUES: pack_uint32,
    BOND_PAIRS: pack_uint32,
    RESIDUE_CHAINS: pack_uint32,
    BOX_VECTORS: pack_float32,
}


def convert_frame(frame: FrameData):
    data = {}

    for key in frame.value_keys:
        converter = converters.get(key, lambda value: value)
        data[key] = converter(frame.values[key])

    for key in frame.array_keys:
        converter = converters.get(key, list)
        data[key] = converter(frame.arrays[key])

    return data
