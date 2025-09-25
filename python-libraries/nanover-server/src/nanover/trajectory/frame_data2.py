from dataclasses import dataclass
from inspect import getmembers
from typing import Any

import numpy as np
import numpy.typing as npt

from nanover.trajectory.frame_data import (
    PARTICLE_POSITIONS,
    BOND_PAIRS,
    BOND_ORDERS,
    PARTICLE_VELOCITIES,
    PARTICLE_FORCES,
    PARTICLE_FORCES_SYSTEM,
    PARTICLE_ELEMENTS,
    PARTICLE_NAMES,
    PARTICLE_RESIDUES,
    PARTICLE_COUNT,
    RESIDUE_NAMES,
    RESIDUE_IDS,
    RESIDUE_CHAINS,
    RESIDUE_COUNT,
    CHAIN_NAMES,
    CHAIN_COUNT,
    KINETIC_ENERGY,
    POTENTIAL_ENERGY,
    USER_ENERGY,
    SYSTEM_TEMPERATURE,
    USER_FORCES_SPARSE,
    USER_FORCES_INDEX,
    USER_WORK_DONE,
    BOX_VECTORS,
    SIMULATION_TIME,
    SIMULATION_COUNTER,
    SIMULATION_EXCEPTION,
    SERVER_TIMESTAMP,
    MissingDataError,
    FRAME_INDEX,
)

FrameDict = dict[str, Any]

EnumArray = npt.NDArray[np.uint8]
IndexArray = npt.NDArray[np.uint32]
FloatArray = npt.NDArray[np.float32]
StringArray = list[str]


@dataclass(kw_only=True)
class _Shortcut:
    key: str

    def make_property(shortcut):
        def getter(self: "FrameData"):
            try:
                return self.frame_dict[shortcut.key]
            except KeyError:
                raise MissingDataError(shortcut.key) from None

        def setter(self: "FrameData", value):
            self.frame_dict[shortcut.key] = value

        return property(fget=getter, fset=setter)


def _shortcut(*, key: str) -> Any:
    return _Shortcut(key=key)


def merge_frame_dicts(a: dict, b: dict, ignore_reset=False):
    merged = {}
    if b.get(FRAME_INDEX, None) != 0 or ignore_reset:
        merged.update(a)
    merged.update(b)
    return merged


def replace_shortcuts(cls):
    for name, attribute in getmembers(
        cls, lambda attribute: isinstance(attribute, _Shortcut)
    ):
        setattr(cls, name, attribute.make_property())
    return cls


@replace_shortcuts
class FrameData:
    _shortcuts: dict[str, _Shortcut]

    def __init__(self, frame_dict: FrameDict | None = None):
        self.frame_dict = frame_dict or {}

    def __bool__(self):
        return bool(self.frame_dict)

    def __contains__(self, key: str):
        return key in self.frame_dict

    def __getitem__(self, key: str):
        return self.frame_dict[key]

    def __setitem__(self, key: str, value: Any):
        self.frame_dict[key] = value

    def copy(self):
        return FrameData(self.frame_dict.copy())

    def update(self, other: "FrameData", ignore_reset=False):
        self.frame_dict = merge_frame_dicts(
            self.frame_dict, other.frame_dict, ignore_reset=ignore_reset
        )

    frame_index: int = _shortcut(key=FRAME_INDEX)

    box_vectors: FloatArray = _shortcut(key=BOX_VECTORS)

    bond_pairs: IndexArray = _shortcut(key=BOND_PAIRS)
    bond_orders: EnumArray = _shortcut(key=BOND_ORDERS)

    particle_count: int = _shortcut(key=PARTICLE_COUNT)
    particle_positions: FloatArray = _shortcut(key=PARTICLE_POSITIONS)
    particle_velocities: FloatArray = _shortcut(key=PARTICLE_VELOCITIES)
    particle_forces: FloatArray = _shortcut(key=PARTICLE_FORCES)
    particle_forces_system: FloatArray = _shortcut(key=PARTICLE_FORCES_SYSTEM)
    particle_elements: EnumArray = _shortcut(key=PARTICLE_ELEMENTS)
    particle_names: StringArray = _shortcut(key=PARTICLE_NAMES)
    particle_residues: IndexArray = _shortcut(key=PARTICLE_RESIDUES)

    residue_count: int = _shortcut(key=RESIDUE_COUNT)
    residue_names: StringArray = _shortcut(key=RESIDUE_NAMES)
    residue_ids: StringArray = _shortcut(key=RESIDUE_IDS)
    residue_chains: IndexArray = _shortcut(key=RESIDUE_CHAINS)

    chain_count: int = _shortcut(key=CHAIN_COUNT)
    chain_names: StringArray = _shortcut(key=CHAIN_NAMES)

    kinetic_energy: float = _shortcut(key=KINETIC_ENERGY)
    potential_energy: float = _shortcut(key=POTENTIAL_ENERGY)
    system_temperature: float = _shortcut(key=SYSTEM_TEMPERATURE)

    user_energy: float = _shortcut(key=USER_ENERGY)
    user_work_done: float = _shortcut(key=USER_WORK_DONE)
    user_forces_sparse: FloatArray = _shortcut(key=USER_FORCES_SPARSE)
    user_forces_index: IndexArray = _shortcut(key=USER_FORCES_INDEX)

    simulation_time: float = _shortcut(key=SIMULATION_TIME)
    simulation_counter: float = _shortcut(key=SIMULATION_COUNTER)
    simulation_exception: float = _shortcut(key=SIMULATION_EXCEPTION)

    server_timestamp: float = _shortcut(key=SERVER_TIMESTAMP)
