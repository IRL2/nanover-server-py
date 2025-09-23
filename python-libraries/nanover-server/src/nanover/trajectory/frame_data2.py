from dataclasses import dataclass
from typing import Any, no_type_check

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
        def get(self: "FrameData"):
            return self.frame_dict[shortcut.key]

        def set(self: "FrameData", value):
            self.frame_dict[shortcut.key] = value

        return property(fget=get, fset=set)


class _FrameDataMeta(type):
    """
    Metaclass that adds shortcuts to the :class:`FrameData` class.

    The shortcuts are defined as a tuple of :class:`_Shortcut` named tuples
    under the :attr:`_shortcuts` class attribute.
    """

    _shortcuts: dict[str, _Shortcut] = {}

    def __init__(cls, name, bases, nmspc):
        shortcuts = {}
        super().__init__(name, bases, nmspc)
        for attribute_name, attribute in nmspc.items():
            if isinstance(attribute, _Shortcut):
                shortcuts[attribute_name] = attribute
                setattr(cls, attribute_name, attribute.make_property())
        cls._shortcuts = shortcuts


@no_type_check
class FrameData(metaclass=_FrameDataMeta):
    _shortcuts: dict[str, _Shortcut]

    def __init__(self, frame_dict: FrameDict):
        self.frame_dict = frame_dict

    box_vectors: FloatArray = _Shortcut(key=BOX_VECTORS)

    bond_pairs: IndexArray = _Shortcut(key=BOND_PAIRS)
    bond_orders: EnumArray = _Shortcut(key=BOND_ORDERS)

    particle_positions: FloatArray = _Shortcut(key=PARTICLE_POSITIONS)
    particle_velocities: FloatArray = _Shortcut(key=PARTICLE_VELOCITIES)
    particle_forces: FloatArray = _Shortcut(key=PARTICLE_FORCES)
    particle_forces_system: FloatArray = _Shortcut(key=PARTICLE_FORCES_SYSTEM)
    particle_elements: EnumArray = _Shortcut(key=PARTICLE_ELEMENTS)
    particle_names: StringArray = _Shortcut(key=PARTICLE_NAMES)
    particle_residues: IndexArray = _Shortcut(key=PARTICLE_RESIDUES)
    particle_count: int = _Shortcut(key=PARTICLE_COUNT)

    residue_names: StringArray = _Shortcut(key=RESIDUE_NAMES)
    residue_ids: StringArray = _Shortcut(key=RESIDUE_IDS)
    residue_chains: IndexArray = _Shortcut(key=RESIDUE_CHAINS)
    residue_count: int = _Shortcut(key=RESIDUE_COUNT)

    chain_names: StringArray = _Shortcut(key=CHAIN_NAMES)
    chain_count: int = _Shortcut(key=CHAIN_COUNT)

    kinetic_energy: float = _Shortcut(key=KINETIC_ENERGY)
    potential_energy: float = _Shortcut(key=POTENTIAL_ENERGY)
    system_temperature: float = _Shortcut(key=SYSTEM_TEMPERATURE)

    user_energy: float = _Shortcut(key=USER_ENERGY)
    user_forces_sparse: FloatArray = _Shortcut(key=USER_FORCES_SPARSE)
    user_forces_index: IndexArray = _Shortcut(key=USER_FORCES_INDEX)
    user_work_done: float = _Shortcut(key=USER_WORK_DONE)

    simulation_time: float = _Shortcut(key=SIMULATION_TIME)
    simulation_counter: float = _Shortcut(key=SIMULATION_COUNTER)
    simulation_exception: float = _Shortcut(key=SIMULATION_EXCEPTION)

    server_timestamp: float = _Shortcut(key=SERVER_TIMESTAMP)
