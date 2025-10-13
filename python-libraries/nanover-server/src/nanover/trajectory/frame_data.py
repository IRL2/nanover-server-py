from dataclasses import dataclass
from inspect import getmembers
from typing import Any

import numpy as np
import numpy.typing as npt

from .convert import frame_dict_packer
from . import keys

FrameDict = dict[str, Any]

EnumArray = npt.NDArray[np.uint8]
IndexArray = npt.NDArray[np.uint32]
FloatArray = npt.NDArray[np.float32]
StringArray = list[str]


class MissingDataError(KeyError):
    """
    A shortcut does not contain data to return.
    """

    pass


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
    a_reset = a.get(keys.FRAME_INDEX, None) == 0
    b_reset = b.get(keys.FRAME_INDEX, None) == 0 and not ignore_reset

    merged = {}

    if not b_reset:
        merged.update(a)

    merged.update(b)

    if a_reset or b_reset:
        merged[keys.FRAME_INDEX] = 0

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

    @classmethod
    def empty(cls):
        return cls()

    @classmethod
    def from_dict(cls, frame_dict: FrameDict):
        """
        Return a new FrameData from a dict of unpacked data.
        """
        return cls(frame_dict)

    @classmethod
    def unpack_from_dict(cls, frame_dict: FrameDict):
        """
        Return a new FrameData from a dict of packed data.
        """
        return cls.from_dict(frame_dict_packer.unpack(frame_dict))

    def pack_to_dict(self):
        """
        Return a dict of packed data from this frame.
        """
        return frame_dict_packer.pack(self.frame_dict)

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

    def __delitem__(self, key: str):
        del self.frame_dict[key]

    def copy(self):
        return FrameData(self.frame_dict.copy())

    def update(self, other: "FrameData", ignore_reset=False):
        self.frame_dict = merge_frame_dicts(
            self.frame_dict, other.frame_dict, ignore_reset=ignore_reset
        )

    frame_index: int = _shortcut(key=keys.FRAME_INDEX)

    box_vectors: FloatArray = _shortcut(key=keys.BOX_VECTORS)

    bond_pairs: IndexArray = _shortcut(key=keys.BOND_PAIRS)
    bond_orders: EnumArray = _shortcut(key=keys.BOND_ORDERS)

    particle_count: int = _shortcut(key=keys.PARTICLE_COUNT)
    particle_positions: FloatArray = _shortcut(key=keys.PARTICLE_POSITIONS)
    particle_velocities: FloatArray = _shortcut(key=keys.PARTICLE_VELOCITIES)
    particle_forces: FloatArray = _shortcut(key=keys.PARTICLE_FORCES)
    particle_forces_system: FloatArray = _shortcut(key=keys.PARTICLE_FORCES_SYSTEM)
    particle_elements: EnumArray = _shortcut(key=keys.PARTICLE_ELEMENTS)
    particle_names: StringArray = _shortcut(key=keys.PARTICLE_NAMES)
    particle_residues: IndexArray = _shortcut(key=keys.PARTICLE_RESIDUES)

    residue_count: int = _shortcut(key=keys.RESIDUE_COUNT)
    residue_names: StringArray = _shortcut(key=keys.RESIDUE_NAMES)
    residue_ids: StringArray = _shortcut(key=keys.RESIDUE_IDS)
    residue_chains: IndexArray = _shortcut(key=keys.RESIDUE_CHAINS)

    chain_count: int = _shortcut(key=keys.CHAIN_COUNT)
    chain_names: StringArray = _shortcut(key=keys.CHAIN_NAMES)

    kinetic_energy: float = _shortcut(key=keys.KINETIC_ENERGY)
    potential_energy: float = _shortcut(key=keys.POTENTIAL_ENERGY)
    system_temperature: float = _shortcut(key=keys.SYSTEM_TEMPERATURE)

    user_energy: float = _shortcut(key=keys.USER_ENERGY)
    user_work_done: float = _shortcut(key=keys.USER_WORK_DONE)
    user_forces_sparse: FloatArray = _shortcut(key=keys.USER_FORCES_SPARSE)
    user_forces_index: IndexArray = _shortcut(key=keys.USER_FORCES_INDEX)

    simulation_time: float = _shortcut(key=keys.SIMULATION_TIME)
    simulation_counter: float = _shortcut(key=keys.SIMULATION_COUNTER)
    simulation_exception: float = _shortcut(key=keys.SIMULATION_EXCEPTION)

    server_timestamp: float = _shortcut(key=keys.SERVER_TIMESTAMP)
