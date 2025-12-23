from dataclasses import dataclass
from inspect import getmembers
from typing import Any

import numpy as np
import numpy.typing as npt

from .frame_dict import frame_dict_packer, merge_frame_dicts, FrameDict
from . import keys

EnumArray = npt.NDArray[np.uint8]
IndexArray = npt.NDArray[np.uint32]
FloatArray = npt.NDArray[np.float32]
StringArray = list[str]


class MissingDataError(KeyError):
    """
    A shortcut does not contain data to return.
    """

    pass


@dataclass
class _Shortcut:
    key: str

    def make_property(shortcut) -> property:
        def getter(self: "FrameData"):
            try:
                return self.frame_dict[shortcut.key]
            except KeyError:
                raise MissingDataError(shortcut.key) from None

        def setter(self: "FrameData", value):
            self.frame_dict[shortcut.key] = value

        return property(fget=getter, fset=setter)


def _shortcut(key: str) -> Any:
    return _Shortcut(key)


def replace_shortcuts(cls):
    for name, attribute in getmembers(
        cls, lambda attribute: isinstance(attribute, _Shortcut)
    ):
        setattr(cls, name, attribute.make_property())
    return cls


@replace_shortcuts
class FrameData:
    """
    Wrapper around a frame dict, providing shortcut for common fields and convenience methods.
    """

    frame_dict: FrameDict

    frame_index: int = _shortcut(keys.FRAME_INDEX)

    box_vectors: FloatArray = _shortcut(keys.BOX_VECTORS)

    bond_pairs: IndexArray = _shortcut(keys.BOND_PAIRS)
    bond_orders: EnumArray = _shortcut(keys.BOND_ORDERS)

    particle_count: int = _shortcut(keys.PARTICLE_COUNT)
    particle_positions: FloatArray = _shortcut(keys.PARTICLE_POSITIONS)
    particle_velocities: FloatArray = _shortcut(keys.PARTICLE_VELOCITIES)
    particle_forces: FloatArray = _shortcut(keys.PARTICLE_FORCES)
    particle_forces_system: FloatArray = _shortcut(keys.PARTICLE_FORCES_SYSTEM)
    particle_elements: EnumArray = _shortcut(keys.PARTICLE_ELEMENTS)
    particle_names: StringArray = _shortcut(keys.PARTICLE_NAMES)
    particle_residues: IndexArray = _shortcut(keys.PARTICLE_RESIDUES)

    residue_count: int = _shortcut(keys.RESIDUE_COUNT)
    residue_names: StringArray = _shortcut(keys.RESIDUE_NAMES)
    residue_ids: StringArray = _shortcut(keys.RESIDUE_IDS)
    residue_chains: IndexArray = _shortcut(keys.RESIDUE_CHAINS)

    chain_count: int = _shortcut(keys.CHAIN_COUNT)
    chain_names: StringArray = _shortcut(keys.CHAIN_NAMES)

    kinetic_energy: float = _shortcut(keys.KINETIC_ENERGY)
    potential_energy: float = _shortcut(keys.POTENTIAL_ENERGY)
    system_temperature: float = _shortcut(keys.SYSTEM_TEMPERATURE)

    user_energy: float = _shortcut(keys.USER_ENERGY)
    user_work_done: float = _shortcut(keys.USER_WORK_DONE)
    user_forces_sparse: FloatArray = _shortcut(keys.USER_FORCES_SPARSE)
    user_forces_index: IndexArray = _shortcut(keys.USER_FORCES_INDEX)

    simulation_time: float = _shortcut(keys.SIMULATION_TIME)
    simulation_counter: float = _shortcut(keys.SIMULATION_COUNTER)
    simulation_exception: float = _shortcut(keys.SIMULATION_EXCEPTION)

    server_timestamp: float = _shortcut(keys.SERVER_TIMESTAMP)

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

    def copy(self):
        """
        Return a shallow copy this FrameData; an independent container of uncloned data.
        """
        return FrameData(self.frame_dict.copy())

    def update(self, other: "FrameData"):
        """
        Update this FrameData with another containing newer values. Ignores the previous data when frame index equals 0
        (indicating frame reset).
        """
        self.frame_dict = merge_frame_dicts(self.frame_dict, other.frame_dict)

    def __init__(self, frame_dict: FrameDict | None = None):
        self.frame_dict = frame_dict or {}

    def __repr__(self):
        if not self.frame_dict:
            return "<FrameData (empty)>"

        def format_key(key, value):
            if isinstance(value, np.ndarray):
                return f"{key}{list(value.shape)}"
            if isinstance(value, list):
                return f"{key}[{len(value)}]"
            return key

        key_text = ", ".join(
            format_key(key, self.frame_dict[key])
            for key in sorted(self.frame_dict.keys())
        )

        return f"<FrameData of {key_text}>"

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
