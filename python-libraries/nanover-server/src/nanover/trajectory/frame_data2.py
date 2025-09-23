from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt

from nanover.trajectory.frame_data import PARTICLE_POSITIONS

FrameDict = dict[str, Any]


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


class FrameData(metaclass=_FrameDataMeta):
    _shortcuts: dict[str, _Shortcut]

    def __init__(self, frame_dict: FrameDict):
        self.frame_dict = frame_dict

    particle_positions: npt.NDArray[np.float32] = _Shortcut(key=PARTICLE_POSITIONS)
