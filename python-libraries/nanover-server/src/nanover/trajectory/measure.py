from abc import ABCMeta, abstractmethod
from dataclasses import dataclass
from enum import IntEnum, auto

from numbers import Real
from typing import Any, Self


# Convience for users to allow for wildcard importing without pulling abstract classes.
__all__ = ("Measure", "Distance", "Angle", "Dihedral")


class UpdateStatus(IntEnum):
    NEW = auto()
    UPDATED = auto()
    DELETE = auto()


@dataclass(init=False)
class BaseMeasure(metaclass=ABCMeta):
    name: str
    value: float
    _update_status: UpdateStatus

    def __init__(self, name: str, value: float) -> None:
        self.name = name
        self.value = value
        self._update_status = UpdateStatus.NEW

    @abstractmethod
    def __hash__(self) -> int: ...

    def _to_comparible_type(self, value: Any) -> bool:
        """Converts `value` to a type usable for comparison to current `Measure` value."""
        if isinstance(value, type(self)):
            return value.value
        elif isinstance(value, Real):
            return value
        else:
            raise TypeError(
                f"Only comparisons of {type(self)} to numeric or other "
                f"{type(self)} types is supported, not {type(self)} "
                f"and {type(value)}. "
            )

    def __eq__(self, value: Any) -> bool:
        """Compares the value of the measure to the given `value`."""
        return self._to_comparible_type(value) == self.value

    def __ne__(self, value: Any) -> bool:
        return not self == value

    def __lt__(self, value: Any) -> bool:
        return self.value < self._to_comparible_type(value)

    def __gt__(self, value: Any) -> bool:
        return self.value > self._to_comparible_type(value)

    def __le__(self, value: Any) -> bool:
        return self.value <= self._to_comparible_type(value)

    def __ge__(self, value: Any) -> bool:
        return self.value >= self._to_comparible_type(value)

    def update(self, value: Self) -> None:
        if not isinstance(value, (BaseMeasure)):
            raise TypeError(
                f"Can only update measurements with numeric values, not {type(value)}."
            )
        self._update_status = UpdateStatus.UPDATED
        self.value = value


class Measure(BaseMeasure):
    """Measurement class to handle scalar value measures."""

    def __hash__(self) -> int:
        return hash(self.name)  # Nothing else to use for hashing here.


class Distance(BaseMeasure):
    """Measurement class to handle distances between atoms."""

    atom1: int
    atom2: int

    def __hash__(self) -> int:
        return hash((self.atom1, self.atom2))


class Angle(BaseMeasure):
    """Measurement class to handle angles between atoms."""

    atom1: int
    atom2: int
    atom3: int

    def __hash__(self) -> int:
        return hash((self.atom1, self.atom2, self.atom3))


class Dihedral(BaseMeasure):
    """Measurement class to handle dihedrals (torsions) between atoms."""

    atom1: int
    atom2: int
    atom3: int
    atom4: int

    def __hash__(self) -> int:
        return hash((self.atom1, self.atom2, self.atom3, self.atom4))
