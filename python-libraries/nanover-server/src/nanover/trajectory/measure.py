from abc import ABCMeta, abstractmethod
from dataclasses import dataclass
from enum import IntEnum, auto

from numbers import Real
from typing import TypeVar, Hashable, Any, Self


# Convience for users to allow for wildcard importing without pulling abstract classes.
__all__ = ("Measure", "Distance", "Angle", "Dihedral")


MeasureKey = TypeVar("MeasureKey", bound=Hashable)


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

    @property
    @abstractmethod
    def key(self) -> MeasureKey:
        """Returns hashable datastructure to be used as keys for mapping objects."""
        ...

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
        self.value = value.value

    @abstractmethod
    def to_fields(self) -> tuple[Any, ...]: ...


class Measure(BaseMeasure):
    """Measurement class to handle scalar value measures."""

    @property
    def key(self) -> MeasureKey:
        return self.name  # Nothing else to use for hashing here.

    def to_fields(self) -> tuple[Any, ...]:
        return self.name, self.value


class Distance(BaseMeasure):
    """Measurement class to handle distances between atoms."""

    atom1: int
    atom2: int

    @property
    def key(self) -> MeasureKey:
        return (self.atom1, self.atom2)

    def to_fields(self) -> tuple[Any, ...]:
        return self.name, [self.atom1, self.atom2], self.value


class Angle(BaseMeasure):  # TODO handle angles in radians/degrees + periodicity
    """Measurement class to handle angles between atoms."""

    atom1: int
    atom2: int
    atom3: int

    @property
    def key(self) -> MeasureKey:
        return (self.atom1, self.atom2, self.atom3)

    def to_fields(self) -> tuple[Any, ...]:
        return self.name, [self.atom1, self.atom2, self.atom3], self.value


class Dihedral(BaseMeasure):
    """Measurement class to handle dihedrals (torsions) between atoms."""

    atom1: int
    atom2: int
    atom3: int
    atom4: int

    @property
    def key(self) -> MeasureKey:
        return (self.atom1, self.atom2, self.atom3, self.atom4)

    def to_fields(self) -> tuple[Any, ...]:
        return self.name, [self.atom1, self.atom2, self.atom3, self.atom4], self.value
