from abc import ABCMeta, abstractmethod
from dataclasses import dataclass
from enum import IntEnum, auto

from typing import Hashable, Any, Self

from openmm.unit import unit as omunit


# Convience for users to allow for wildcard importing without pulling abstract classes.
__all__ = ("Measure", "Distance", "Angle", "Dihedral")


MeasureKey = Hashable


class UpdateStatus(IntEnum):
    NEW = auto()
    UPDATED = auto()
    DELETE = auto()


@dataclass(init=False)
class BaseMeasure(metaclass=ABCMeta):
    name: str
    value: float
    unit: str
    _update_status: UpdateStatus

    def __init__(
        self, name: str, value: float, unit: omunit.Unit | str | None = None
    ) -> None:
        self.name = name
        self.value = value
        if isinstance(unit, omunit.Unit):
            self.unit = unit.get_symbol()
        else:
            self.unit = unit or ""
        self._update_status = UpdateStatus.NEW

    @abstractmethod
    def __str__(self) -> str: ...

    @property
    @abstractmethod
    def key(self) -> MeasureKey:
        """Returns hashable datastructure to be used as keys for mapping objects."""
        ...

    def _to_comparible_type(self, value: Any) -> float:
        """Converts `value` to a type usable for comparison to current `Measure` value."""
        if isinstance(value, type(self)):
            return value.value
        elif isinstance(value, (int, float)):
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

    def __str__(self) -> str:
        unit = f" {self.unit if self.unit is not None else ''}"
        return f"{self.name}: {self.value}{unit}"

    @property
    def key(self) -> MeasureKey:
        return self.name  # Nothing else to use for hashing here.

    def to_fields(self) -> tuple[Any, ...]:
        return self.name, self.value


class Distance(BaseMeasure):
    """Measurement class to handle distances between atoms."""

    atom1: int
    atom2: int

    def __init__(
        self, name: str, atom1_index: int, atom2_index: int, distance: float, unit=None
    ) -> None:
        super().__init__(name, distance, unit)
        self.atom1, self.atom2 = atom1_index, atom2_index

    def __str__(self) -> str:
        unit = f" {self.unit if self.unit is not None else ''}"
        return f"(Distance) {self.name}: {self.value}{unit}"

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

    def __init__(
        self,
        name: str,
        atom1_index: int,
        atom2_index: int,
        atom3_index: int,
        angle: float,
        radians: bool = False,
    ) -> None:
        unit = "radians" if radians else "degrees"
        super().__init__(name, angle, unit)
        self.atom1, self.atom2, self.atom3 = atom1_index, atom2_index, atom3_index

    def __str__(self) -> str:
        unit = f" {self.unit if self.unit is not None else ''}"
        return f"(Angle) {self.name}: {self.value}{unit}"

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

    def __init__(
        self,
        name: str,
        atom1_index: int,
        atom2_index: int,
        atom3_index: int,
        atom4_index: int,
        angle: float,
        radians: bool = False,
    ) -> None:
        unit = "radians" if radians else "degrees"
        super().__init__(name, angle, unit)
        self.atom1, self.atom2, self.atom3, self.atom4 = (
            atom1_index,
            atom2_index,
            atom3_index,
            atom4_index,
        )

    def __str__(self) -> str:
        unit = f" {self.unit if self.unit is not None else ''}"
        return f"(Dihedral) {self.name}: {self.value}{unit}"

    @property
    def key(self) -> MeasureKey:
        return (self.atom1, self.atom2, self.atom3, self.atom4)

    def to_fields(self) -> tuple[Any, ...]:
        return self.name, [self.atom1, self.atom2, self.atom3, self.atom4], self.value
