from abc import ABCMeta, abstractmethod
from dataclasses import dataclass
import copy

from typing import (
    Hashable,
    Any,
    Self,
    overload,
    Literal,
    TypeVar,
    Sequence,
    Iterable,
    Iterator,
)

import numpy as np
from openmm.unit import unit as omunit


# Convience for users to allow for wildcard importing without pulling abstract classes.
__all__ = ("Scalar", "Distance", "Angle", "Dihedral", "MultiMeasure")

MeasureKey = Hashable


@dataclass(init=False)
class BaseMeasure(metaclass=ABCMeta):
    name: str
    value: float
    unit: str

    def __init__(
        self, name: str, value: float, unit: omunit.Unit | str | None = None
    ) -> None:
        self.name = name
        self.value = value
        if isinstance(unit, omunit.Unit):
            self.unit = unit.get_symbol()
        else:
            self.unit = unit or ""

    @abstractmethod
    def __str__(self) -> str: ...

    @property
    @abstractmethod
    def key(self) -> MeasureKey:
        """Returns hashable datastructure to be used as keys for mapping objects."""
        ...

    @overload
    def _to_comparible_type(self, value: Any, strict: Literal[True]) -> float: ...

    @overload
    def _to_comparible_type(
        self, value: Any, strict: Literal[False]
    ) -> float | None: ...

    @overload
    def _to_comparible_type(self, value: Self | float | int) -> float: ...

    def _to_comparible_type(self, value: Any, strict: bool = True) -> float | None:
        """Converts `value` to a type usable for comparison to current `Measure` value."""
        if isinstance(value, type(self)):
            return value.value
        elif isinstance(value, (int, float)):
            return value

        if strict:
            raise TypeError(f"Could not compare {value} of type `{type(value)}`")
        return None

    def __eq__(self, value: Any) -> bool:
        """Compares the value of the measure to the given `value`."""
        return self._to_comparible_type(value, False) == self.value

    def __ne__(self, value: Any) -> bool:
        return not self == value

    def __lt__(self, value: Any) -> bool:
        return value is not None and self.value < self._to_comparible_type(value)

    def __gt__(self, value: Any) -> bool:
        return value is not None and self.value > self._to_comparible_type(value)

    def __le__(self, value: Any) -> bool:
        return value is not None and self.value <= self._to_comparible_type(value)

    def __ge__(self, value: Any) -> bool:
        return value is not None and self.value >= self._to_comparible_type(value)

    def update(self, value: Self | float | int) -> Self:
        if (new_val := self._to_comparible_type(value, False)) is None:
            raise TypeError(
                f"Can only update measurements with numeric values, not {type(value)}."
            )

        self.value = new_val
        return self

    @abstractmethod
    def to_fields(self) -> tuple[Any, ...]: ...


class Scalar(BaseMeasure):
    """Measurement class to handle scalar value measures."""

    def __str__(self) -> str:
        unit = f" {self.unit if self.unit is not None else ''}"
        return f"<(Scalar) {self.name}: {self.value}{unit}>"

    def __eq__(self, value):
        if isinstance(value, Scalar):
            return self.name == value.name and super().__eq__(value)
        return super().__eq__(value)

    @property
    def key(self) -> MeasureKey:
        return self.name  # Nothing else to use for hashing here.

    def to_fields(self) -> tuple[Any, ...]:
        return self.name, self.value, self.unit


class Distance(BaseMeasure):
    """Measurement class to handle distances between atoms."""

    atom1: int  # TODO ensure idx are correct/different atoms + >0
    atom2: int

    def __init__(
        self,
        name: str,
        atom1_index: int,
        atom2_index: int,
        distance: float,  # TODO distance > 0
        unit: str = "nm",  # TODO check unit + auto conversion
    ) -> None:
        super().__init__(name, distance, unit)
        self.atom1, self.atom2 = atom1_index, atom2_index

    def __str__(self) -> str:
        unit = f" {self.unit if self.unit is not None else ''}"
        return f"<(Distance) {self.name}: {self.value}{unit}>"

    def __eq__(self, value):
        if isinstance(value, Distance):
            return (self.atom1, self.atom2) == (
                value.atom1,
                value.atom2,
            ) and super().__eq__(value)
        return super().__eq__(value)

    @property
    def key(self) -> MeasureKey:
        return (self.atom1, self.atom2)

    def to_fields(self) -> tuple[Any, ...]:
        return self.name, [self.atom1, self.atom2], self.value, self.unit

    @property
    def distance(self) -> float:
        return self.value

    @distance.setter
    def distance(self, distance: float) -> None:
        self.value = distance


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
        unit = "rad" if radians else "deg"
        super().__init__(name, angle, unit)
        self.atom1, self.atom2, self.atom3 = atom1_index, atom2_index, atom3_index

    def __str__(self) -> str:
        unit = f" {self.unit if self.unit is not None else ''}"
        return f"<(Angle) {self.name}: {self.value}{unit}>"

    def __eq__(self, value):
        if isinstance(value, Angle):
            return (self.atom1, self.atom2, self.atom3) == (
                value.atom1,
                value.atom2,
                value.atom3,
            ) and super().__eq__(value)
        return super().__eq__(value)

    @property
    def key(self) -> MeasureKey:
        return (self.atom1, self.atom2, self.atom3)

    def to_fields(self) -> tuple[Any, ...]:
        return self.name, [self.atom1, self.atom2, self.atom3], self.value, self.unit

    @property
    def angle(self) -> float:
        return self.value

    @angle.setter
    def angle(self, angle: float) -> None:
        self.update(angle)

    def update(self, value: Self | float | int) -> Self:
        periodicity = 2 * np.pi if self.unit == "rad" else 360.0
        value = self._to_comparible_type(value) % periodicity
        return super().update(value)

    def to_radians(self) -> Self:
        if self.unit == "rad":
            return self
        self.unit = "rad"
        self.value = np.deg2rad(self.value)

        return self

    def to_degrees(self) -> Self:
        if self.unit == "deg":
            return self
        self.unit = "deg"
        self.value = np.rad2deg(self.value)

        return self


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
        unit = "rad" if radians else "deg"
        super().__init__(name, angle, unit)
        self.atom1, self.atom2, self.atom3, self.atom4 = (
            atom1_index,
            atom2_index,
            atom3_index,
            atom4_index,
        )

    def __str__(self) -> str:
        unit = f" {self.unit if self.unit is not None else ''}"
        return f"<(Dihedral) {self.name}: {self.value}{unit}>"

    def __eq__(self, value):
        if isinstance(value, Dihedral):
            return (self.atom1, self.atom2, self.atom3, self.atom4) == (
                value.atom1,
                value.atom2,
                value.atom3,
                value.atom4,
            ) and super().__eq__(value)
        return super().__eq__(value)

    @property
    def key(self) -> MeasureKey:
        return (self.atom1, self.atom2, self.atom3, self.atom4)

    def to_fields(self) -> tuple[Any, ...]:
        return (
            self.name,
            [self.atom1, self.atom2, self.atom3, self.atom4],
            self.value,
            self.unit,
        )

    @property
    def dihedral(self) -> float:
        return self.value

    @dihedral.setter
    def dihedral(self, dihedral: float) -> None:
        self.update(dihedral)

    def update(self, value: Self | float | int) -> Self:
        periodicity = 2 * np.pi if self.unit == "rad" else 360.0
        value = self._to_comparible_type(value) % periodicity
        return super().update(value)

    def to_radians(self) -> Self:
        if self.unit == "rad":
            return self
        self.unit = "rad"
        self.value = np.deg2rad(self.value)

        return self

    def to_degrees(self) -> Self:
        if self.unit == "deg":
            return self
        self.unit = "deg"
        self.value = np.rad2deg(self.value)

        return self


_M = TypeVar("_M", bound=BaseMeasure)


class MultiMeasure(Sequence[_M]):
    def __init__(
        self, initial_measure: _M, all_values: Iterable[float], copy: bool = False
    ) -> None:
        """
        Constructs a new `MultiMeasure` class.

        :param: initial_measure: Template `Measure` from which all subsequent values will be derived using.
        :param: all_values: Iterable containing all values to draw from.
        :param: copy: Whether to return a new instance of the `Measure` type when iterating over the data.
        """
        super().__init__()
        if not isinstance(initial_measure, BaseMeasure):
            raise TypeError(
                f"`initial_measure` should be of type `BaseMeasure` not `{type(initial_measure)}`"
            )

        self._data = np.array(all_values)
        self._measure_template = initial_measure
        self.copy = copy

    def __len__(self) -> int:
        return len(self._data)

    def __iter__(self) -> Iterator[_M]:
        for el in self._data:
            current_measure = self._measure_template.update(el)
            if not self.copy:
                yield current_measure
            else:
                copy.copy(current_measure)
