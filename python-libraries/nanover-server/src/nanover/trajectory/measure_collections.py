import operator
import io
import itertools

from typing import Iterable, Callable, Any, TypeVar, overload, Iterator

import numpy as np

from nanover.trajectory import keys
from nanover.trajectory.measure import (
    MeasureKey,
    BaseMeasure,
    Scalar,
    Distance,
    Angle,
    Dihedral,
)
from nanover.trajectory import FrameData
from nanover.trajectory.frame_dict import FrameDict


BM = TypeVar("BM", bound=BaseMeasure)
MeasureMap = dict[MeasureKey, BM]


FRAMEDATA_MEASURE_LABELS: dict[type[BaseMeasure], str] = {
    Scalar: keys.MEASURE_SCALAR,
    Distance: keys.MEASURE_DISTANCE,
    Angle: keys.MEASURE_ANGLE,
    Dihedral: keys.MEASURE_DIHEDRAL,
}
FRAMEDATA_MEASURE_FIELD_KEYS: dict[type[BaseMeasure], tuple[str, ...]] = {
    Scalar: ("name", "value", "unit"),
    Distance: ("name", "atom_indices", "value", "unit"),
    Angle: ("name", "atom_indices", "value", "unit"),
    Dihedral: ("name", "atom_indices", "value", "unit"),
}


def _create_unitype_measuremap(
    measures: Iterable[BM] | None, check_type: type | None = None
) -> MeasureMap[BM]:
    """Creates a `MeasureMap`, ensuring all elements are the same (`check_type`)."""
    map_ = {el.key: el for el in measures} if measures is not None else MeasureMap()

    if check_type is None and map_:
        check_type = type(next(iter(map_)))

    if not all(type(el) is check_type for el in map_.values()):
        raise TypeError(
            f"Inconsistent type, all elements should be of type `{check_type}`."
        )
    return map_


def _flatten(fields: Iterable[Iterable[Any] | Any]) -> Iterable[Any]:
    """Flattens and iterable of values (containing list/array) into a single level of values."""
    for el in fields:
        if isinstance(el, (np.ndarray, list)):
            yield from el
        else:
            yield el


def _measures_from_framedata(frame: FrameData, measure_type: type[BM]) -> Iterable[BM]:
    if (
        measure_type not in FRAMEDATA_MEASURE_LABELS
        or measure_type not in FRAMEDATA_MEASURE_FIELD_KEYS
    ):
        raise TypeError(f"Unrecognised measure type: `{measure_type}`.")

    fields = [
        f"{FRAMEDATA_MEASURE_LABELS.get(measure_type)}.{el}"
        for el in FRAMEDATA_MEASURE_FIELD_KEYS.get(measure_type, [])
    ]
    field_iterator = itertools.zip_longest(
        *[frame.frame_dict.get(f, []) for f in fields], fillvalue=None
    )

    for values in field_iterator:
        yield measure_type(*_flatten(values))


class MeasureCollection:
    """Container class to handle sets of related measurements."""

    def __init__(
        self,
        scalars: Iterable[Scalar] | None = None,
        distances: Iterable[Distance] | None = None,
        angles: Iterable[Angle] | None = None,
        dihedrals: Iterable[Dihedral] | None = None,
    ):
        self.scalars = _create_unitype_measuremap(scalars, Scalar)
        self.distances = _create_unitype_measuremap(distances, Distance)
        self.angles = _create_unitype_measuremap(angles, Angle)
        self.dihedrals = _create_unitype_measuremap(dihedrals, Dihedral)

        self._type_mapping: dict[
            type[BaseMeasure], Callable[[MeasureCollection], MeasureMap]
        ] = {
            Scalar: operator.attrgetter("scalars"),
            Distance: operator.attrgetter("distances"),
            Angle: operator.attrgetter("angles"),
            Dihedral: operator.attrgetter("dihedrals"),
        }

    @classmethod
    def from_frame(cls, frame: FrameData) -> "MeasureCollection":
        scalars = _measures_from_framedata(frame, Scalar)
        distances = _measures_from_framedata(frame, Distance)
        angles = _measures_from_framedata(frame, Angle)
        dihedrals = _measures_from_framedata(frame, Dihedral)

        return cls(
            scalars=scalars, distances=distances, angles=angles, dihedrals=dihedrals
        )

    def __repr__(self) -> str:
        with io.StringIO() as str_io:
            str_io.write(f"{type(self).__name__} containing:\n")
            m: dict[str, MeasureMap] = {
                "Scalars": self.scalars,
                "Distances": self.distances,
                "Angles": self.angles,
                "Dihedrals": self.dihedrals,
            }

            for k, v in m.items():
                if not v:
                    continue
                str_io.write(f"\t{len(v)} {k}; {', '.join(map(str, v.values()))}\n")

            out_str = str_io.getvalue()
        return out_str

    def __str__(self) -> str:
        num_measures = [
            len(el)
            for el in (self.scalars, self.distances, self.angles, self.dihedrals)
        ]
        return "{} containing: {} scalar, {} distance, {} angle, and {} dihedral measurements.".format(
            type(self).__name__, *num_measures
        )

    def __iter__(self) -> Iterator[BaseMeasure]:
        all_measures = itertools.chain.from_iterable(
            el(self).values() for el in self._type_mapping.values()
        )
        yield from all_measures

    def _get_measure_from_name(self, name: str) -> BaseMeasure | None:
        """Finds the first corresponding `Measure` (if it exists) from the underlying `MeasureMaps`.

        Will search from the list of `Scalars, Distances, Angles, Dihedrals` in that order.
        """
        for mapping in self._type_mapping.values():
            for measure in mapping(self).values():
                if measure.name == name:
                    return measure

        return None

    def __getitem__(self, key: Any) -> BaseMeasure:
        # Check if looking for measure based on name first.
        if isinstance(key, str):
            if (value := self._get_measure_from_name(key)) is None:
                raise KeyError(f'Could not find "{key}" in collections.')
            return value

        # Now check for exact match as a provided measure.
        if isinstance(key, BaseMeasure):
            if (target_getter := self._type_mapping.get(type(key), None)) is None:
                raise KeyError(f"Invalid {key}, only accepts `Measure` types or `str`.")
            if (value := target_getter(self).get(key.key, None)) is None:
                raise KeyError(f'Could not find "{key}" in collections.')
            return value

        # Lastly check if provided `key` is set of atom indices.
        if isinstance(key, Iterable):
            key = tuple(key)
            for mapping in self._type_mapping.values():
                if (value := mapping(self).get(key, None)) is not None:
                    return value
            else:
                raise KeyError(
                    f'Could not find measure with matching parameters of "{key}"'
                )

        # Nothing could be found.
        raise KeyError(f'Could not find measure with given data: "{key}"')

    @overload
    def get_measure(self, measure: str) -> BaseMeasure:
        """Returns relevant measure by matching against its name."""
        ...

    @overload
    def get_measure(self, measure: tuple[int]) -> BaseMeasure:
        """
        Returns relevant measure by matching against the atom indices of a measure.
        Only applicable for geometric measurement types, e.g. for a `Distance`.
        """
        ...

    @overload
    def get_measure(self, measure: BaseMeasure) -> BaseMeasure:
        """Returns relevant measure by matching against the provided `key`."""
        ...

    def get_measure(self, measure: str | tuple[int] | BaseMeasure) -> BaseMeasure:
        """
        Returns relevant "measure" from the collection.

        If `measure` is a:
        - `str`, will return the first matching measure with the same "name".
        - `tuple[int]`, will return the measure with the given set of (ordered) atom indices.

        All matching will check for identity from the sets of: Scalars, Distances, Angles, Dihedrals (respectively).
        """
        return self[measure]

    def _measure_iterator(
        self, measurements: Iterable[BaseMeasure]
    ) -> Iterable[tuple[BaseMeasure, MeasureMap]]:
        for measure in measurements:
            target_getter = self._type_mapping.get(type(measure), None)
            if target_getter is None:
                continue

            yield measure, target_getter(self)

    def update(self, measurements: Iterable[BaseMeasure] | BaseMeasure) -> None:
        """Updates existing stored measurements with `measurements`."""
        if not isinstance(measurements, Iterable):
            measurements = [measurements]

        for measure, existing_measures in self._measure_iterator(measurements):
            if (m := existing_measures.get(measure.key, None)) is not None:
                m.update(measure)
            else:
                existing_measures[measure.key] = measure

    def remove(self, measurements: Iterable[BaseMeasure] | BaseMeasure) -> None:
        """Removes measurement(s) from stored list, if they exist."""
        if not isinstance(measurements, Iterable):
            measurements = [measurements]

        for measure, existing_measures in self._measure_iterator(measurements):
            existing_measures.pop(measure.key, None)

    def _add_measureset_to_framedict(
        self, measurements: MeasureMap, frame_dict: FrameDict
    ) -> None:
        measure_type = type(next(iter(measurements.values())))
        if (field_keys := FRAMEDATA_MEASURE_FIELD_KEYS.get(measure_type)) is None:
            raise KeyError(
                f"Unsupported measurement type ({measure_type}), "
                f"only {', '.join(map(str, FRAMEDATA_MEASURE_FIELD_KEYS))} are supported."
            )

        measure_fields = (el.to_fields() for el in measurements.values())
        for field_name, data in zip(
            field_keys,
            zip(*measure_fields),
        ):
            arr = (
                np.array(data, dtype=np.uint32)
                if field_name.endswith("atom_indices")
                else np.array(data)
            )
            frame_dict.update(
                {
                    f"{FRAMEDATA_MEASURE_LABELS.get(measure_type)}.{field_name}": arr,
                }
            )

    def add_to_framedata(self, framedata: FrameData) -> FrameData:
        """Adds currently stored measurements to the given `framedata`."""
        new_data: FrameDict = {}

        for s in (self.scalars, self.distances, self.angles, self.dihedrals):
            if not s:
                continue  # Skip empty mappings.
            self._add_measureset_to_framedict(s, new_data)

        framedata.update(FrameData(new_data))
        return framedata
