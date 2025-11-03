import operator
import io

from typing import Iterable, Callable, Any

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


MeasureMap = dict[MeasureKey, BaseMeasure]


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
    measures: Iterable[BaseMeasure] | None, check_type: type | None = None
) -> MeasureMap:
    """Creates a `MeasureMap`, ensuring all elements are the same (`check_type`)."""
    map_ = {el.key: el for el in measures} if measures is not None else MeasureMap()

    if check_type is None and map_:
        check_type = type(next(iter(map_)))

    if not all(type(el) is check_type for el in map_.values()):
        raise TypeError(
            f"Inconsistent type, all elements should be of type `{check_type}`."
        )
    return map_


class MeasureCollection:
    """Container class to handle sets of related measurements."""

    def __init__(
        self,
        scalars: Iterable[Scalar] | None = None,
        distances: Iterable[Distance] | None = None,
        angles: Iterable[Angle] | None = None,
        dihedrals: Iterable[Dihedral] | None = None,
    ):
        self.scalars: MeasureMap = _create_unitype_measuremap(scalars, Scalar)
        self.distances: MeasureMap = _create_unitype_measuremap(distances, Distance)
        self.angles: MeasureMap = _create_unitype_measuremap(angles, Angle)
        self.dihedrals: MeasureMap = _create_unitype_measuremap(dihedrals, Dihedral)

        self._type_mapping: dict[type[BaseMeasure], Callable[[Any], MeasureMap]] = {
            Scalar: operator.attrgetter("scalars"),
            Distance: operator.attrgetter("distances"),
            Angle: operator.attrgetter("angles"),
            Dihedral: operator.attrgetter("dihedrals"),
        }

    def __repr__(self) -> str:
        with io.StringIO() as str_io:
            str_io.write(f"{type(self).__name__} containing:\n")
            m = {
                "Scalars": self.scalars,
                "Distances": self.distances,
                "Angles": self.angles,
                "Dihedrals": self.dihedrals,
            }

            for k, v in m.items():
                if not v:
                    continue
                str_io.write(f"  {len(v)} {k}; {', '.join(map(str, v.values()))}\n")

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

    def __getitem__(self, key: BaseMeasure | Any) -> BaseMeasure:
        """Returns relevant `key` from underlying `MeasureMap`s."""
        if (target_getter := self._type_mapping.get(type(key), None)) is None:
            raise KeyError(f"Invalid {key}, only accepts Measure types. ")
        if (value := target_getter(self).get(key.key, None)) is None:
            raise KeyError(f'Could not find "{key}" in collections.')

        return value

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
                np.array(data, dtype=np.float32)
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
