import operator
import io

from typing import Iterable, Callable, Any

from nanover.trajectory.measure import (
    MeasureKey,
    BaseMeasure,
    Measure,
    Distance,
    Angle,
    Dihedral,
)
from nanover.trajectory import FrameData
from nanover.trajectory.frame_dict import FrameDict


MeasureMap = dict[MeasureKey, BaseMeasure]


FRAMEDATA_MEASURE_FIELD_PREFIX = "measure"
FRAMEDATA_MEASURE_LABELS: dict[type[BaseMeasure], str] = {
    Measure: "scalar",
    Distance: "distance",
    Angle: "angle",
    Dihedral: "dihedral",
}
FRAMEDATA_MEASURE_FIELD_KEYS: dict[type[BaseMeasure], tuple[str, ...]] = {
    Measure: ("name", "value", "unit"),
    Distance: ("name", "atom_indices", "value", "unit"),
    Angle: ("name", "atom_indices", "value", "unit"),
    Dihedral: ("name", "atom_indices", "value", "unit"),
}


class MeasureCollection:
    """Container class to handle sets of related measurements."""

    def __init__(
        self,
        scalars: Iterable[Measure] | None = None,
        distances: Iterable[Distance] | None = None,
        angles: Iterable[Angle] | None = None,
        dihedrals: Iterable[Dihedral] | None = None,
    ):
        self.scalars: MeasureMap = (
            {el.key: el for el in scalars} if scalars is not None else MeasureMap()
        )
        self.distances: MeasureMap = (
            {el.key: el for el in distances} if distances is not None else MeasureMap()
        )
        self.angles: MeasureMap = (
            {el.key: el for el in angles} if angles is not None else MeasureMap()
        )
        self.dihedrals: MeasureMap = (
            {el.key: el for el in dihedrals} if dihedrals is not None else MeasureMap()
        )

        self._type_mapping: dict[type[BaseMeasure], Callable[[Any], MeasureMap]] = {
            Measure: operator.attrgetter("scalars"),
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

    def _measureset_to_tuples(self, measurements: MeasureMap) -> Iterable[tuple[Any]]:
        """Yields each element in set of `Measure`s as relevant FrameData parameters."""
        for el in measurements.values():
            yield el.to_fields()

    def _add_measureset_to_framedict(
        self, measurements: MeasureMap, frame_dict: FrameDict
    ) -> None:
        measure_type = type(next(iter(measurements.values())))
        if (field_keys := FRAMEDATA_MEASURE_FIELD_KEYS.get(measure_type)) is None:
            raise KeyError(
                f"Unsupported measurement type ({measure_type}), "
                f"only {', '.join(map(str, FRAMEDATA_MEASURE_FIELD_KEYS))} are supported."
            )

        for field_name, data in zip(
            field_keys,
            zip(*self._measureset_to_tuples(measurements)),
        ):
            frame_dict.update(
                {
                    f"{FRAMEDATA_MEASURE_FIELD_PREFIX}.{FRAMEDATA_MEASURE_LABELS.get(measure_type)}.{field_name}": list(
                        data
                    ),
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
