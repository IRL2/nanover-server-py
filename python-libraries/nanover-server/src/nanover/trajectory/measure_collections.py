import operator

from typing import Iterable, Callable, Any, TypeVar

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
FRAMEDATA_MEASURE_FIELD_KEYS: dict[BaseMeasure, tuple[str, ...]] = {
    Measure: ("name", "value"),
    Distance: ("name", "atom_indices", "value"),
    Angle: ("name", "atom_indices", "value"),
    Dihedral: ("name", "atom_indices", "value"),
}


class MeasureCollection:
    def __init__(
        self,
        scalars: Iterable[Measure] | None = None,
        distances: Iterable[Distance] | None = None,
        angles: Iterable[Angle] | None = None,
        dihedrals: Iterable[Dihedral] | None = None,
    ):
        self.scalars: MeasureMap = (
            {el.key: el for el in scalars} if scalars is not None else set()
        )
        self.distances: MeasureMap = (
            {el.key: el for el in distances} if distances is not None else set()
        )
        self.angles: MeasureMap = (
            {el.key: el for el in angles} if angles is not None else set()
        )
        self.dihedrals: MeasureMap = (
            {el.key: el for el in dihedrals} if dihedrals is not None else set()
        )

        self._type_mapping: dict[BaseMeasure, Callable[[], MeasureMap]] = {
            Measure: operator.attrgetter("scalars"),
            Distance: operator.attrgetter("distances"),
            Angle: operator.attrgetter("angles"),
            Dihedral: operator.attrgetter("dihedrals"),
        }

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
            print(list(hash(el) for el in existing_measures))
            print(hash(measure.key))
            print(measure.key in existing_measures)
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
        for el in measurements:
            yield el.to_fields()

    def _add_measureset_to_framedict(
        self, measurements: set[BaseMeasure], frame_dict: FrameDict
    ) -> None:
        measure_type = type(next(iter(measurements)))
        if measure_type not in FRAMEDATA_MEASURE_FIELD_KEYS:
            raise KeyError(
                f"Unsupported measurement type ({measure_type}), "
                f"only {', '.join(FRAMEDATA_MEASURE_FIELD_KEYS)} are supported."
            )

        for field_name, data in zip(
            FRAMEDATA_MEASURE_FIELD_KEYS.get(measure_type),
            zip(*self._measureset_to_tuples(measurements)),
        ):
            frame_dict.update(
                {
                    f"{FRAMEDATA_MEASURE_FIELD_PREFIX}.{str(measure_type).lower()}.{field_name}": data,
                }
            )

    def add_to_framedata(self, framedata: FrameData) -> FrameData:
        """Adds currently stored measurements to the given `framedata`."""
        new_data: FrameDict = {}

        for s in (self.scalars, self.distances, self.angles, self.dihedrals):
            print(s)
            if not s:
                continue  # Skip empty mappings.
            self._add_measureset_to_framedict(s, new_data)

        framedata.update(FrameData(new_data))
        return framedata
