import operator

from typing import Iterable, Callable

from nanover.trajectory.measure import BaseMeasure, Measure, Distance, Angle, Dihedral


class MeasureCollection:
    def __init__(
        self,
        scalars: Iterable[Measure] | None = None,
        distances: Iterable[Distance] | None = None,
        angles: Iterable[Angle] | None = None,
        dihedrals: Iterable[Dihedral] | None = None,
    ):
        self.scalars = set(scalars) if scalars is not None else set()
        self.distances = set(distances) if distances is not None else set()
        self.angles = set(angles) if angles is not None else set()
        self.dihedrals = set(dihedrals) if dihedrals is not None else set()

        self._type_mapping: dict[BaseMeasure, Callable[[], set[BaseMeasure]]] = {
            Measure: operator.attrgetter("scalars"),
            Distance: operator.attrgetter("distances"),
            Angle: operator.attrgetter("angles"),
            Dihedral: operator.attrgetter("dihedrals"),
        }

    def _measure_iterator(
        self, measurements: Iterable[BaseMeasure]
    ) -> Iterable[tuple[BaseMeasure, set[BaseMeasure]]]:
        for measure in measurements:
            target_getter = self._type_mapping.get(type(measure), None)
            if target_getter is None:
                pass

            yield measure, target_getter(self)

    def update(self, measurements: Iterable[BaseMeasure] | BaseMeasure) -> None:
        """Updates existing stored measurements with `measurements`."""
        if not isinstance(measurements, Iterable):
            measurements = [measurements]

        for measure, existing_measures in self._measure_iterator(measurements):
            if measure in existing_measures:
                existing_measures.pop(measure).update(measure)
            existing_measures.add(measure)

    def remove(self, measurements: Iterable[BaseMeasure] | BaseMeasure) -> None:
        """Removes measurement(s) from stored list, if they exist."""
        if not isinstance(measurements, Iterable):
            measurements = [measurements]

        for measure, existing_measures in self._measure_iterator(measurements):
            existing_measures.pop(measure)
