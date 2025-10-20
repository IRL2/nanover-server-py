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

    def update(self, new_measurements: Iterable[BaseMeasure]) -> None:
        """Updates existing stored measurements with new data from `new_measurements`."""

        for measure in new_measurements:
            target_getter = self._type_mapping.get(type(measure), None)
            if target_getter is None:
                pass

            existing_measures = target_getter(self)
            if measure in existing_measures:
                existing_measures.pop(measure).update(measure)
            existing_measures.add(measure)
