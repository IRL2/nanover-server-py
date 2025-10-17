from typing import Iterable

from nanover.trajectory.measure import Measure, Distance, Angle, Dihedral


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

    def update(self) -> None: ...
