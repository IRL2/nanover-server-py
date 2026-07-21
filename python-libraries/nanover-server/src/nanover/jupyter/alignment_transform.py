import numpy as np
import numpy.typing as npt
from MDAnalysis import Universe

from .utilities import SceneObjectsUtility
from .frame_listener import FrameListener
from nanover.core import AppServerMinimalImd
from nanover.trajectory import FrameData
from nanover.utilities.transforms import (
    Transform,
    find_transformation_between_point_patterns,
)

from MDAnalysis.lib.transformations import (
    quaternion_from_matrix,
    decompose_matrix,
)


class AlignmentTransform(FrameListener):
    key: str | None = None
    objects: SceneObjectsUtility | None = None

    def __init__(self, app_server: AppServerMinimalImd):
        super().__init__(app_server)

        self.atoms: list[int] = []
        self.positions: npt.NDArray = np.array([])
        self.transform = Transform.identity()

    def config(self, *, key: str, objects: SceneObjectsUtility):
        self.key = key
        self.objects = objects

    def set_atoms_from_positions(self, *, positions: npt.NDArray, atoms: list[int]):
        self.atoms = atoms
        self.positions = positions

    def set_atoms_from_framedata(self, *, frame: FrameData, atoms: list[int]):
        self.atoms = atoms
        self.positions = frame.particle_positions[atoms]

    def set_atoms_from_universe(self, *, universe: Universe, atoms: list[int]):
        self.atoms = atoms
        self.positions = universe.atoms.positions[atoms] / 10  # angstrom -> nm

    def on_frame_update(self, full_frame: FrameData, frame_update: FrameData):
        if self.key is not None:
            self.update_from_framedata(full_frame)

    def update_from_matrix(self, matrix: npt.NDArray):
        self.transform = Transform.from_parent_to_local_matrix(matrix)

        r = quaternion_from_matrix(matrix)
        r = [*r[1:], r[0]]
        s, _, _, t, _ = decompose_matrix(matrix)

        if self.key is not None and self.objects is not None:
            self.objects.update_object(
                f"transform.{self.key}",
                dict(
                    parent="simulation",
                    transform=[*t, *r, *s],
                ),
            )

    def update_from_framedata(self, frame: FrameData):
        self.update_from_positions(frame.particle_positions[self.atoms])

    def update_from_universe(self, universe: Universe):
        self.update_from_positions(universe.atoms.positions[self.atoms])

    def update_from_positions(self, positions: npt.NDArray):
        matrix = find_transformation_between_point_patterns(self.positions, positions)
        self.update_from_matrix(matrix)
