from collections.abc import Iterable, Sequence

import numpy as np
import numpy.typing as npt
from MDAnalysis import AtomGroup, Universe
from MDAnalysis.lib import transformations
from nanover.trajectory import FrameData

# position[x,y,z], rotation[x,y,z,w], scale[x,y,z]
STATE_TRANSFORM_IDENTITY = (0, 0, 0, 0, 0, 0, 1, 1, 1, 1)


def unpack_partial_state_transform(transform: Iterable[float]):
    """Pad partial state transform iterable up to full length with components of identity state transform."""
    i = -1
    for i, value in enumerate(transform):
        yield value
    for j in range(i + 1, len(STATE_TRANSFORM_IDENTITY)):
        yield STATE_TRANSFORM_IDENTITY[j]


def matrix_from_state_transform(transform: Iterable[float]):
    tx, ty, tz, rx, ry, rz, rw, sx, sy, sz = unpack_partial_state_transform(transform)

    translation = transformations.translation_matrix((tx, ty, tz))
    rotation = transformations.quaternion_matrix((rw, rx, ry, rz))
    scale = np.diagflat((sx, sy, sz, 1.0))

    return translation @ rotation @ scale


def state_transform_from_matrix(matrix: npt.NDArray) -> Sequence[float]:
    s, _, a, t, _ = transformations.decompose_matrix(matrix)
    w, x, y, z = transformations.quaternion_from_euler(*a)
    return *t, x, y, z, w, *s


class Transform:
    """Convenience wrapper around transformation matrix."""

    @classmethod
    def identity(cls):
        return cls.from_local_to_parent_matrix(np.identity(4))

    @classmethod
    def from_state_transform(cls, transform: Sequence[float]):
        return cls.from_local_to_parent_matrix(matrix_from_state_transform(transform))

    @classmethod
    def from_state_cursor(cls, cursor):
        return cls.from_state_transform((*cursor["position"], *cursor["rotation"]))

    @classmethod
    def from_local_to_parent_matrix(cls, local_to_parent: npt.NDArray):
        return cls(
            local_to_parent=local_to_parent,
            parent_to_local=np.linalg.inv(local_to_parent),
        )

    @classmethod
    def from_parent_to_local_matrix(cls, parent_to_local: npt.NDArray):
        return cls(
            parent_to_local=parent_to_local,
            local_to_parent=np.linalg.inv(parent_to_local),
        )

    @property
    def local_to_parent_matrix(self):
        return self._local_to_parent

    @property
    def parent_to_local_matrix(self):
        return self._parent_to_local

    def __init__(
        self,
        *,
        local_to_parent: npt.NDArray,
        parent_to_local: npt.NDArray,
    ):
        self._local_to_parent = local_to_parent
        self._parent_to_local = parent_to_local

    def __repr__(self):
        s, _, a, t, _ = transformations.decompose_matrix(self.local_to_parent_matrix)
        parts = []

        if not np.allclose(t, [0, 0, 0]):
            parts.append(f"translate[{t[0]:.1f}, {t[1]:.1f}, {t[2]:.1f}]")
        if not np.allclose(a, [0, 0, 0]):
            a = np.rad2deg(a)
            parts.append(f"rotate[{a[0]:.1f}deg, {a[1]:.1f}deg, {a[2]:.1f}deg]")
        if not np.allclose(s, [1, 1, 1]):
            parts.append(f"scale[{s[0]:.1f}, {s[1]:.1f}, {s[2]:.1f}]")

        return f"<Transform of {' '.join(parts)}>"

    def to_state_transform(self) -> Sequence[float]:
        return state_transform_from_matrix(self.local_to_parent_matrix)

    def point_local_to_parent(self, point):
        return _transform_vec3(self._local_to_parent, point)

    def points_local_to_parent(self, points):
        return _transform_vec3s(self._local_to_parent, points)

    def point_parent_to_local(self, point):
        return _transform_vec3(self._parent_to_local, point)

    def points_parent_to_local(self, points):
        return _transform_vec3s(self._parent_to_local, points)


def _transform_vec3(matrix: npt.NDArray, vector) -> npt.NDArray:
    return (matrix @ np.array([*vector[:3], 1]).reshape(4, 1)).reshape(4)[:3]


def _transform_vec3s(matrix: npt.NDArray, vectors) -> npt.NDArray:
    v = np.asarray(vectors).T
    expanded = np.vstack((v, np.ones([1, v.shape[1]], v.dtype)))
    return (matrix @ expanded)[:-1].T


class StructureAlignment:
    @classmethod
    def from_atom_group(cls, atom_group: AtomGroup):
        return cls(
            atoms=atom_group.indices,
            positions=atom_group.positions / 10,  # angstrom -> nm
        )

    @classmethod
    def from_atoms_and_framedata(cls, atoms: list[int], framedata: FrameData):
        return cls(
            atoms=atoms,
            positions=framedata.particle_positions[atoms],
        )

    def __init__(self, *, atoms: list[int], positions: npt.NDArray):
        self.atoms = atoms
        self.positions = positions

    def calculate_transform_to_universe(self, universe: Universe):
        return self.calculate_transform_to_positions(
            universe.atoms.positions[self.atoms] / 10  # angstrom -> nm
        )

    def calculate_transform_to_framedata(self, framedata: FrameData):
        return self.calculate_transform_to_positions(
            framedata.particle_positions[self.atoms]
        )

    def calculate_transform_to_positions(self, positions: npt.NDArray):
        return Transform.from_local_to_parent_matrix(
            find_transformation_between_point_patterns(self.positions, positions)
        )


def find_transformation_between_point_patterns(
    start_points: npt.NDArray,
    final_points: npt.NDArray,
):
    """
    Return the 4x4 transformation matrix that best maps the start points to the final points.
    """
    # https://web.stanford.edu/class/cs273/refs/umeyama.pdf
    # https://doi.org/10.1109/34.88573
    p = start_points.T
    q = final_points.T

    # centroids
    cen_P = np.mean(p, axis=1).reshape(-1, 1)
    cen_Q = np.mean(q, axis=1).reshape(-1, 1)

    # centered vectors
    X = p - cen_P
    Y = q - cen_Q

    # svd
    U, _sigma, Vt = np.linalg.svd(X @ Y.T)

    # reflection correction
    d = np.identity(Vt.T.shape[1])
    d[-1, -1] = np.linalg.det(Vt.T @ U.T)

    # rotation and translation
    R = Vt.T @ d @ U.T
    t = cen_Q - R @ cen_P

    # convert to TRS matrix4x4
    transform = np.identity(4)
    transform[:3, 3] = t.reshape(-1)
    transform[:3, :3] = R

    return transform
