import pytest
import types

from nanover.trajectory.measure import Scalar, Distance, Angle, Dihedral
from nanover.trajectory.measure_collections import MeasureCollection


def _create_unique_measures(class_, num, offset=0):
    _submodule = types.SimpleNamespace()
    _submodule.Scalar, _submodule.Distance, _submodule.Angle, _submodule.Dihedral = (
        Scalar,
        Distance,
        Angle,
        Dihedral,
    )

    match class_:
        case _submodule.Scalar:
            gen = (Scalar(f"test{i + offset}", i) for i in range(num))
        case _submodule.Distance:
            gen = (
                Distance("", *map(lambda x: x + i + offset, range(2)), i)
                for i in range(num)
            )
        case _submodule.Angle:
            gen = (
                Angle("", *map(lambda x: x + i + offset, range(3)), i)
                for i in range(num)
            )
        case _submodule.Dihedral:
            gen = (
                Dihedral("", *map(lambda x: x + i + offset, range(4)), i)
                for i in range(num)
            )

    yield from gen


@pytest.fixture
def create_measure_collection():
    def collection(num_scalars=0, num_distances=0, num_angles=0, num_dihedral=0):
        return MeasureCollection(
            _create_unique_measures(Scalar, num_scalars),
            _create_unique_measures(Distance, num_distances),
            _create_unique_measures(Angle, num_angles),
            _create_unique_measures(Dihedral, num_dihedral),
        )

    return collection


def test_measure_collection_ctr_empty():
    collection = MeasureCollection()

    assert all(
        len(el) == 0
        for el in [
            collection.scalars,
            collection.distances,
            collection.angles,
            collection.dihedrals,
        ]
    )


def test_measure_collection_ctr_single_type():
    scalars = list(_create_unique_measures(Scalar, 3))

    collection = MeasureCollection(scalars=scalars)
    assert len(collection.scalars) == 3
    assert all(
        len(el) == 0
        for el in [collection.distances, collection.angles, collection.dihedrals]
    )

    # Also test invalid measure types
    with pytest.raises(TypeError):
        collection = MeasureCollection(distances=scalars)


def test_measure_collection_ctr_multiple_types():
    scalars = list(_create_unique_measures(Scalar, 4))
    distances = list(_create_unique_measures(Distance, 1))
    angles = list(_create_unique_measures(Angle, 2))
    dihedrals = list(_create_unique_measures(Dihedral, 3))

    collection = MeasureCollection(scalars, distances, angles, dihedrals)
    assert len(collection.scalars) == 4
    assert len(collection.distances) == 1
    assert len(collection.angles) == 2
    assert len(collection.dihedrals) == 3

    # Also check invalid instantiation
    with pytest.raises(TypeError):
        collection = MeasureCollection(scalars, scalars, angles, dihedrals)


def test_measure_collection_item_access(create_measure_collection):
    collection: MeasureCollection = create_measure_collection(3)
    scalar = next(_create_unique_measures(Scalar, 1))
    distance = next(_create_unique_measures(Distance, 1))

    assert len(collection.scalars) == 3
    assert collection[scalar] == scalar
    with pytest.raises(KeyError):
        collection[distance]  # Invalid entry should fail.


def test_measure_collection_update():
    set1 = _create_unique_measures(Scalar, 3)
    set2 = _create_unique_measures(Scalar, 2, offset=5)

    collection = MeasureCollection(set1)
    assert len(collection.scalars) == 3

    # Updating with different measure should not add to first set.
    collection.update(_create_unique_measures(Distance, 2))
    assert len(collection.scalars) == 3
    assert len(collection.distances) == 2

    # Adding more unique scalars should add new keys.
    collection.update(set2)
    assert len(collection.scalars) == 5
    assert len(collection.distances) == 2

    # Adding scalars with the same names should update values but keep the same length.
    different_values = list(_create_unique_measures(Scalar, 2))
    for el in different_values:
        el.update(10)
    collection.update(different_values)
    assert len(collection.scalars) == 5
    assert all(collection[el].value == 10 for el in different_values)


def test_measure_collection_removal(create_measure_collection):
    collection: MeasureCollection = create_measure_collection(1, 2, 3, 4)
    scalar = next(_create_unique_measures(Scalar, 1))

    assert len(collection.scalars) == 1
    assert len(collection.distances) == 2
    collection.remove(scalar)
    assert len(collection.scalars) == 0
    assert len(collection.distances) == 2

    collection.remove(scalar)  # Removal of non-existant key should not throw an error
    collection.remove(
        next(_create_unique_measures(Distance, 1, offset=10))
    )  # Same with non-epty collection.

    # Check removal of iterable is also OK.
    collection.remove(_create_unique_measures(Angle, 3))
    assert len(collection.angles) == 0

    # Removal of iterable with non-entry should also work
    collection.remove(_create_unique_measures(Dihedral, 10))
    assert len(collection.dihedrals) == 0
