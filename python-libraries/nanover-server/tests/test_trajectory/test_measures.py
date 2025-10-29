import pytest
from hypothesis import given, example, strategies as st

from openmm.unit import angstrom

from nanover.trajectory import measure

# TODO
# test ctr measures
# test angle conversions


@st.composite
def scalar_measure(draw, name="test_scalar", value=None):
    value = draw(st.floats()) if value is None else value
    return measure.Scalar(name, value)


@st.composite
def distance_measure(draw, name="test_distance", atom1=None, atom2=None, value=None):
    a1, a2 = draw(st.lists(st.integers(0), min_size=2, max_size=2, unique=True))

    value = draw(st.floats()) if value is None else value
    atom1 = a1 if atom1 is None else atom1
    atom2 = a2 if atom2 is None else atom2
    return measure.Distance(name, atom1, atom2, value, angstrom)


@st.composite
def angle_measure(
    draw, name="test_angle", atom1=None, atom2=None, atom3=None, value=None
):
    a1, a2, a3 = draw(st.lists(st.integers(0), min_size=3, max_size=3, unique=True))

    value = draw(st.floats()) if value is None else value
    atom1 = a1 if atom1 is None else atom1
    atom2 = a2 if atom2 is None else atom2
    atom3 = a3 if atom3 is None else atom3
    return measure.Angle(name, atom1, atom2, atom3, value)


@st.composite
def dihedral_measure(
    draw,
    name="test_dihedral",
    atom1=None,
    atom2=None,
    atom3=None,
    atom4=None,
    value=None,
):
    a1, a2, a3, a4 = draw(st.lists(st.integers(0), min_size=4, max_size=4, unique=True))

    value = draw(st.floats()) if value is None else value
    atom1 = a1 if atom1 is None else atom1
    atom2 = a2 if atom2 is None else atom2
    atom3 = a3 if atom3 is None else atom3
    atom4 = a4 if atom4 is None else atom4
    return measure.Dihedral(name, atom1, atom2, atom3, atom4, value)


def _comparison_set(comparible1, comparible2, value1, value2):
    assert (comparible1 == comparible2) == (value1 == value2)
    assert (comparible1 > comparible2) == (value1 > value2)
    assert (comparible1 < comparible2) == (value1 < value2)
    assert (comparible1 >= comparible2) == (value1 >= value2)
    assert (comparible1 <= comparible2) == (value1 <= value2)


def _explicit_comparison_setup(class_, num_atoms):
    atoms = range(num_atoms)

    # Check numbers
    assert class_("", *atoms, 42.0) == 42.0
    assert class_("", *atoms, 42.0) != 0

    # Check against like class - indices
    assert class_("", *atoms, 42.0) == class_("", *atoms, 42.0)
    assert class_("", *atoms, 42.0) != class_("", *reversed(atoms), 42.0)
    assert class_("", *atoms, 42.0) != class_("", *map(lambda x: x + 1, atoms), 42.0)

    # Check against like class - values
    assert class_("", *atoms, 42.0) == class_("", *atoms, 42.0)
    assert class_("", *atoms, 42.0) != class_("", *atoms, 0)

    # Check name is ignored, rather only atom indices matter
    assert class_("test", *atoms, 42.0) == class_("different", *atoms, 42.0)


def test_scalar_explicit_examples():
    # Check numbers
    assert measure.Scalar("", 42.0) == 42.0
    assert measure.Scalar("", 42.0) != 0

    # Check against like class
    assert measure.Scalar("", 42.0) == measure.Scalar("", 42.0)
    assert measure.Scalar("", 42.0) != measure.Scalar("", 0)

    # Check name also fails
    assert measure.Scalar("test", 42.0) != measure.Scalar("different", 42.0)


@given(scalar_measure(), st.floats())
def test_scalar_comparisons_number(scalar, value):
    _comparison_set(scalar, value, scalar.value, value)


@given(scalar_measure(), scalar_measure())
def test_scalar_equality_scalar(first_measure, second_measure):
    _comparison_set(
        first_measure, second_measure, first_measure.value, second_measure.value
    )


@pytest.mark.parametrize(
    "class_,num_atoms",
    [(measure.Distance, 2), (measure.Angle, 3), (measure.Dihedral, 4)],
)
def test_distance_explicit_examples(class_, num_atoms):
    _explicit_comparison_setup(class_, num_atoms)


@given(distance_measure(), st.floats())
def test_distance_equality_number(distance, value):
    _comparison_set(distance, value, distance.value, value)


@given(distance_measure(), distance_measure())
def test_distance_equality_distance(distance, value):
    if (distance.atom1, distance.atom2) == (value.atom1, value.atom2):
        _comparison_set(distance, value, distance.value, value.value)
    else:
        assert distance != value


@given(angle_measure(), st.floats())
def test_angle_equality_number(angle, value):
    _comparison_set(angle, value, angle.value, value)


@given(angle_measure(), angle_measure())
def test_angle_equality_angle(angle, value):
    if (angle.atom1, angle.atom2, angle.atom3) == (
        value.atom1,
        value.atom2,
        value.atom3,
    ):
        _comparison_set(angle, value, angle.value, value.value)
    else:
        assert angle != value


@given(dihedral_measure(), st.floats())
def test_dihedral_equality_number(dihedral, value):
    _comparison_set(dihedral, value, dihedral.value, value)


@given(dihedral_measure(), dihedral_measure())
def test_dihedral_equality_dihedral(dihedral, value):
    if (dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4) == (
        value.atom1,
        value.atom2,
        value.atom3,
        value.atom4,
    ):
        _comparison_set(dihedral, value, dihedral.value, value.value)
    else:
        assert dihedral != value
