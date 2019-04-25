# TODO: tests for containers as values
# TODO: tests for None as value

import sys
import numpy as np

import pytest
from hypothesis import strategies as st
from hypothesis import given

from narupa.protocol.trajectory import FrameData as GrpcFrameData
from narupa.trajectory.frame_data import (
    FrameData, RecordView, PYTHON_TYPES_TO_GRPC_VALUE_ATTRIBUTE,
)
from narupa.trajectory import frame_data

MAX_DOUBLE = sys.float_info.max
MIN_DOUBLE = sys.float_info.min
MAX_FLOAT32 = np.finfo(np.float32).max
MIN_FLOAT32 = np.finfo(np.float32).min
MAX_UINT32 = np.iinfo(np.uint32).max

# This strategy generates a single value (i.e. not a container) that is valid
# in as value in a FrameData, and that can be safely compared with "==" (i.e.
# not numbers as ints are stored as doubles).
EXACT_SINGLE_VALUE_STRATEGY = st.one_of(
    st.text(), st.booleans(),
)

NUMBER_SINGLE_VALUE_STRATEGY = st.one_of(
    st.floats(min_value=float(MIN_DOUBLE), max_value=float(MAX_DOUBLE)),
    st.integers(min_value=MIN_DOUBLE, max_value=MAX_DOUBLE),
)

ARRAYS_STRATEGIES = {
    'string_values': st.lists(st.text(), min_size=1),
    'index_values': st.lists(st.integers(min_value=0, max_value=MAX_UINT32), min_size=1),
    'float_values': st.lists(st.floats(min_value=float(MIN_FLOAT32), max_value=float(MAX_FLOAT32)), min_size=1),
}



@st.composite
def raw_frame_with_single_value(draw, value_strategy):
    """
    Strategy to generate tuple of a GRPC FrameData with a single value
    and that value.

    The value is stored under the key "sample.value".

    The return value is a tuple containing the GRP FrameData, and the
    value by itself.

    The value in the GRPC FrameData is taken from the strategy provided as
    "value_strategy".
    """
    value = draw(value_strategy)
    raw = GrpcFrameData()
    attribute_value = PYTHON_TYPES_TO_GRPC_VALUE_ATTRIBUTE[type(value)]
    setattr(raw.values['sample.value'], attribute_value, value)
    return raw, value


@st.composite
def raw_frame_with_single_array(draw, value_strategy, field):
    """
    Strategy to generate tuple of a GRPC FrameData with a single array
    and the corresponding list.

    The array is stored under the key "sample.array".

    The return value is a tuple containing the GRP FrameData, and the
    list corresponding to the stored array.

    The array in the GRPC FrameData is taken from the strategy provided as
    "value_strategy". The "field" argument tells under which field from the
    ValueArray object the array must be stored.
    """
    value = draw(value_strategy)
    raw = GrpcFrameData()
    getattr(raw.arrays['sample.array'], field).values.extend(value)
    return raw, value


@pytest.fixture
def simple_frame():
    """
    Create a FrameData with some filled fields.
    """
    raw = GrpcFrameData()
    raw.values['sample.number'].number_value = 4.6
    raw.values['sample.string'].string_value = 'foo bar'
    raw.values['sample.true'].bool_value = True
    raw.values['sample.false'].bool_value = False
    raw.arrays['array.index'].index_values.values.extend(range(18, 0, -3))
    raw.arrays['array.float'].float_values.values.extend([2.3, 4.5, 6.7])
    raw.arrays['array.string'].string_values.values.extend(['foo', 'bar', 'toto'])

    raw.arrays[frame_data.POSITIONS].float_values.values.extend((
        1.0, 2.1, 3.2,
        4.3, 5.4, 6.5,
        7.6, 8.7, 9.8,
        0.9, 1.1, 2.2,
    ))
    raw.arrays[frame_data.BONDS].index_values.values.extend((
        0, 1,
        1, 2,
        2, 3,
    ))
    raw.arrays[frame_data.ELEMENTS].index_values.values.extend((10, 12, 14, 16))
    return FrameData(raw)


@given(raw_frame_with_single_value(st.floats()))
def test_get_single_numeric_value(raw_frame_and_expectation):
    """
    We can get a numeric value by key from a FrameData.

    Numeric values are stored as doubles by protobuf. We therefore cannot
    compare their exact values.
    """
    raw_frame, expected_value = raw_frame_and_expectation
    frame = FrameData(raw_frame)
    assert pytest.approx(frame.values['sample.value'], expected_value)


@given(raw_frame_with_single_value(EXACT_SINGLE_VALUE_STRATEGY))
def test_get_single_exact_value(raw_frame_and_expectation):
    """
    We can get a value by key from a FrameData.

    This test covers all non-container types that can be compared safely with
    an equality check (i.e. not numeric).
    """
    raw_frame, expected_value = raw_frame_and_expectation
    frame = FrameData(raw_frame)
    assert frame.values['sample.value'] == expected_value


# Arguments for the @given decorator cannot be parametrized. So we wrote a
# separate test for each strategy. This function contains the logic of these
# tests.
def _test_get_exact_array(raw_frame_and_expectation):
    raw_frame, expected_value = raw_frame_and_expectation
    frame = FrameData(raw_frame)
    assert frame.arrays['sample.array'] == expected_value


@given(raw_frame_with_single_array(ARRAYS_STRATEGIES['string_values'], 'string_values'))
def test_get_string_array(raw_frame_and_expectation):
    """
    We can get an array of strings from a FrameData.
    """
    _test_get_exact_array(raw_frame_and_expectation)


@given(raw_frame_with_single_array(ARRAYS_STRATEGIES['index_values'], 'index_values'))
def test_get_index_array(raw_frame_and_expectation):
    """
    We can get an array of indices from a FrameData.
    """
    _test_get_exact_array(raw_frame_and_expectation)


@given(raw_frame_with_single_array(ARRAYS_STRATEGIES['float_values'], 'float_values'))
def test_get_float_array(raw_frame_and_expectation):
    """
    We can get an array of floats from a FrameData.
    """
    raw_frame, expected_value = raw_frame_and_expectation
    frame = FrameData(raw_frame)
    assert all(
        pytest.approx(actual, expected)
        for actual, expected in zip(frame.arrays['sample.array'], expected_value)
    )


@pytest.mark.parametrize('record_name', ('values', 'arrays'))
def test_missing_value_key(record_name):
    """
    We get a KeyError when requesting a key that does not exist.
    """
    frame = FrameData()
    with pytest.raises(KeyError):
        getattr(frame, record_name)['missing.key']


@pytest.mark.parametrize('key, record_type, expected', (
    ('sample.number', 'values', True),
    ('sample.string', 'values', True),
    ('sample.true', 'values', True),
    ('sample.false', 'values', True),
    ('array.index', 'arrays', True),
    ('array.float', 'arrays', True),
    ('array.string', 'arrays', True),
    ('missing.key', 'values', False),
    ('missing.key', 'arrays', False),
))
def test_record_view_contains(simple_frame, key, record_type, expected):
    """
    The "in" operator returns the expected value on ValuesView and ArraysView.
    """
    assert (key in getattr(simple_frame, record_type)) == expected


@pytest.mark.parametrize('key, expected', (
    ('sample.number', True),
    ('sample.string', True),
    ('sample.true', True),
    ('sample.false', True),
    ('array.index', True),
    ('array.float', True),
    ('array.string', True),
    ('missing.key', False),
))
def test_frame_data_contains(simple_frame, key, expected):
    """
    The "in" operator returns the expected value on FrameData.
    """
    assert (key in simple_frame) == expected


@pytest.mark.parametrize('key, record_type, expected', (
    ('sample.string', 'values', 'foo bar'),
    ('array.string', 'arrays', ['foo', 'bar', 'toto']),
    ('missing', 'values', 'default'),  # Not in the FrameData
    ('missing', 'arrays', 'default'),  # Not in the FrameData
))
def test_RecordView_get(simple_frame, key, record_type, expected):
    """
    The "get" method of ValuesView and ArraysView return the expected value.
    """
    assert getattr(simple_frame, record_type).get(key, default='default') == expected


def test_base_class_init_fails():
    """
    Instantiating a non-subclassed RecordView fails with a NotImplemented Error.
    """
    with pytest.raises(NotImplementedError):
        RecordView(GrpcFrameData())


def test_partial_view_fails_singular():
    """
    If a subclass of RecordView is missing the "singular" attribute, it fails.
    """
    class DummyView(RecordView):
        record_name = 'values'

    with pytest.raises(NotImplementedError):
        DummyView(GrpcFrameData())


def test_partial_view_fails_record_name():
    """
    If a subclass of RecordView is missing the "record_type" attribute, it fails.
    """
    class DummyView(RecordView):
        singular = 'dummy'

    with pytest.raises(NotImplementedError):
        DummyView(GrpcFrameData())


def test_partial_view_fails_converter():
    """
    If a subclass of RecordView is missing the '_convert_to_python' method, it fails.
    """
    class DummyView(RecordView):
        singular = 'dummy'
        record_name = 'values'

    raw = GrpcFrameData()
    raw.values['sample'].string_value = 'foobar'

    dummy = DummyView(raw)
    with pytest.raises(NotImplementedError):
        dummy.get('sample')


def test_partial_view_fails_set():
    """
    If a subclass of RecordView is missing the 'set' method, it fails.
    """
    class DummyView(RecordView):
        singular = 'dummy'
        record_name = 'values'

    raw = GrpcFrameData()
    raw.values['sample'].string_value = 'foobar'

    dummy = DummyView(raw)
    with pytest.raises(NotImplementedError):
        dummy.set('sample', 0)


def test_positions_shortcut_get(simple_frame):
    """
    Test that the "positions" shortcut of FrameData returns the expected value.

    Because "positions" contains floats, we need to test it separately from the
    other shortcuts that can be compared exactly.
    """
    positions = simple_frame.positions
    expected = [
        [1.0, 2.1, 3.2],
        [4.3, 5.4, 6.5],
        [7.6, 8.7, 9.8],
        [0.9, 1.1, 2.2],
    ]
    # pytest.approx can compare lists, but it does not support nested structures.
    # Here we compare the 2D lists row per row.
    assert all(
        pytest.approx(positions_row, expected_row)
        for positions_row, expected_row in zip(positions, expected)
    )


@pytest.mark.parametrize('key, expected', (
    ('bonds', [[0, 1], [1, 2], [2, 3]]),
    ('elements', [10, 12, 14, 16]),
))
def test_exact_shortcuts_get(simple_frame, key, expected):
    """
    The shortcuts that return exact values return the expected one.
    """
    assert getattr(simple_frame, key) == expected


@pytest.mark.parametrize('value, expected', (
    (np.array([[3.4, 2.5, 1.2], [5.6, 2.1, 6.7]]), [3.4, 2.5, 1.2, 5.6, 2.1, 6.7]),
    ([[3.4, 2.5, 1.2], [5.6, 2.1, 6.7]], [3.4, 2.5, 1.2, 5.6, 2.1, 6.7]),
))
def test_positions_shortcuts_set(value, expected):
    """
    The "positions" shortcut can be set from a list or from a numpy array.
    """
    frame = FrameData()
    frame.positions = value
    assert pytest.approx(
        frame.raw.arrays[frame_data.POSITIONS].float_values.values, expected
    )


@pytest.mark.parametrize('key, raw_key, value, expected', (
    ('bonds', frame_data.BONDS, [[3, 4], [2, 3]], [3, 4, 2, 3]),
    ('elements', frame_data.ELEMENTS, [2, 3, 5], [2, 3, 5]),
))
def test_exact_shortcuts_set(key, raw_key, value, expected):
    """
    The shortcuts with exact values can be set.

    This test can only cover the shortcut that set "index_values".
    """
    frame = FrameData()
    setattr(frame, key, value)
    assert frame.raw.arrays[raw_key].index_values.values == expected


@pytest.mark.parametrize('exception', (KeyError, frame_data.MissingDataError))
def test_missing_shortcut(exception):
    """
    If there is no data for a given shortcut, it can be caught by looking for
    a KeyError or a ad-hoc MissingDataError.
    """
    frame = FrameData()
    with pytest.raises(exception):
        frame.positions


@given(EXACT_SINGLE_VALUE_STRATEGY)
def test_set_new_exact_value(value):
    """
    Values can be set. Assumes exact comparison.
    """
    frame = FrameData()
    frame.values['sample.new'] = value
    type_attribute = PYTHON_TYPES_TO_GRPC_VALUE_ATTRIBUTE[type(value)]
    assert 'sample.new' in frame.raw.values
    assert getattr(frame.raw.values['sample.new'], type_attribute) == value


@given(NUMBER_SINGLE_VALUE_STRATEGY)
def test_set_new_number_value(value):
    """
    Numbers can be set as values.
    """
    frame = FrameData()
    frame.values['sample.new'] = value
    type_attribute = PYTHON_TYPES_TO_GRPC_VALUE_ATTRIBUTE[type(value)]
    assert 'sample.new' in frame.raw.values
    assert pytest.approx(getattr(frame.raw.values['sample.new'], type_attribute), value)


@given(
    initial_value=NUMBER_SINGLE_VALUE_STRATEGY,
    new_value=EXACT_SINGLE_VALUE_STRATEGY,
)
def test_set_existing_value_number_to_exact(initial_value, new_value):
    """
    An existing value can be modified by a value of a different type.
    """
    frame = FrameData()
    frame.values['sample.new'] = initial_value
    assert pytest.approx(frame.values['sample.new'], initial_value)
    frame.values['sample.new'] = new_value
    assert frame.values['sample.new'] == new_value


@given(
    initial_value=EXACT_SINGLE_VALUE_STRATEGY,
    new_value=NUMBER_SINGLE_VALUE_STRATEGY,
)
def test_set_existing_value_exact_to_number(initial_value, new_value):
    """
    An existing value can be modified by a value of a different type.

    This test exchanges the order of the types from the initial value ato the
    new value compared to `test_set_existing_value_number_to_exact`. In case
    the order of the fields in the raw FrameData have an influence.
    """
    frame = FrameData()
    frame.values['sample.new'] = initial_value
    assert frame.values['sample.new'] == initial_value
    frame.values['sample.new'] = new_value
    assert pytest.approx(frame.values['sample.new'], new_value)


@given(st.one_of(
    ARRAYS_STRATEGIES['index_values'],
    ARRAYS_STRATEGIES['string_values'],
))
def test_set_new_exact_array(value):
    """
    A new array of exact values can be set.
    """
    frame = FrameData()
    frame.arrays['sample.new'] = value
    assert frame.arrays['sample.new'] == value


@given(ARRAYS_STRATEGIES['float_values'])
def test_set_new_float_array(value):
    """
    A new array of numbers can be set.
    """
    frame = FrameData()
    frame.arrays['sample.new'] = value
    assert pytest.approx(frame.arrays['sample.new'], value)


@given(
    init_value=ARRAYS_STRATEGIES['index_values'],
    new_value=ARRAYS_STRATEGIES['string_values'],
)
def test_set_existing_array(init_value, new_value):
    """
    An array can be replaced by an array of a different type.
    """
    frame = FrameData()
    frame.arrays['sample.new'] = init_value
    assert frame.arrays['sample.new'] == init_value
    frame.arrays['sample.new'] = new_value
    assert frame.arrays['sample.new'] == new_value


@pytest.mark.parametrize('value', (
    [],  # Empty list, type cannot be guessed
    21,  # Not an sequence
    [None, None],  # Not a valid type
))
def test_set_wrong_type_array_fails(value):
    """
    Setting an array with a broken value fails with a ValueError.
    """
    frame = FrameData()
    with pytest.raises(ValueError):
        frame.arrays['sample.net'] = value
