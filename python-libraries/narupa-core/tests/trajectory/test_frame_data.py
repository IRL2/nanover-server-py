# TODO: tests for containers as values
# TODO: tests for None as value

import sys
import numpy as np

import pytest
from hypothesis import strategies as st
from hypothesis import given

from narupa.protocol.trajectory import FrameData as GrpcFrameData
from narupa.trajectory.frame_data import FrameData

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
    st.floats(min_value=MIN_DOUBLE, max_value=MAX_DOUBLE),
    st.integers(min_value=MIN_DOUBLE, max_value=MAX_DOUBLE),
)

ARRAYS_STRATEGIES = {
    'string_values': st.lists(st.text()),
    'index_values': st.lists(st.integers(min_value=0, max_value=MAX_UINT32)),
    'float_values': st.lists(st.floats(min_value=MIN_FLOAT32, max_value=MAX_FLOAT32)),
}

PYTHON_TYPES_TO_GRPC_VALUE_ATTRIBUTE = {
    int: 'number_value', float: 'number_value', str: 'string_value',
    bool: 'bool_value',
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
