import sys

import pytest
from hypothesis import strategies as st
from hypothesis import given

from narupa.protocol.trajectory import FrameData as GrpcFrameData
from narupa.trajectory.frame_data import FrameData

MAX_DOUBLE = sys.float_info.max
MIN_DOUBLE = sys.float_info.min

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

PYTHON_TYPES_TO_GRPC_ATTRIBUTE = {
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
    attribute_value = PYTHON_TYPES_TO_GRPC_ATTRIBUTE[type(value)]
    setattr(raw.values['sample.value'], attribute_value, value)
    return raw, value


@given(raw_frame_with_single_value(st.floats()))
def test_get_single_float_value(raw_frame_and_expectation):
    """
    We can get a float value by key from a FrameData.
    """
    raw_frame, expected_value = raw_frame_and_expectation
    frame = FrameData(raw_frame)
    assert pytest.approx(frame['sample.value'], expected_value)


@given(raw_frame_with_single_value(EXACT_SINGLE_VALUE_STRATEGY))
def test_get_single_exact_value(raw_frame_and_expectation):
    """
    We can get a value by key from a FrameData.

    This test covers all non-container types that can be compared safely with
    an equality check (i.e. not floats).
    """
    raw_frame, expected_value = raw_frame_and_expectation
    frame = FrameData(raw_frame)
    assert frame['sample.value'] == expected_value


def test_missing_key():
    """
    We get a KeyError when requesting a key that does not exist.
    """
    frame = FrameData()
    with pytest.raises(KeyError):
        frame['missing.key']
