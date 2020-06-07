import numbers
from math import inf, nan

import pytest
from google.protobuf.struct_pb2 import Value
from narupa.utilities.protobuf_utilities import dict_to_struct, object_to_value, value_to_object
from hypothesis import strategies as st, given


def assert_value_equal(proto_value, python_value):
    """
    Checks if a protobuf value and a python value are approximately equal for numeric types,
    or exactly equal otherwise.
    """
    if isinstance(python_value, numbers.Number):
        assert pytest.approx(proto_value, python_value)
    else:
        assert proto_value == python_value


def assert_list_value_equal(list_value, python_list):
    """
    Recursively checks that a protobuf list value, which may consist of values, lists or structs
    is equal to a python list consisting of values, lists or dictionaries.
    """
    assert len(list_value) == len(python_list)
    proto_val: Value
    for proto_val, python_val in zip(list_value, python_list):
        if isinstance(python_val, list):
            assert_struct_dictionary_equal(proto_val, python_val)
        elif isinstance(python_val, dict):
            assert_list_value_equal(proto_val, python_val)
        else:
            assert_value_equal(proto_val, python_val)


def assert_struct_dictionary_equal(struct, dictionary):
    """
    Recursively checks that a protobuf struct value, which may consist of values, lists or structs
    is equal to a python dictionary consisting of values, lists or dictionaries.
    """
    assert len(dictionary) == len(struct)
    for key in dictionary:
        if isinstance(dictionary[key], list):
            assert_list_value_equal(struct[key], dictionary[key])
        elif isinstance(dictionary[key], dict):
            assert_struct_dictionary_equal(struct[key], dictionary[key])
        else:
            assert_value_equal(struct[key], dictionary[key])


@pytest.mark.parametrize('dictionary',
                         [
                             {'int': 1, 'bool': True, 'list': [1.0, 2.4, 3.8], 'str': 'hello'},
                             {'list': [True, False, True]},
                             {'list': ['a', 'b', 'c']},
                             {'float': 4.3},
                             {'small': 1e-12},
                             {'inf': inf},
                             {'nan': nan},
                             {},
                             {'recursive': {'int': 5}},
                             {'recursive_list': {'list': [5, 3, 2]}},
                             {'different_list_types': [5, True, 'hello']},
                             {'struct_list': [{'test': 1}, {'test': 2}, {'test': 3}]},
                         ])
def test_dict_to_struct(dictionary):
    struct = dict_to_struct(dictionary)
    assert_struct_dictionary_equal(struct, dictionary)


@pytest.mark.parametrize('dictionary',
                         [
                             {object(): 'test'},
                             {'test': object()},
                         ])
def test_dict_to_struct_invalid(dictionary):
    with pytest.raises(ValueError):
        _ = dict_to_struct(dictionary)


@st.composite
def grpc_value(draw, value_strategy):
    value = draw(value_strategy)
    return object_to_value(value), value


MIXED_LIST_STRATEGIES = st.lists(st.one_of(st.text(), st.booleans(), st.none()), min_size=1)
MIXED_DICT_STRATEGIES = st.dictionaries(st.text(), st.one_of(st.text(), st.booleans(), st.none()), min_size=1)


@given(grpc_value(st.floats()))
def test_to_and_from_value_float(protobuf_and_expected_value):
    protobuf_value, expected_value = protobuf_and_expected_value
    assert pytest.approx(value_to_object(protobuf_value), expected_value)


@given(grpc_value(st.integers()))
def test_to_and_from_value_integer(protobuf_and_expected_value):
    protobuf_value, expected_value = protobuf_and_expected_value
    assert pytest.approx(value_to_object(protobuf_value), expected_value)


@given(grpc_value(st.text()))
def test_to_and_from_value_string(protobuf_and_expected_value):
    protobuf_value, expected_value = protobuf_and_expected_value
    assert value_to_object(protobuf_value) == expected_value


@given(grpc_value(st.booleans()))
def test_to_and_from_value_bool(protobuf_and_expected_value):
    protobuf_value, expected_value = protobuf_and_expected_value
    assert value_to_object(protobuf_value) == expected_value


@given(grpc_value(MIXED_LIST_STRATEGIES))
def test_to_and_from_value_mixed_lists(protobuf_and_expected_value):
    protobuf_value, expected_value = protobuf_and_expected_value
    assert value_to_object(protobuf_value) == expected_value


@given(grpc_value(MIXED_DICT_STRATEGIES))
def test_to_and_from_value_mixed_dicts(protobuf_and_expected_value):
    protobuf_value, expected_value = protobuf_and_expected_value
    assert value_to_object(protobuf_value) == expected_value
