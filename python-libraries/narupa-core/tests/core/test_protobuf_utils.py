import pytest
from google.protobuf.struct_pb2 import ListValue, Struct

from narupa.core.protobuf_utilities import dict_to_struct


def assert_struct_dictionary_equal(struct, dictionary):
    assert len(dictionary) == len(struct)
    for key in dictionary:
        if isinstance(struct[key], ListValue):
            assert pytest.approx(struct[key]) == dictionary[key]
        elif isinstance(struct[key], Struct):
            assert_struct_dictionary_equal(struct[key], dictionary[key])
        else:
            assert struct[key] == dictionary[key]


@pytest.mark.parametrize('dictionary',
                         [
                             {'int': 1, 'bool': True, 'list': [1.0, 2.4, 3.8], 'str': 'hello'},
                             {},
                             {'recursive': {'int': 5}}
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