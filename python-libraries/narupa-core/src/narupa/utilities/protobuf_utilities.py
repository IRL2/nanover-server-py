from typing import Dict, List, Iterable, Mapping, Union, Any

from google.protobuf.internal.well_known_types import _SetStructValue, _GetStructValue  # type: ignore
from google.protobuf.json_format import MessageToDict
from google.protobuf.struct_pb2 import Struct, Value, ListValue

# Mypy does not support recursive types, yet.
# See <https://github.com/python/mypy/issues/731>.
# This means we cannot write the type in a straight forward way such as:
# Serializable = Union[
#     None, str, int, float, bool,
#     Iterable['Serializable'],
#     Mapping[str, 'Serializable'],
# ]
# Instead, I separate the nested levels in their own types and fall back to Any
# after a couple of levels. This should be OK as we barely use nesting.
_SerializablePrimitive = Union[None, str, int, float, bool]
_TerminalIterable = Iterable[Any]
_TerminalMapping = Mapping[str, Any]
_Level0Iterable = Iterable[Union[_SerializablePrimitive, _TerminalIterable, _TerminalMapping]]
_Level0Mapping = Mapping[str, Union[_SerializablePrimitive, _TerminalIterable, _TerminalMapping]]
Serializable = Union[
    _SerializablePrimitive,
    _Level0Iterable,
    _Level0Mapping,
]


def dict_to_struct(dictionary: Dict[str, Serializable]) -> Struct:
    """
    Converts a python dictionary to a protobuf :class:`Struct`.
    The dictionary must consist of types that can be serialised.

    :param dictionary: Dictionary to convert.
    :return: :class:`Struct` containing copies of all the items of the dictionary.
    """
    struct = Struct()
    try:
        struct.update(dictionary)
    except (ValueError, TypeError, AttributeError):
        raise ValueError(
            'Could not convert object into a protobuf struct. The object to '
            'be converted must be a dictionary with string keys containing '
            'only numbers, strings, booleans and collections of those types. '
        )
    return struct


def struct_to_dict(struct: Struct) -> Dict[str, Serializable]:
    """
    Converts a protobuf :class:`Struct` to a python dictionary.

    :param struct: :class:`Struct` to convert.
    :return: A dictionary containing copies of all the items in the struct.
    """
    return {key: value_to_object(value) for key, value in struct.fields.items()}


def list_value_to_list(list: ListValue) -> List[Serializable]:
    """
    Converts a protobuf :class:`ListValue` to a python dictionary.

    :param list: :class:`ListValue` to convert.
    :return: A list containing copies of all the items in the list value.
    """
    return [value_to_object(value) for value in list.values]


def object_to_value(obj: Serializable) -> Value:
    """
    Convert a python object in an equivalent protobuf Value.
    :param obj: A python object.
    :return: A protobuf Value equivalent to the given object.
    """
    value = Value()
    _SetStructValue(value, obj)
    return value


def value_to_object(value: Value) -> Serializable:
    """
    Converts a protobuf Value into an equivalent python object.
    :param value: A protobuf Value to convert.
    :return: A python object equivalent to the given protobuf Value.
    """
    expanded = _GetStructValue(value)
    if isinstance(expanded, ListValue):
        return list_value_to_list(expanded)
    if isinstance(expanded, Struct):
        return struct_to_dict(expanded)
    return expanded


def deep_copy_serializable_dict(dictionary: Dict[str, Serializable]) -> Dict[str, Serializable]:
    """
    Makes a deep copy of a dictionary by converting it to a protobuf Struct and
    back. Only protobuf serializable elements will be preserved.
    """
    return struct_to_dict(dict_to_struct(dictionary))
