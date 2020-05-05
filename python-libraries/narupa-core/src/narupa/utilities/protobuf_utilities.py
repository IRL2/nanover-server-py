from typing import Dict

from google.protobuf.internal.well_known_types import _SetStructValue
from google.protobuf.json_format import MessageToDict
from google.protobuf.struct_pb2 import Struct, Value


def dict_to_struct(dictionary: Dict[str, object]) -> Struct:
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
        raise ValueError("Unable to construct serialise dictionary into a protobuf struct. "
                         "Only value types such as numbers, strings, booleans, and collections of those types"
                         "can be serialised.")
    return struct


def struct_to_dict(struct: Struct) -> Dict[str, object]:
    """
    Converts a protobuf :class:`Struct` to a python dictionary.

    :param struct: :class:`Struct` to convert.
    :return: A dictionary containing copies of all the items in the struct.
    """
    return MessageToDict(struct)


def object_to_value(obj: object) -> Value:
    """
    Convert a python object in an equivalent protobuf Value.
    :param obj: A python object.
    :return: A protobuf Value equivalent to the given object.
    """
    value = Value()
    _SetStructValue(value, obj)
    return value


def value_to_object(value: Value) -> object:
    """
    Converts a protobuf Value into an equivalent python object.
    :param value: A protobuf Value to convert.
    :return: A python object equivalent to the given protobuf Value.
    """
    return MessageToDict(value)


def deep_copy_serializable_dict(dictionary: Dict[str, object]) -> Dict[str, object]:
    """
    Makes a deep copy of a dictionary by converting it to a protobuf Struct and
    back. Only protobuf serializable elements will be preserved.
    """
    return struct_to_dict(dict_to_struct(dictionary))
