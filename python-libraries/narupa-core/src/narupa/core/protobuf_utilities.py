from typing import Dict

from google.protobuf.json_format import MessageToDict
from google.protobuf.struct_pb2 import Struct


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
    except (ValueError, TypeError):
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