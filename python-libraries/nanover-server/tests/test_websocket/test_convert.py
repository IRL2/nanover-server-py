from hypothesis import given, strategies as st
from nanover.websocket.convert import (
    pack_dict_frame,
    unpack_dict_frame,
)
from data_strategies import dict_frames


@given(dict_frames())
def test_pack_unpack_dict_frames(dict_frame):
    """
    Test that frames using both arbitrary and known frame fields can be packed and then unpacked consistently.
    """
    packed = pack_dict_frame(dict_frame)
    unpacked = unpack_dict_frame(packed)
    assert unpacked == dict_frame
