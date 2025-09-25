from hypothesis import given

from nanover.testing.utilities import simplify_numpy
from nanover.trajectory.convert import (
    pack_dict_frame,
    unpack_dict_frame,
)
from nanover.testing.strategies import dict_frames


@given(dict_frames())
def test_pack_unpack_dict_frames(dict_frame):
    """
    Test that frames using both arbitrary and known frame fields can be packed and then unpacked consistently.
    """
    packed = pack_dict_frame(dict_frame)
    unpacked = unpack_dict_frame(packed)

    simplify_numpy(unpacked)

    assert unpacked == dict_frame
