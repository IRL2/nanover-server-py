from hypothesis import given
from nanover.testing.utilities import simplify_numpy
from nanover.trajectory.frame_dict import frame_dict_packer
from nanover.testing.strategies import dict_frames


@given(dict_frames())
def test_pack_unpack_dict_frames(dict_frame):
    """
    Test that frames using both arbitrary and known frame fields can be packed and then unpacked consistently.
    """
    packed = frame_dict_packer.pack(dict_frame)
    unpacked = frame_dict_packer.unpack(packed)

    assert simplify_numpy(unpacked) == simplify_numpy(dict_frame)
