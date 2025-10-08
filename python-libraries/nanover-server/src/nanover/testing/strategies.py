import numpy as np
from hypothesis import strategies as st
from nanover.trajectory import FrameData
import nanover.trajectory.keys as keys
from nanover.trajectory.convert import converters, pack_dict_frame


def uint8s():
    return st.integers(min_value=0, max_value=2**8 - 1)


def uint32s():
    return st.integers(min_value=0, max_value=2**32 - 1)


def float32s():
    return st.floats(width=32, allow_nan=False)


@st.composite
def arrays(draw, values, dtype=None):
    value = draw(st.lists(values))
    return np.array(value, dtype=dtype)


@st.composite
def arrays2d(draw, values, size, dtype=None):
    value = draw(arrays(st.lists(values, min_size=size, max_size=size)))
    return np.array(value, dtype=dtype)


def enum_arrays():
    return arrays(uint8s(), dtype=np.uint8)


def index_arrays():
    return arrays(uint32s(), dtype=np.uint32)


def vec3_arrays():
    return arrays2d(float32s(), dtype=np.float32, size=3)


def index2_arrays():
    return arrays2d(uint32s(), dtype=np.uint32, size=2)


known_types = {
    keys.FRAME_INDEX: st.integers(min_value=1),
    keys.BOX_VECTORS: vec3_arrays(),
    keys.PARTICLE_POSITIONS: vec3_arrays(),
    keys.PARTICLE_ELEMENTS: enum_arrays(),
    keys.PARTICLE_RESIDUES: index_arrays(),
    keys.BOND_PAIRS: index2_arrays(),
    keys.RESIDUE_CHAINS: index_arrays(),
}


@st.composite
def frames(draw):
    return FrameData(draw(dict_frames()))


@st.composite
def packed_frame_dicts(draw):
    frame_dict = draw(dict_frames())
    packed_dict = pack_dict_frame(frame_dict)
    return packed_dict


@st.composite
def dict_frames(draw):
    custom = st.dictionaries(
        keys=custom_frame_keys(),
        values=packable_structures(),
    )

    known = st.fixed_dictionaries(
        mapping={},
        optional=known_types,
    )

    frame = {}
    frame.update(draw(custom))
    frame.update(draw(known))

    return frame


@st.composite
def state_updates(draw):
    arguments = st.dictionaries(
        keys=dictionary_keys().filter(lambda key: not key.startswith("interaction.")),
        values=packable_structures(),
    )

    return draw(arguments)


@st.composite
def command_arguments(draw):
    arguments = st.dictionaries(
        keys=st.text(max_size=32),
        values=packable_structures(),
    )

    return draw(arguments)


@st.composite
def packable_structures(draw):
    structures = st.recursive(
        packable_primitives(),
        lambda children: st.one_of(
            st.lists(children),
            st.dictionaries(dictionary_keys(), values=children),
        ),
        max_leaves=10,
    )

    return draw(structures)


@st.composite
def custom_frame_keys(draw):
    return draw(dictionary_keys().filter(lambda key: key not in converters))


@st.composite
def dictionary_keys(draw):
    return draw(st.text(max_size=16))


@st.composite
def packable_primitives(draw):
    primitive = st.one_of(
        st.none(),
        st.booleans(),
        st.integers(min_value=-(2**63), max_value=2**63),
        st.floats(allow_nan=False),
        st.text(max_size=16),
    )

    return draw(primitive)
