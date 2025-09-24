import numpy as np
from hypothesis import strategies as st
from nanover.trajectory.frame_data import (
    PARTICLE_POSITIONS,
    PARTICLE_ELEMENTS,
    PARTICLE_RESIDUES,
    BOND_PAIRS,
    RESIDUE_CHAINS,
    BOX_VECTORS,
)
from nanover.websocket.convert import converters, convert_dict_frame_to_grpc_frame


def uint8s():
    return st.integers(min_value=0, max_value=2**8 - 1)


def uint32s():
    return st.integers(min_value=0, max_value=2**32 - 1)


def float32s():
    return st.floats(width=32, allow_nan=False)


@st.composite
def arrays(draw, values, shape=(-1, 1), dtype=None):
    size = shape[1]
    value = draw(st.lists(st.lists(values, min_size=size, max_size=size)))
    return np.array(value, dtype=dtype)


def np_vec3s():
    return arrays(float32s(), dtype=np.float32, shape=(-1, 3))


known_types = {
    BOX_VECTORS: np_vec3s(),
    PARTICLE_POSITIONS: np_vec3s(),
    PARTICLE_ELEMENTS: arrays(uint8s()),
    PARTICLE_RESIDUES: arrays(uint32s()),
    BOND_PAIRS: arrays(uint32s()),
    RESIDUE_CHAINS: arrays(uint32s()),
}


@st.composite
def grpc_frames(draw):
    known = st.fixed_dictionaries(
        mapping={},
        optional=known_types,
    )

    dict_frame = {}
    dict_frame.update(draw(known))

    grpc_frame = convert_dict_frame_to_grpc_frame(dict_frame)
    return grpc_frame


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
