from hypothesis import strategies as st
from nanover.trajectory.frame_data import (
    PARTICLE_POSITIONS,
    PARTICLE_ELEMENTS,
    PARTICLE_RESIDUES,
    BOND_PAIRS,
    RESIDUE_CHAINS,
    BOX_VECTORS,
)
from nanover.websocket.convert import converters


def uint8s():
    return st.integers(min_value=0, max_value=2**8 - 1)


def uint32s():
    return st.integers(min_value=0, max_value=2**32 - 1)


def float32s():
    return st.floats(width=32, allow_nan=False)


known_types = {
    BOX_VECTORS: st.lists(float32s()),
    PARTICLE_POSITIONS: st.lists(float32s()),
    PARTICLE_ELEMENTS: st.lists(uint8s()),
    PARTICLE_RESIDUES: st.lists(uint32s()),
    BOND_PAIRS: st.lists(uint32s()),
    RESIDUE_CHAINS: st.lists(uint32s()),
}


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
