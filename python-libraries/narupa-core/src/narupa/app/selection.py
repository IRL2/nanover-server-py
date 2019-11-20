# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from typing import Dict, Iterable, Set

INTERACTION_SINGLE = 'single'
INTERACTION_GROUP = 'group'
INTERACTION_RESTRAINT = 'restraint'

RENDERER_LIQUORICE = 'liquorice'

KEY_SELECTION_ID = 'id'
KEY_SELECTION_NAME = 'name'
KEY_SELECTION_SELECTED = 'selected'
KEY_SELECTED_PARTICLE_IDS = 'particle_ids'
KEY_SELECTION_PROPERTIES = 'properties'
KEY_PROPERTY_INTERACTION_METHOD = 'narupa.interaction.method'
KEY_PROPERTY_VELOCITY_RESET = 'narupa.interaction.velocity_reset'
KEY_PROPERTY_RENDERER = 'narupa.rendering.renderer'

INTERACTION_METHOD_DEFAULT = INTERACTION_SINGLE
VELOCITY_RESET_DEFAULT = False
RENDERER_DEFAULT = RENDERER_LIQUORICE
SELECTED_PARTICLE_IDS_DEFAULT = None


class NarupaImdSelection:
    """
    A selection of a group of particles in a Narupa simulation.

    :ivar selection_id: Blah
    """

    # The unique selection ID
    selection_id: str

    # User readable name of the selection
    selection_name: str

    # Set of particles indices that are in this selection
    selected_particle_ids: Set[int]

    # The interaction method for Narupa iMD
    interaction_method: str

    # Should the velocities be reset for this interaction
    velocity_reset: bool

    # The renderer to be used for this selection
    rendering_renderer = str

    @classmethod
    def from_dictionary(cls, dict: Dict):
        """
        Decode a dictionary into a selection.

        :param dict: A dictionary, such as json or a protobuf value.
        :return:
        """
        selection = cls(dict[KEY_SELECTION_ID], dict[KEY_SELECTION_NAME])
        selection.set_particles(
            get_nested_or_default(
                dict,
                SELECTED_PARTICLE_IDS_DEFAULT,
                KEY_SELECTION_SELECTED,
                KEY_SELECTED_PARTICLE_IDS,
            )
        )
        selection.interaction_method = get_nested_or_default(
            dict,
            INTERACTION_METHOD_DEFAULT,
            KEY_SELECTION_PROPERTIES,
            KEY_PROPERTY_INTERACTION_METHOD,
        )
        selection.velocity_reset = get_nested_or_default(
            dict,
            VELOCITY_RESET_DEFAULT,
            KEY_SELECTION_PROPERTIES,
            KEY_PROPERTY_VELOCITY_RESET,
        )
        selection.rendering_renderer = get_nested_or_default(
            dict,
            RENDERER_DEFAULT,
            KEY_SELECTION_PROPERTIES,
            KEY_PROPERTY_RENDERER,
        )

        return selection

    def __init__(self, id: str, name: str = 'Unnamed Selection'):
        """
        Create a new selection.

        :param id: The unique ID of the selection.
        :param name: A user-friendly name for the selection.
        """
        self.selection_id = id
        self.selection_name = name
        self.selected_particle_ids = set()

        self.interaction_method = INTERACTION_METHOD_DEFAULT
        self.velocity_reset = VELOCITY_RESET_DEFAULT
        self.rendering_renderer = RENDERER_DEFAULT

    def clear_particles(self):
        """
        Clear all particles in this selection.

        :return:
        """
        self.selected_particle_ids.clear()

    def set_particles(self, particle_ids: Iterable[int] = None):
        """
        Set the particles in this selection, replacing the previous selection.

        :param particle_ids:
        :return:
        """
        if particle_ids is None:
            particle_ids = set()
        self.selected_particle_ids = set(particle_ids)

    def add_particles(self, particle_ids: Iterable[int] = None):
        """
        Add more particles to this selection, appending the previous selection.

        :param particle_ids:
        :return:
        """
        if particle_ids is None:
            particle_ids = set()
        self.selected_particle_ids.update(particle_ids)

    def to_dictionary(self) -> Dict:
        """
        Encode this selection into a nested dictionary.

        :return: A dictionary representation of this selection.
        """
        return {
            KEY_SELECTION_ID: self.selection_id,
            KEY_SELECTION_NAME: self.selection_name,
            KEY_SELECTION_SELECTED: {
                KEY_SELECTED_PARTICLE_IDS: list(self.selected_particle_ids),
            },
            KEY_SELECTION_PROPERTIES: {
                KEY_PROPERTY_INTERACTION_METHOD: self.interaction_method,
                KEY_PROPERTY_VELOCITY_RESET: self.velocity_reset,
                KEY_PROPERTY_RENDERER: self.rendering_renderer
            }
        }


def get_nested_or_default(dict: Dict, default, *keys: Iterable[str]):
    """
    Iterate down a nested dictionary by accesssing subsequent keys, returning the default if at any point a key is not found.
    :param dict: The dictionary to iterate.
    :param default: The default value if a key is not found.
    :param keys: The keys to look up recursively in the dictionary.
    :return: The value found in the dictionary, or the default
    """
    for key in keys:
        try:
            dict = dict[key]
        except KeyError:
            return default
    return dict
