# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from contextlib import contextmanager
from typing import Dict, Iterable, Set, Union, TypeVar, Optional, cast

from narupa.utilities.event import Event

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
KEY_PROPERTY_HIDE = 'narupa.rendering.hide'

INTERACTION_METHOD_DEFAULT = INTERACTION_SINGLE
VELOCITY_RESET_DEFAULT = False
RENDERER_DEFAULT = RENDERER_LIQUORICE
SELECTED_PARTICLE_IDS_DEFAULT = None


class RenderingSelection:
    """

    A local representation of a selection in a NarupaIMD simulation.

    Not intended to be created directly, but rather via the :class:`NarupaImdClient`, which handles
    synchronising the state of the selections with the server.

    Changes to this selection are only applied when a call to
    flush_changes() is made, or modifications are made within
    a modify(): context manager.
    """

    @classmethod
    def from_dictionary(cls, dict: Dict) -> 'RenderingSelection':
        """
        Decode a dictionary into a selection.

        :param dict: A dictionary, such as json or a protobuf value.
        :return: A selection decoded from the given dictionary
        """
        try:
            selection = cls(dict[KEY_SELECTION_ID], dict[KEY_SELECTION_NAME])
        except KeyError:
            raise ValueError("Selection is missing ID or name")
        selection.read_from_dictionary(dict)

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
        self.renderer = RENDERER_DEFAULT
        self.hide = False

        self._updated = Event()
        self._removed = Event()

    @contextmanager
    def modify(self):
        """
        Gives a context in which the selection can have multiple modifications made to it, and which calls
        flush_changes() when the context is left
        """
        yield self
        self.flush_changes()

    def flush_changes(self):
        """
        Flushes all changes by notifying any listeners that this selection has changed.
        """
        self.updated.invoke(self)

    def remove(self):
        """
        Remove this selection from the server.
        """
        self.removed.invoke(self)

    def clear_particles(self):
        """
        Clear all particles in this selection.

        """
        self.selected_particle_ids.clear()

    def set_particles(self, particle_ids: Optional[Iterable[int]] = None):
        """
        Set the particles in this selection, replacing the previous selection.

        :param particle_ids:
        """
        if particle_ids is None:
            self.selected_particle_ids = set()
        else:
            self.selected_particle_ids = set(particle_ids)

    def add_particles(self, particle_ids: Iterable[int]):
        """
        Add more particles to this selection, appending the previous selection.

        :param particle_ids:
        """
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
                KEY_SELECTED_PARTICLE_IDS: list(self.selected_particle_ids) \
                    if self.selected_particle_ids is not None \
                    else None,
            },
            KEY_SELECTION_PROPERTIES: {
                KEY_PROPERTY_INTERACTION_METHOD: self.interaction_method,
                KEY_PROPERTY_VELOCITY_RESET: self.velocity_reset,
                KEY_PROPERTY_RENDERER: self.renderer,
                KEY_PROPERTY_HIDE: self.hide
            }
        }

    @property
    def selection_id(self) -> str:
        """
        The unique selection ID
        """
        return self._selection_id

    @selection_id.setter
    def selection_id(self, value: str):
        self._selection_id = value

    @property
    def selection_name(self) -> str:
        """
        User readable name of the selection
        """
        return self._selection_name

    @selection_name.setter
    def selection_name(self, value: str):
        self._selection_name = value

    @property
    def selected_particle_ids(self) -> Set[int]:
        """
        Set of particles indices that are in this selection
        """
        return self._selected_particle_ids

    @selected_particle_ids.setter
    def selected_particle_ids(self, value: Set[int]):
        self._selected_particle_ids = set(value)

    @property
    def updated(self) -> Event:
        """
        Event which is invoked when modifications to this selections are applied
        """
        return self._updated

    @property
    def removed(self) -> Event:
        """
        Event which is invoked when this selection is removed
        """
        return self._removed

    @property
    def interaction_method(self) -> str:
        """
        The interaction method for Narupa iMD
        """
        return self._interaction_method

    @interaction_method.setter
    def interaction_method(self, value: str):
        self._interaction_method = value

    @property
    def velocity_reset(self) -> bool:
        """
        Should the velocities be reset for this interaction
        """
        return self._velocity_reset

    @velocity_reset.setter
    def velocity_reset(self, value: bool):
        self._velocity_reset = value

    @property
    def hide(self) -> bool:
        """
        Should this renderer be hidden
        """
        return self._hide

    @hide.setter
    def hide(self, value: bool):
        self._hide = value

    @property
    def renderer(self) -> Union[str, Dict]:
        """
        The renderer to be used for this selection. Either a string name of a predefined visualiser,
        or a dictionary describing one.
        """
        return self._renderer

    @renderer.setter
    def renderer(self, value: Union[str, Dict]):
        self._renderer = value

    def read_from_dictionary(self, dict: Dict):
        self.read_selected_particles_from_dictionary(dict)
        self.read_interaction_method_from_dictionary(dict)
        self.read_velocity_reset_from_dictionary(dict)
        self.read_renderer_from_dictionary(dict)
        self.read_hide_from_dictionary(dict)

    def read_selected_particles_from_dictionary(self, dict: Dict):
        self.set_particles(get_nested_or_default(
            dict,
            SELECTED_PARTICLE_IDS_DEFAULT,
            KEY_SELECTION_SELECTED,
            KEY_SELECTED_PARTICLE_IDS,
        ))

    def read_interaction_method_from_dictionary(self, dict: Dict):
        self.interaction_method = get_nested_or_default(
            dict,
            INTERACTION_METHOD_DEFAULT,
            KEY_SELECTION_PROPERTIES,
            KEY_PROPERTY_INTERACTION_METHOD,
        )

    def read_velocity_reset_from_dictionary(self, dict: Dict):
        self.velocity_reset = get_nested_or_default(
            dict,
            VELOCITY_RESET_DEFAULT,
            KEY_SELECTION_PROPERTIES,
            KEY_PROPERTY_VELOCITY_RESET,
        )

    def read_renderer_from_dictionary(self, dict: Dict):
        self.renderer = get_nested_or_default(
            dict,
            RENDERER_DEFAULT,
            KEY_SELECTION_PROPERTIES,
            KEY_PROPERTY_RENDERER,
        )

    def read_hide_from_dictionary(self, dict: Dict):
        self.hide = get_nested_or_default(
            dict,
            False,
            KEY_SELECTION_PROPERTIES,
            KEY_PROPERTY_HIDE,
        )


T = TypeVar('T')


def get_nested_or_default(dict: Dict, default: T, *keys: Iterable[str]) -> T:
    """
    Iterate down a nested dictionary by accessing subsequent keys, returning the default if at any point a key is not found.

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
    return cast(T, dict)
