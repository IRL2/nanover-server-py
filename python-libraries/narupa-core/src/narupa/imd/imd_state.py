# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing methods for storing ParticleInteractions in a StateDictionary.
"""
from typing import Dict, Any

from narupa.state.state_dictionary import StateDictionary
from narupa.utilities.change_buffers import DictionaryChange
from narupa.imd.particle_interaction import ParticleInteraction
from narupa.utilities.protobuf_utilities import Serializable

IMD_SERVICE_NAME = "imd"
INTERACTION_PREFIX = 'interaction.'
VELOCITY_RESET_KEY = 'imd.velocity_reset_available'


class ImdStateWrapper:
    """
    A wrapper around a StateDictionary that provides convenient methods for
    accessing and modifying ParticleInteractions.

    :param state_dictionary: The state dictionary to wrap.
    :param velocity_reset_available: Whether the dynamics this service is being
        used in allows velocity reset.
    """
    _interactions: Dict[str, ParticleInteraction]

    def __init__(
            self,
            state_dictionary: StateDictionary,
            velocity_reset_available=False,
    ):
        self.state_dictionary = state_dictionary
        self.velocity_reset_available = velocity_reset_available

        self.state_dictionary.update_locks(
            self,
            acquire={VELOCITY_RESET_KEY: None},
        )
        self.state_dictionary.update_state(
            self,
            change=DictionaryChange(updates={
                VELOCITY_RESET_KEY: velocity_reset_available
            }),
        )

        self.state_dictionary.content_updated.add_callback(self._on_state_updated)
        self._interactions = {}

    @property
    def velocity_reset_available(self):
        with self.state_dictionary.lock_content() as state:
            return state[VELOCITY_RESET_KEY]

    @velocity_reset_available.setter
    def velocity_reset_available(self, value: bool):
        change = DictionaryChange(updates={VELOCITY_RESET_KEY: value})
        self.state_dictionary.update_state(self, change)

    def insert_interaction(self, interaction_id: str, interaction: ParticleInteraction):
        assert interaction_id.startswith(INTERACTION_PREFIX)
        change = DictionaryChange(updates={
            interaction_id: interaction_to_dict(interaction),
        })
        self.state_dictionary.update_state(None, change)

    def remove_interaction(self, interaction_id: str):
        assert interaction_id.startswith(INTERACTION_PREFIX)
        change = DictionaryChange(removals=set([interaction_id]))
        self.state_dictionary.update_state(None, change)

    @property
    def active_interactions(self) -> Dict[str, ParticleInteraction]:
        """
        The current dictionary of active interactions, keyed by interaction id.

        :return: A copy of the dictionary of active interactions.
        """
        return self._interactions.copy()

    def _on_state_updated(self, access_token, change: DictionaryChange):
        for removed_key in change.removals:
            if removed_key.startswith(INTERACTION_PREFIX) and removed_key in self._interactions:
                del self._interactions[removed_key]
        for key, value in change.updates.items():
            if key.startswith(INTERACTION_PREFIX):
                self._interactions[key] = dict_to_interaction(value)


def interaction_to_dict(interaction: ParticleInteraction) -> Dict[str, Serializable]:
    try:
        # properties with the same key as the builtins will be discarded
        # dicussion: https://gitlab.com/intangiblerealities/narupa-protocol/-/merge_requests/182#note_374156050
        return {
            **interaction.properties,
            "position": [float(f) for f in interaction.position],
            "particles": [int(i) for i in interaction.particles],
            "interaction_type": interaction.interaction_type,
            "scale": interaction.scale,
            "mass_weighted": interaction.mass_weighted,
            "reset_velocities": interaction.reset_velocities,
            "max_force": interaction.max_force,
        }
    except AttributeError as e:
        raise TypeError from e


def dict_to_interaction(dictionary: Dict[str, Any]) -> ParticleInteraction:
    kwargs = dict(**dictionary)
    if 'particles' in kwargs:
        kwargs['particles'] = [int(i) for i in kwargs['particles']]
    return ParticleInteraction(**kwargs)
