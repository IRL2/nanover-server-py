# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing an implementation of an IMD service.
"""
from typing import Dict, Any

from narupa.state.state_dictionary import StateDictionary
from narupa.utilities.change_buffers import DictionaryChange
from narupa.imd.particle_interaction import ParticleInteraction

IMD_SERVICE_NAME = "imd"
INTERACTION_PREFIX = 'interaction.'
VELOCITY_RESET_KEY = 'imd.velocity_reset_enabled'


class ImdService:
    """
    An implementation of an IMD service, that keeps track of interactions
    produced by clients that can be consumed by an interactive molecular
    dynamics simulation.

    :param state_dictionary: The state dictionary to store interactions in.
    :param velocity_reset_enabled: Whether the dynamics this service is being
        used in supports velocity reset.
    :param number_of_particles: If set, allows the IMD service to validate the
        particle indices in an interaction are valid. If not set, it is up to
        consumers of the interactions to validate.
    """

    def __init__(
            self,
            state_dictionary: StateDictionary,
            velocity_reset_enabled=False,
    ):
        self.state_dictionary = state_dictionary
        self.velocity_reset_enabled = velocity_reset_enabled

        # prototype attempt at faster ImdServer.interactions
        def on_content_updated(change: DictionaryChange, **_):
            self._interactions.update({
                key[len(INTERACTION_PREFIX):]: dict_to_interaction(value)
                for key, value in change.updates.items()
                if key.startswith(INTERACTION_PREFIX)
            })
            for key in change.removals:
                self._interactions.pop(key[len(INTERACTION_PREFIX):], None)

        self._interactions: Dict[str, ParticleInteraction] = {}
        self.state_dictionary.content_updated.add_callback(on_content_updated)

        self.state_dictionary.update_locks(
            self,
            acquire={VELOCITY_RESET_KEY: None},
        )
        self.state_dictionary.update_state(
            self,
            change=DictionaryChange(updates={
                VELOCITY_RESET_KEY: velocity_reset_enabled
            }),
        )

    @property
    def velocity_reset_enabled(self):
        with self.state_dictionary.lock_content() as state:
            return state[VELOCITY_RESET_KEY]

    @velocity_reset_enabled.setter
    def velocity_reset_enabled(self, value: bool):
        change = DictionaryChange(updates={VELOCITY_RESET_KEY: value})
        self.state_dictionary.update_state(self, change)

    def close(self):
        pass

    def insert_interaction(self, interaction_id: str, interaction: ParticleInteraction):
        change = DictionaryChange(updates={
            INTERACTION_PREFIX + interaction_id: interaction_to_dict(interaction),
        })
        self.state_dictionary.update_state(None, change)

    def remove_interaction(self, interaction_id: str):
        change = DictionaryChange(removals=[INTERACTION_PREFIX + interaction_id])
        self.state_dictionary.update_state(None, change)

    @property
    def active_interactions(self) -> Dict[str, ParticleInteraction]:
        """
        The current dictionary of active interactions, keyed by interaction id.

        :return: A copy of the dictionary of active interactions.
        """
        return self._interactions


def interaction_to_dict(interaction: ParticleInteraction):
    try:
        return {
            "position": [float(f) for f in interaction.position],
            "particles": [int(i) for i in interaction.particles],
            "interaction_type": interaction.type,
            "scale": interaction.scale,
            "mass_weighted": interaction.mass_weighted,
            "reset_velocities": interaction.reset_velocities,
            "properties": dict(**interaction.properties),
            "max_force": interaction.max_force,
        }
    except AttributeError as e:
        raise TypeError from e


def dict_to_interaction(dictionary: Dict[str, Any]) -> ParticleInteraction:
    kwargs = dict(**dictionary)
    kwargs['particles'] = [int(i) for i in kwargs['particles']]
    return ParticleInteraction(**kwargs)
