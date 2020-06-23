# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing an implementation of an IMD service.
"""
from typing import Dict, Callable, Optional, Iterable

from narupa.imd.imd_client import _interaction_to_dict, _dict_to_interaction
from narupa.state.state_dictionary import StateDictionary
from narupa.utilities.change_buffers import DictionaryChange
from narupa.imd.particle_interaction import ParticleInteraction
from narupa.protocol.imd import (
    InteractiveMolecularDynamicsServicer,
    add_InteractiveMolecularDynamicsServicer_to_server,
)

IMD_SERVICE_NAME = "imd"


class ImdService(InteractiveMolecularDynamicsServicer):
    """
    An implementation of an IMD service, that keeps track of interactions
    produced by clients that can be consumed by an interactive molecular
    dynamics simulation.

    :param state_dictionary: If provided, interactions will be stored in this
        existing state dictionary instead of a newly created one.
    :param velocity_reset_enabled: Whether the dynamics this service is being
        used in supports velocity reset.
    :param number_of_particles: If set, allows the IMD service to validate the
        particle indices in an interaction are valid. If not set, it is up to
        consumers of the interactions to validate.
    """

    velocity_reset_enabled: bool
    number_of_particles: Optional[int]

    def __init__(
            self,
            state_dictionary: Optional[StateDictionary] = None,
            velocity_reset_enabled=False,
            number_of_particles=None,
    ):
        self.state_dictionary = state_dictionary if state_dictionary is not None else StateDictionary()
        self.name: str = IMD_SERVICE_NAME
        self.add_to_server_method: Callable = add_InteractiveMolecularDynamicsServicer_to_server
        self.velocity_reset_enabled = velocity_reset_enabled
        self.number_of_particles = number_of_particles

    def close(self):
        pass

    def insert_interaction(self, interaction_id: str, interaction: ParticleInteraction):
        key = 'interaction.'+interaction_id
        value = _interaction_to_dict(interaction)
        change = DictionaryChange(updates={key: value})
        self.state_dictionary.update_state(None, change)

    def remove_interaction(self, interaction_id: str):
        self.remove_interactions_by_ids([interaction_id])

    def remove_interactions_by_ids(self, interaction_ids: Iterable[str]):
        change = DictionaryChange(
            removals=['interaction.'+key for key in interaction_ids],
        )
        self.state_dictionary.update_state(None, change)

    @property
    def active_interactions(self) -> Dict[str, ParticleInteraction]:
        """
        The current dictionary of active interactions, keyed by interaction id.

        :return: A copy of the dictionary of active interactions.
        """
        with self.state_dictionary.lock_content() as content:
            return {
                key: _dict_to_interaction(value)
                for key, value in content.items()
                if key.startswith('interaction.')
            }
