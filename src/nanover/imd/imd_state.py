"""
Module providing methods for storing ParticleInteractions in a StateDictionary.
"""

from collections.abc import Mapping
from logging import getLogger
from typing import Any

from nanover.imd.particle_interaction import ParticleInteraction
from nanover.utilities.change_buffers import DictionaryChange
from nanover.utilities.event import Event
from nanover.utilities.state_dictionary import StateDictionary

INTERACTION_PREFIX = "interaction."
VELOCITY_RESET_KEY = "imd.velocity_reset_available"


class ImdStateWrapper:
    """
    A wrapper around a StateDictionary that provides convenient methods for
    accessing and modifying ParticleInteractions.

    :param state_dictionary: The state dictionary to wrap.
    :param velocity_reset_available: Whether the dynamics this service is being
        used in allows velocity reset.
    """

    _interactions: dict[str, ParticleInteraction]

    def __init__(
        self,
        state_dictionary: StateDictionary,
        velocity_reset_available=False,
    ):
        self.state_dictionary = state_dictionary
        self.velocity_reset_available = velocity_reset_available

        self.state_dictionary.update_state(
            self,
            change=DictionaryChange(
                updates={VELOCITY_RESET_KEY: velocity_reset_available}
            ),
        )

        self.state_dictionary.content_updated.add_callback(self._on_state_updated)
        self._interactions = {}

        self.interaction_started = Event()
        self.interaction_stopped = Event()
        self.interaction_updated = Event()

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
        change = DictionaryChange(
            updates={
                interaction_id: interaction_to_dict(interaction),
            }
        )
        self.state_dictionary.update_state(None, change)

    def remove_interaction(self, interaction_id: str):
        assert interaction_id.startswith(INTERACTION_PREFIX)
        change = DictionaryChange(removals={interaction_id})
        self.state_dictionary.update_state(None, change)

    def clear_interactions(self):
        for interaction_id in self.active_interactions:
            self.remove_interaction(interaction_id)

    @property
    def active_interactions(self) -> dict[str, ParticleInteraction]:
        """
        The current dictionary of active interactions, keyed by interaction id.

        :return: A copy of the dictionary of active interactions.
        """
        return self._interactions.copy()

    def _on_state_updated(self, access_token, change: DictionaryChange):
        for removed_key in change.removals:
            if removed_key.startswith(INTERACTION_PREFIX):
                interaction = self._interactions.pop(removed_key, None)
                if interaction is not None:
                    self.interaction_stopped.invoke(
                        key=removed_key, interaction=interaction
                    )
        for key, value in change.updates.items():
            if key.startswith(INTERACTION_PREFIX):
                try:
                    interaction = dict_to_interaction(value)
                    started = key not in self._interactions
                    self._interactions[key] = interaction
                    if started:
                        self.interaction_started.invoke(
                            key=key, interaction=interaction
                        )
                    self.interaction_updated.invoke(key=key, interaction=interaction)
                except Exception:
                    getLogger().exception(
                        f"Didn't understand '{value}' ({key}) as interaction"
                    )


def interaction_to_dict(interaction: ParticleInteraction) -> dict[str, Any]:
    try:
        # properties with the same key as the builtins will be discarded
        # discussion: https://gitlab.com/intangiblerealities/nanover-server-py/-/merge_requests/182#note_374156050
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


def dict_to_interaction(dictionary: Mapping[str, Any]) -> ParticleInteraction:
    kwargs: dict = dict(**dictionary)
    if "particles" in kwargs:
        kwargs["particles"] = [int(i) for i in kwargs["particles"]]
    return ParticleInteraction(**kwargs)
