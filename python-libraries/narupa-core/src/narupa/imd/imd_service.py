# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing an implementation of an IMD service.
"""
from typing import Dict, Callable, Optional, Iterable

import grpc
from narupa.state.state_dictionary import StateDictionary
from narupa.state.state_service import StateService
from narupa.utilities.grpc_utilities import (
    RpcAlreadyTerminatedError,
    subscribe_rpc_termination,
)
from narupa.utilities.change_buffers import DictionaryChangeMultiView, \
    DictionaryChange

from narupa.imd.particle_interaction import ParticleInteraction
from narupa.protocol.imd import (
    InteractiveMolecularDynamicsServicer,
    InteractionEndReply,
    SubscribeInteractionsRequest,
    InteractionsUpdate,
    ParticleInteraction as ParticleInteractionMessage,
    add_InteractiveMolecularDynamicsServicer_to_server)

IMD_SERVICE_NAME = "imd"


class ImdService(InteractiveMolecularDynamicsServicer):
    """
    An implementation of an IMD service, that keeps track of interactions
    produced by clients that can be consumed by an interactive molecular
    dynamics simulation.

    :param callback: A callback to be used whenever an interaction is published
        or updated.
    :param velocity_reset_enabled: Whether the dynamics this service is being
        used in supports velocity reset.
    :param number_of_particles: If set, allows the IMD service to validate the
        particle indices in an interaction are valid. If not set, it is up to
        consumers of the interactions to validate.
    """

    _interaction_updated_callback: Optional[
        Callable[[ParticleInteraction], None]]
    velocity_reset_enabled: bool
    number_of_particles: Optional[int]

    def __init__(
            self,
            state_dictionary: Optional[StateDictionary] = None,
            callback=None,
            velocity_reset_enabled=False,
            number_of_particles=None,
    ):
        self._state_dictionary = state_dictionary or StateDictionary()
        self.name: str = IMD_SERVICE_NAME
        self.add_to_server_method: Callable = add_InteractiveMolecularDynamicsServicer_to_server
        self._interaction_updated_callback = callback
        self.velocity_reset_enabled = velocity_reset_enabled
        self.number_of_particles = number_of_particles

    def close(self):
        pass

    def PublishInteraction(self, request_iterator, context):
        """
        Publishes an interaction to be consumed by an IMD simulation.

        A stream represents the lifetime of a single interaction, and
        is used to gracefully determine when to halt an interaction.
        An interaction will continue until the stream is closed, at which point
        it will be deleted.

        The implementation is stateless, and so each Interaction message should
        contain the full set of fields.

        :param request_iterator: Stream of interaction requests.
        :param context: gRPC context.
        :return: :class: InteractionEndReply, indicating the successful interaction.
        """
        try:
            self._publish_interaction_stream(request_iterator)
        except KeyError as e:
            context.set_code(grpc.StatusCode.PERMISSION_DENIED)
            message = str(e)
            context.set_details(message)
        except ValueError as e:
            context.set_code(grpc.StatusCode.INVALID_ARGUMENT)
            message = str(e)
            context.set_details(message)
        return InteractionEndReply()

    def SubscribeInteractions(
            self,
            request: SubscribeInteractionsRequest,
            context,
    ):
        """
        Provides a stream of updates to interactions.
        """
        interval = request.update_interval
        with self._state_dictionary.get_change_buffer() as change_buffer:
            try:
                subscribe_rpc_termination(context, change_buffer.freeze)
            except RpcAlreadyTerminatedError:
                return
            for changes, removals in change_buffer.subscribe_changes(interval):
                yield _changes_to_interactions_update_message(changes, removals)

    def insert_interaction(self, interaction: ParticleInteraction):
        change = DictionaryChange(
            updates={'interaction.'+interaction.interaction_id: interaction},
            removals=[],
        )
        self._state_dictionary.update_state(None, change)

        if self._interaction_updated_callback is not None:
            self._interaction_updated_callback(interaction)

    def remove_interaction(self, interaction: ParticleInteraction):
        change = DictionaryChange(
            updates={},
            removals=['interaction.'+interaction.interaction_id],
        )
        self._state_dictionary.update_state(None, change)

    def remove_interactions_by_ids(self, interaction_ids: Iterable[str]):
        change = DictionaryChange(
            updates={},
            removals=['interaction.'+key for key in interaction_ids],
        )
        self._state_dictionary.update_state(None, change)

    @property
    def active_interactions(self) -> Dict[str, ParticleInteraction]:
        """
        The current dictionary of active interactions, keyed by interaction id.

        :return: A copy of the dictionary of active interactions.
        """
        with self._state_dictionary.lock_content() as state:
            return {
                key: value
                for key, value in state.items()
                if key.startswith('interaction.')
            }

    def set_interaction_updated_callback(self, callback: Callable[
        [ParticleInteraction], None]):
        """
        Sets the callback to be used whenever an interaction is received.

        :param callback: Method to be called, taking the received particle interaction as an argument.
        """
        self._interaction_updated_callback = callback

    def _interaction_id_exists(self, interaction_id: str) -> bool:
        with self._state_dictionary.lock_content() as content:
            return 'interaction.'+interaction_id in content

    def _publish_interaction_stream(self, request_iterator):
        interactions_in_request = set()

        def validate_interaction_id(interaction: ParticleInteractionMessage):
            id_owned = interaction.interaction_id in interactions_in_request
            id_exists = self._interaction_id_exists(interaction.interaction_id)
            if not id_owned and id_exists:
                raise KeyError('The given player_id and interaction_id is '
                               'already in use in another interaction.')

        try:
            for interaction_message in request_iterator:
                interaction = ParticleInteraction.from_proto(
                    interaction_message)

                validate_interaction_id(interaction_message)
                self._validate_interaction(interaction)
                interactions_in_request.add(interaction_message.interaction_id)

                self.insert_interaction(interaction)
        finally:
            self.remove_interactions_by_ids(interactions_in_request)

    def _validate_interaction(self, interaction):
        self._validate_particle_range(interaction)
        self._validate_velocity_reset(interaction)

    def _validate_velocity_reset(self, interaction):
        if interaction.reset_velocities and not self.velocity_reset_enabled:
            raise ValueError(
                'This simulation does not support interactive velocity reset.')

    def _validate_particle_range(self, interaction):
        if interaction.particles is None or len(interaction.particles) == 0:
            return
        max_index = max(interaction.particles)
        out_of_range = self.number_of_particles is not None and max_index >= self.number_of_particles
        if out_of_range:
            raise ValueError('Range of particles selected invalid.')


def _changes_to_interactions_update_message(changes, removals):
    """
    Creates a protobuf InteractionsUpdate message from changes and removals to
    an interactions dictionary.
    """
    message = InteractionsUpdate()
    protos = [interaction.proto for interaction in changes.values()]
    message.updated_interactions.extend(protos)
    message.removals.extend(interaction_id[len('interaction.'):] for interaction_id in removals)
    return message
