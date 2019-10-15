# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing an implementation of an IMD service.
"""
from threading import Lock
from typing import Iterable, List, Dict, Tuple

import grpc

from narupa.imd.particle_interaction import ParticleInteraction
from narupa.protocol.imd import InteractiveMolecularDynamicsServicer, InteractionEndReply


class ImdService(InteractiveMolecularDynamicsServicer):
    """
    An implementation of an IMD service, that keeps track of interactions
    produced by clients that can be consumed by an interactive molecular
    dynamics simulation.

    :param callback: A callback to be used whenever an interaction is published
    or updated.
    :param velocity_reset_enabled: Whether the dynamics this service is being used in
    supports velocity reset.
    :param number_of_particles: If set, allows the IMD service to validate the particle indices
    in an interaction are valid. If not set, it is up to consumers of the
    interactions to validate.

    """

    def __init__(self,
                 callback=None,
                 velocity_reset_enabled=False,
                 number_of_particles=None):
        self._interactions = {}
        self._callback = callback
        self.velocity_reset_enabled = velocity_reset_enabled
        self.number_of_particles = number_of_particles
        self._interaction_lock = Lock()

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
        interactions_in_request = set()
        try:
            self._publish_interaction(
                interactions_in_request,
                request_iterator,
                context,
            )
        except KeyError as e:
            context.set_code(grpc.StatusCode.PERMISSION_DENIED)
            message = str(e)
            context.set_details(message)
        except ValueError as e:
            context.set_code(grpc.StatusCode.INVALID_ARGUMENT)
            message = str(e)
            context.set_details(message)
        finally:
            # clean up all the interactions after the stream
            with self._interaction_lock:
                for key in interactions_in_request:
                    del self._interactions[key]
        return InteractionEndReply()

    @property
    def active_interactions(self) -> Dict[Tuple[str, str], ParticleInteraction]:
        """
        The current dictionary of active interactions, keyed by player id and interaction id.
        :return A copy of the dictionary of active interactions.
        """
        with self._interaction_lock:
            return dict(self._interactions)

    @staticmethod
    def get_key(interaction):
        """
        Gets the key to use with an interaction.
        :param interaction:
        :return: key formed of a tuple of player_id and interaction_id.
        """
        return interaction.player_id, interaction.interaction_id

    def set_callback(self, callback):
        """
        Sets the callback to be used whenever an interaction is received.
        :param callback: Method to be called
        """
        self._callback = callback

    def _publish_interaction(self, interactions_in_request, request_iterator,
                             context):
        for interaction in request_iterator:
            key = self.get_key(interaction)
            with self._interaction_lock:
                if key not in interactions_in_request and key in self._interactions:
                    raise KeyError('The given player_id and interaction_id '
                                   'is already in use in another '
                                   'interaction.')
                # convert to wrapped interaction type.
                interaction = ParticleInteraction.from_proto(interaction)
                self._validate_interaction(interaction)

                interactions_in_request.add(key)
                self._interactions[key] = interaction
            if self._callback is not None:
                self._callback()

    def _validate_interaction(self, interaction):
        self._validate_particle_range(interaction)
        self._validate_velocity_reset(interaction)

    def _validate_velocity_reset(self, interaction):
        if interaction.reset_velocities and not self.velocity_reset_enabled:
            raise ValueError('This simulation does not support interactive velocity reset.')

    def _validate_particle_range(self, interaction):
        if interaction.particles is None or len(interaction.particles) == 0:
            return
        max_index = max(interaction.particles)
        out_of_range = self.number_of_particles is not None and max_index >= self.number_of_particles
        if out_of_range:
            raise ValueError('Range of particles selected invalid.')
