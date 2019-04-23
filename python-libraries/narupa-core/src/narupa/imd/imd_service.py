# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing an implementation of an IMD service.
"""
import grpc

from narupa.protocol.imd import InteractiveMolecularDynamicsServicer, InteractionEndReply


class ImdService(InteractiveMolecularDynamicsServicer):
    """
    An implementation of an IMD service, that keeps track of interactions
    produced by clients that can be consumed by an interactive molecular dynamics
    simulation.

    :param callback: A callback to be used whenever an interaction is published or updated.

    """

    def __init__(self, callback=None):
        self._interactions = {}
        self._callback = callback

    def PublishInteraction(self, request_iterator, context):
        """
        Publishes an interaction to be consumed by an IMD simulation.

        A stream represents the lifetime of a single interaction, and
        is used to gracefully determine when to halt an interaction.
        An interaction will continue until the stream is closed, at which point it
        will be deleted.

        The implementation is stateless, and so each Interaction message should contain
        the full set of fields.

        :param request_iterator:
        :param context:
        :return:
        """
        active_interactions = set()
        try:
            self._process_interactions(active_interactions, request_iterator, context)
        except KeyError:
            context.set_code(grpc.StatusCode.PERMISSION_DENIED)
            message = "Tried to create an interaction with a player ID and device ID combination that's already in use"
            context.set_details(message)
        finally:
            # clean up all the interactions after the stream
            for key in active_interactions:
                del self.interactions[key]
        return InteractionEndReply()

    @property
    def interactions(self):
        """
        The current dictionary of interactions, keyed by a tuple consisting of player_id and interaction_id.
        :return:
        """
        return self._interactions

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
        :return:
        """
        self._callback = callback

    def _process_interactions(self, active_interactions, request_iterator, context):
        for interaction in request_iterator:
            key = ImdService.get_key(interaction)
            if key not in active_interactions and key in self.interactions:
                raise KeyError
            active_interactions.add(key)
            self.interactions[key] = interaction
            if self._callback is not None:
                self._callback()
