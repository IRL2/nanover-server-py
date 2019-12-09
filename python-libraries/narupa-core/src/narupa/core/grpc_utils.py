# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Utilities for working with GRPC, particularly wrappers with more descriptive
names.
"""

from typing import Callable
import grpc


class RpcAlreadyTerminatedError(Exception):
    pass


def subscribe_channel_connectivity_change(
        channel: grpc.Channel,
        callback: Callable[[grpc.ChannelConnectivity], None],
        force_connection: bool = False,
):
    """
    Subscribe to channel connectivity changes with a callback that is called
    with the channel's latest connectivity status. Optionally force the channel
    to begin connecting instead of waiting for an RPC attempt.

    .. warning::
        On channel close, all subscriptions are removed without state change.
    """
    channel.subscribe(callback, force_connection)


def subscribe_rpc_termination(
        context: grpc.RpcContext,
        callback: Callable[[], None]
):
    """
    Subscribe the termination of an RPC with the given callback.
    
    :raises RpcAlreadyTerminatedError: if the callback will not be used because
        termination has already occurred.
    """
    added = context.add_callback(callback)
    if not added:
        raise RpcAlreadyTerminatedError()
