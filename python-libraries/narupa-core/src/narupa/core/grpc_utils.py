# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Utilities for working with GRPC, particularly wrappers with more descriptive
names.
"""

from typing import Callable
import grpc


class RpcContextAlreadyTerminatedError(Exception):
    pass


def subscribe_channel_status_change(
        channel: grpc.Channel,
        callback: Callable[[grpc.ChannelConnectivity], None],
        force_connection: bool,
):
    channel.subscribe(callback, force_connection)


def subscribe_rpc_termination(
        context: grpc.RpcContext,
        callback: Callable[[], None]
):
    if not context.add_callback(callback):
        raise RpcContextAlreadyTerminatedError()
