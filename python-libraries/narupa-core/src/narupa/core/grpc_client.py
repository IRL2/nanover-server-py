"""
Module providing a base class for gRPC clients.
"""
# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from concurrent import futures
from typing import Optional

import grpc
from narupa.core import DEFAULT_CONNECT_ADDRESS


class GrpcClient:
    """
    A base class for GRPC clients that handles service connection and client
    closing.

    :param channel: An existing :class:`grpc.Channel` to use to establish a connection.
    :param make_channel_owner: Whether to make this client take ownership of
     ensuring the channel is closed upon disconnection of this client.
    """

    channel: grpc.Channel
    threads: futures.ThreadPoolExecutor
    _channel_owner: bool
    DEFAULT_CONNECTION_PORT: int = 54321

    def __init__(
        self, *, channel: grpc.Channel, make_channel_owner: bool = False, **kwargs
    ):
        # TODO a channel could be wrapped into a more general gRPC connection,
        #  as in the C# implementation.
        self.channel = channel
        self._channel_owner = make_channel_owner
        self.threads = futures.ThreadPoolExecutor(max_workers=10)

    @classmethod
    def insecure_channel(
        cls, *, address: Optional[str] = None, port: Optional[int] = None, **kwargs
    ):
        """
        Create an insecure connection at the given address and port.

        :param address: The URL or IP address of the service to connect to.
        :param port: The port on which to connect.
        :return: An instantiation of a client connected insecurely at the given address and port.
        """
        address = address or DEFAULT_CONNECT_ADDRESS
        port = port or cls.DEFAULT_CONNECTION_PORT
        channel = grpc.insecure_channel(f"{address}:{port}")
        return cls(channel=channel, make_channel_owner=True, **kwargs)

    @property
    def is_channel_owner(self) -> bool:
        """
        Indicates whether this client is responsible for managing the underlying channel.

        :return: True if this client is responsible for managing the underlying channel,
            False otherwise.
        """
        return self._channel_owner

    def close(self):
        """
        Shutdown all threads and close the underlying channel
        if the client has been given that responsibility.
        """
        if self._channel_owner:
            self.channel.close()
        self.threads.shutdown(wait=False)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
