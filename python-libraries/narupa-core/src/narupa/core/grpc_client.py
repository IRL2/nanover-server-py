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

    :param address: The URL or IP address of the service to connect to.
    :param port: The port on which to connect.
    :param stub: The GRPC service stub class.
    """

    def __init__(self, *, address: Optional[str] = None, port: int, stub: Optional = None):
        address = address or DEFAULT_CONNECT_ADDRESS
        self.channel = grpc.insecure_channel(f"{address}:{port}")
        if stub is not None:
            self.stub = stub(self.channel)
        self.threads = futures.ThreadPoolExecutor(max_workers=10)

    def close(self):
        """
        Close the channel and shutdown all threads.
        """
        self.channel.close()
        self.threads.shutdown(wait=False)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
