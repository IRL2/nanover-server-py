# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Simple example of running a server and registering service hub.
"""
import logging

from narupa.essd.server import DiscoveryServer
from narupa.essd.servicehub import ServiceHub

logging.basicConfig(level=logging.DEBUG)

with DiscoveryServer() as server:

    hub = ServiceHub(name="Example Narupa service", address="[::]")
    print(f'Registering hub: {hub}')
    server.register_service(hub)
    try:
        while True:
            pass
    except KeyboardInterrupt:
        print("Closing due to keyboard interrupt.")
