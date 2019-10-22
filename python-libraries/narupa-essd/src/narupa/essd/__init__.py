# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from narupa.essd.server import DiscoveryServer
from narupa.essd.client import DiscoveryClient

"""
Module providing extremely simple service discovery over UDP.

Messages are broadcast over IPv4 UDP to a specific port. The message consists of a json payload encoded with
UTF-8.

An example message is: 

"""
__version__ = '1.0.0'
