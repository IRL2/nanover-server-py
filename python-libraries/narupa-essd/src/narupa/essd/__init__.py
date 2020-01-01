# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from narupa.essd.server import DiscoveryServer
from narupa.essd.client import DiscoveryClient
from narupa.essd.servicehub import ServiceHub

"""
Module providing extremely simple service discovery over UDP.

Messages are broadcast over IPv4 UDP to a specific port. The message consists of a JSON payload encoded with
UTF-8.

An example message is: 

.. code 
  {
    "name": "Example Narupa Service Hub", 
    "address": "localhost", 
    "id": "12345", 
    "essd_version": "1.0.0", 
    "services": 
    {
      "imd": 54322, 
      "trajectory": 54323
    }
  }
  
The :class:`DiscoveryServer` class can be used to broadcast these service definitions, enabling clients, such as 
the :class:`DiscoveryClient` class, to find them by listening to messages on the specified broadcasting port. 

"""
__version__ = '1.0.0'
