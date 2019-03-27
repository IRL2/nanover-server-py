Working space for the new narupa protocol and associated implementations using gRPC.

This is very much a work in progress!

## Getting Started

The project requires python 3.6+ 

The `protocol` folder contains the definitions of the gRPC services. The `compile.sh` script can be used to generate
the python bindings.

The `project` folder contains implementations of the narupa services in the relevant languages.

The `narupy` folder contains some demo server implementations:

* OpenMMTest.py - Downloads a PDB file, runs an OpenMM simulation and streams it to clients.

The `project/csharp/Narupa.Protocol folder contains C# implementations of clients for receiving trajectories and structures.

## Coming Soon

* Interactive MD servers
* Multiplayer
* Lobby
* Trajectory serving
