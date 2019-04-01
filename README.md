Working space for the new narupa protocol and associated implementations using gRPC.

This is very much a work in progress!

## Getting Started

The project requires python 3.6+ 

The `protocol` folder contains the definitions of the gRPC services. The `compile.sh` script can be used to generate
the python bindings.

The `python-libraries` folder contains the library to write narupa clients and
servers in python, as well as the services implemented in python. The
`python-libraries/prototypes` directory contains examples and (sometimes
unmaintained) prototypes using the python libraries.

The `csharp-libraries/Narupa.Protocol` folder contains C# implementations of clients for receiving trajectories and structures.

## Coming Soon

* Interactive MD servers
* Multiplayer
* Lobby
* Trajectory serving

## Running the tests

Running the tests is a crucial part of keeping the code base functional. To run the test of the python libraries, run:

    python -m pytest python-libraries

