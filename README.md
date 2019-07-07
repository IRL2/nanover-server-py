# Narupa 2 Protocol

Repository containing the gRPC protocol and python based implementations 
of servers for Narupa 2. 

## Getting Started

The `protocol` folder contains the definitions of the gRPC services. 

The `python-libraries` folder contains the library to write Narupa clients and
servers in python, as well as the services implemented in python. The
`python-libraries/prototypes` directory contains examples and (sometimes
unmaintained) prototypes using the python libraries.

The `csharp-libraries/Narupa.Protocol` folder contains C# implementations of clients for receiving trajectories and structures.

### Setup narupa-protocol on Mac and Linux

* Install Anaconda (avoid Anaconda 2.7 as it is outdated)
* Install dotnet
* Clone the narupa-protocol repository
* In a terminal, in the repository root:
    * Create a conda environment (here we call the environment "narupa-dev"): `conda create -n narupa-dev python>3.6`
    * Activate the conda environment: `conda activate narupa-dev`
    * Install the required conda package: `conda install -c omnia -c conda-forge openmm MDAnalysis MDAnalysisTests ase`
    * Compile the protocol and install the Narupa libraries in your conda environment: `./compile.sh`. If you plan on
      modifying the python packages, run `./compile.sh --edit` instead.

### Setup narupa-protocol on Windows

* Install Anaconda (avoid Anaconda 2.7 as it is outdated)
* Install the .NET core SDK (see <https://dotnet.microsoft.com/download>)
* Clone the narupa-protocol repository
* In the "Anaconda Powershell Prompt":
    * Create a conda environment (here we call the environment "narupa-dev"): `conda create -n narupa-dev python>3.6`
    * Activate the conda environment: `conda activate narupa-dev`
    * Install the required conda packages: `conda install -c omnia -c conda-forge openmm MDAnalysis MDAnalysisTests ase`
    * Compile the protocol and install the Narupa libraries in your conda environment: `./win_compile.ps1`. If you plan on modifying the python packages, run `./win_compile.ps1 -edit` instead.

## Coming Soon

* Multiplayer
* Lobby

## Running the tests

Running the tests is a crucial part of keeping the code base functional. To run the test of the python libraries, run:

    python -m pytest python-libraries

## Running Example Servers

`narupa.mdanalysis` provides a server for the trajectory service that infinitely loops over the frames of an example
trajectory. To serve the frames on port 54321, run

    python .\python-libraries\narupa-mdanalysis\examples\example.py

`narupa.ase` provides a server that runs a demo ASE simulation using the OpenMM
forcefield. An example usage is:

    python python-libraries/narupa-ase/examples/imd_openmm.py python-libraries/narupa-ase/examples/nanotube.xml 

