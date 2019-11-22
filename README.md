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

### Setup narupa-protocol with Anaconda

* Install Anaconda (avoid Anaconda 2.7 as it is outdated)
* Create a conda environment (here we call the environment "narupa"): `conda create -n narupa "python>3.6"`
* Activate the conda environment: `conda activate narupa`
* Install the Narupa packages: `conda install -c irl -c omnia -c conda-forge narupa-server`

Developers will want the manual install described below.

### Setup narupa-protocol for developers on Mac and Linux

* Install Anaconda (avoid Anaconda 2.7 as it is outdated)
* Install dotnet
* Clone the narupa-protocol repository
* In a terminal, in the repository root:
    * Create a conda environment (here we call the environment "narupa-dev"): `conda create -n narupa-dev "python>3.6"`
    * Activate the conda environment: `conda activate narupa-dev`
    * Install the required conda package: `conda install -c omnia -c conda-forge openmm MDAnalysis MDAnalysisTests ase`
    * Compile the protocol and install the Narupa libraries in your conda environment: `./compile.sh`.  If you do not plan on modifying the python packages, you may run `./compile.sh --no-edit` instead. Otherwise, by default, the narupa packages will be installed in edit mode (`pip install -e`) meaning that changes in the `narupa-protocol` directory will be directly reflected in your python environment.

### Setup narupa-protocol for developers on Windows

* Install Anaconda (avoid Anaconda 2.7 as it is outdated)
* Install the .NET core SDK (see <https://dotnet.microsoft.com/download>)
* Clone the narupa-protocol repository
* In the "Anaconda Powershell Prompt":
    * Create a conda environment (here we call the environment "narupa-dev"): `conda create -n narupa-dev "python>3.6"`
    * Activate the conda environment: `conda activate narupa-dev`
    * Install the required conda packages: `conda install -c omnia -c conda-forge openmm MDAnalysis MDAnalysisTests ase`
    * Compile the protocol and install the Narupa libraries in your conda environment: `./win_compile.ps1`.  If you do not plan on modifying the python packages, run `./win_compile.ps1 -noedit` instead. Otherwise, by default, the narupa packages will be installed in edit mode (`pip install -e`) meaning that changes in the `narupa-protocol` directory will be directly reflected in your python environment.

## Running the tests

Running the tests is a crucial part of keeping the code base functional. To run the test of the python libraries, run:

    python -m pytest python-libraries

## Running the examples

### ASE IMD Simulations 

`narupa.ase` provides a command line interface for running serialised OpenMM simulations. For example, from the `narupa-protocol` directory:

    narupa-omm-ase python-libraries/narupa-ase/examples/nanotube.xml 

The example files are distributed in the directory
`python-library/narupa-ase/examples` from the [git repository](https://gitlab.com/intangiblerealities/narupa-protocol/tree/master/python-libraries/narupa-ase/examples).
The examples include:

* [A carbon nanotube and a methane molecule](https://gitlab.com/intangiblerealities/narupa-protocol/raw/master/python-libraries/narupa-ase/examples/nanotube.xml)
* [A helicene molecule](https://gitlab.com/intangiblerealities/narupa-protocol/raw/master/python-libraries/narupa-ase/examples/helicene.xml)
* [Neuraminidase and tamiflu](https://gitlab.com/intangiblerealities/narupa-protocol/raw/master/python-libraries/narupa-ase/examples/neuraminidase.xml)

#### Jupyter Notebooks 

The [`python-libraries/narupa-ase/examples`](https://gitlab.com/intangiblerealities/narupa-protocol/tree/master/python-libraries/narupa-ase/examples) examples folder also contains several
Jupyter notebooks that demonstrate visualisation and interaction from a notebook.
The [Narupa ASE documentation](python-libraries/narupa-ase/README.md) provides more details on setting up ASE simulations.

### MD Analysis Trajectories

`narupa.mdanalysis` provides a server for the trajectory service that infinitely loops over the frames of an example
trajectory. To serve the frames on port 54321, from the `narupa-protocol` directory, run

    python ./python-libraries/narupa-mdanalysis/examples/example.py