# Narupa 2 Protocol

Repository containing the gRPC protocol and python based implementations 
of servers for Narupa 2. 

This repository is maintained by the Intangible Realities Laboratory, University Of Bristol, 
and distributed under [GPLv3](LICENSE).
See [the list of contributors](CONTRIBUTORS.md) for the individual authors of the project.

## Installing narupa-server with Anaconda

* Install Anaconda (avoid Anaconda 2.7 as it is outdated)
* Create a conda environment (here we call the environment "narupa"): `conda create -n narupa "python>3.6"`
* Activate the conda environment: `conda activate narupa`
* Install the Narupa packages: 

`conda install -c irl -c omnia -c conda-forge narupa-server`

Developers will want the manual install described below.

## Getting Started

### Quick Start

`narupa.ase` provides a command line interface for running OpenMM simulations. For example, from the `narupa-protocol` directory:

    narupa-omm-ase examples/ase/openmm_files/nanotube.xml 

### Tutorials

The [examples](examples) folder contains example notebooks for using Narupa. They
are organised into the following folders: 

* [ase](examples/ase) - Get up and running with interactive simulations with ASE and OpenMM. 
   - [Basic Example](examples/ase/basic_example.ipynb) - Toy example of an interactive simulation. 
   - [Nanotube](examples/ase/openmm_nanotube.ipynb) - Set up an interactive nanotube simulation with OpenMM.
   - [Neuraminidase](examples/ase/openmm_neuraminidase.ipynb) - Set up a ligand-protein binding simulation with OpenMM,
   and experiment with Narupa visualizations. 
   - [Graphene](examples/ase/openmm_graphene.ipynb) - Set up a graphene simulation with physics parameters
   that can be adjusted on the fly. 
* [mdanalysis](examples/mdanalysis) - Visualize static structures and trajectories with MDAnalysis and Narupa. 
    - [Structure](examples/mdanalysis/mdanalysis_lsd.ipynb) - Visualize LSD bound to a receptor in Narupa.
    - [Trajectory](examples/mdanalysis/mdanalysis_trajectory.ipynb) - Build your own trajectory viewer with MDAnalysis
    and Narupa.
* [fundamentals](examples/fundamentals) - Understand how Narupa works, so you can create your own applications.

### Exploring the code  

The `protocol` folder contains the definitions of the gRPC services. 

The `python-libraries` folder contains the library to write Narupa clients and
servers in python, as well as the services implemented in python. The
`python-libraries/prototypes` directory contains examples and (sometimes
unmaintained) prototypes using the python libraries.

The `csharp-libraries/Narupa.Protocol` folder contains C# implementations of clients for receiving trajectories and structures.

## Setting up for development 

The instructions below describe how to set up an development environment with narupa protocol for Windows, Linux 
and Mac.

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
    
Optionally, you can run most of tests in parallel with pytest-xdist:
    
    pytest -m pip install pytest-xdist
    python -m pytest python-libraries -n auto -m 'not serial'
    python -m pytest python-libraries -n0 -m 'serial'

## Citation and External Libraries

If you find this project useful, please cite the following paper: 

M. O’Connor, S.J. Bennie, H.M. Deeks, A. Jamieson-Binnie, A.J. Jones, R.J. Shannon, R. Walters, T. Mitchell, A.J. Mulholland, D.R. Glowacki, [“Interactive molecular dynamics from quantum chemistry to drug binding: an open-source multi-person virtual reality framework”](https://aip.scitation.org/doi/10.1063/1.5092590), J. Chem Phys 150, 224703 (2019)

This project has been made possible by the following open source projects. We gratefully thank them for their efforts, and suggest that you use and cite them:

* [gRPC](https://grpc.io/) (Apache v2) - Communication protocol.
* [ASE](https://wiki.fysik.dtu.dk/ase/) (LGPLv3): Atomic simulation environment used for running simulations ([citation](https://iopscience.iop.org/article/10.1088/1361-648X/aa680e)).
* [OpenMM](http://openmm.org/) (MIT, LGPLv3): GPU accelerated molecular mechanics library ([citation](https://simtk.org/plugins/publications/index.php/?group_id=161)).
* [LAMMPS](https://lammps.sandia.gov/) (GPLv2): Molecular mechanics library ([citation](https://lammps.sandia.gov/cite.html)).
* [MDAnalysis](https://www.mdanalysis.org/) (GPLv2): Molecular dynamics analysis library ([citations](https://www.mdanalysis.org/pages/citations/)).
* [python-osc](https://pypi.org/project/python-osc/) (Public domain) - Open sound control library.
* [MPI4Py](https://mpi4py.readthedocs.io/en/stable/index.html) ([BSD 2-clause license](https://bitbucket.org/mpi4py/mpi4py/src/master/LICENSE.rst)): MPI library for python, used with LAMMPS ([citation](https://mpi4py.readthedocs.io/en/stable/citing.html)).
* [Numpy](https://numpy.org/) (BSD) - Numerical computation library.
* [Netifaces](https://pypi.org/project/netifaces/) (MIT) - Portable library for accessing network interface information.
* [Pytest](https://docs.pytest.org/en/latest/) (MIT) - Python testing framework
* [Hypothesis](https://hypothesis.readthedocs.io/en/latest/) ([Mozilla Public License 2.0](https://github.com/HypothesisWorks/hypothesis/blob/master/hypothesis-python/LICENSE.txt)) - Python testing framework.