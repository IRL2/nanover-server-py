# Narupa 2 Protocol

Repository containing the gRPC protocol and python based implementations
of servers for Narupa 2, providing a framework for developing interactive molecular dynamics simulations.

It is designed to be used with Narupa VR clients, e.g. [Narupa IMD](https://gitlab.com/intangiblerealities/narupa-applications/narupa-imd).

This repository is maintained by the Intangible Realities Laboratory, University Of Bristol,
and distributed under [GPLv3](LICENSE).
See [the list of contributors](CONTRIBUTORS.md) for the individual authors of the project.

**The Narupa iMD VR front-end can be found
[on its own repository](https://gitlab.com/intangiblerealities/narupa-applications/narupa-imd).**

## Getting Started

### Quick Start

`narupa.ase` provides a command line interface for running OpenMM simulations. For example, from the `narupa-protocol` directory:

    narupa-omm-ase examples/ase/openmm_files/nanotube.xml

### Tutorials

The [examples](examples) folder contains [Jupyter notebooks](https://jupyter.org/) for getting started with Narupa. They
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

The tutorials use Jupyter notebooks, [NGLView](https://github.com/arose/nglview) for visualising trajectories, and while not strictly necessary,
assumes you have the [Narupa IMD VR](https://gitlab.com/intangiblerealities/narupa-applications/narupa-imd)
application installed. These can all be installed with conda:

```bash
conda activate narupa
conda install jupyter
conda install nglview
# On Windows only:
conda install -c irl narupa-imd
```

To run the notebooks, download the repository and run jupyter (with [git](https://git-scm.com/) installed):
```bash
git clone https://gitlab.com/intangiblerealities/narupa-protocol.git
cd narupa-protocol
conda activate narupa
jupyter notebook
```


### Exploring the code  

The `protocol` folder contains the definitions of the gRPC services.

The `python-libraries` folder contains the library to write Narupa clients and
servers in python, as well as the services implemented in python. The
`python-libraries/prototypes` directory contains examples and (sometimes
unmaintained) prototypes using the python libraries.

The `csharp-libraries/Narupa.Protocol` folder contains C# implementations of clients for receiving trajectories and structures.

## Setting up for development

* Install Anaconda (avoid Anaconda 2.7 as it is outdated)
* Open the "Anaconda Powershell Prompt" to type the following commands.
* Create a conda environment (here we call the environment "narupa"): `conda create -n narupa "python>3.6"`
* Activate the conda environment: `conda activate narupa`
* Install the Narupa packages: `conda install -c irl -c omnia -c conda-forge narupa-server`

Narupa can interact with the [LAMMPS](https://lammps.sandia.gov/) simulation engine.
If you want to use this specific feature, you need to:

* install LAMMPS with python capabilities
* install mpy4py: `conda install -c conda-forge mpi4py` on Linux and MacOS,
  `python -m pip install mpi4py` on Windows.

Developers will want the manual install described below.

### Setup narupa-protocol for developers on Mac and Linux

* Install Anaconda (avoid Anaconda 2.7 as it is outdated)
* Install dotnet
* Clone the narupa-protocol repository
* In a terminal, in the repository root:
    * Create a conda environment (here we call the environment "narupa-dev"): `conda create -n narupa-dev "python>3.6"`
    * Activate the conda environment: `conda activate narupa-dev`
    * Install the required conda package: `conda install -c omnia -c conda-forge openmm MDAnalysis MDAnalysisTests ase mpi4py`
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
* The `narupa-lammps` module and its tests require MPI to be installed. Download and install Microsoft MPI from https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi

## Running the tests

Running the tests is a crucial part of keeping the code base functional. To run the test of the python libraries, run:

    python -m pytest python-libraries

Optionally, you can run most of tests in parallel with pytest-xdist:

    pytest -m pip install pytest-xdist
    python -m pytest python-libraries -n auto -m 'not serial'
    python -m pytest python-libraries -n0 -m 'serial'

## Running the examples

### ASE IMD Simulations

`narupa.ase` provides a command line interface for running serialised OpenMM simulations. For example, from the `narupa-protocol` directory:

    narupa-omm-ase examples/ase/nanotube.xml

The example files are distributed in the directory
`examples/ase/` from the [git repository](https://gitlab.com/intangiblerealities/narupa-protocol/tree/master/examples/ase).

#### Jupyter Notebooks

The [`python-libraries/narupa-ase/examples`](https://gitlab.com/intangiblerealities/narupa-protocol/tree/master/python-libraries/narupa-ase/examples) examples folder also contains several
Jupyter notebooks that demonstrate visualisation and interaction from a notebook.
The [Narupa ASE documentation](python-libraries/narupa-ase/README.md) provides more details on setting up ASE simulations.

### MD Analysis Trajectories

`narupa.mdanalysis` provides a server for the trajectory service that infinitely loops over the frames of an example
trajectory. To serve the frames on port 54321, from the `narupa-protocol` directory, run

    python ./examples/mdanalysis/example.py

## Troubleshooting

### Autoconnect

If you are having autoconnecting to servers, you can run `narupa-essd-list` to verify which local network servers are visible to your machine.

>>>>>>> master
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
