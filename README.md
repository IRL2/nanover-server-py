# NanoVer Python Server

[![License: MIT](https://img.shields.io/badge/License-MIT-darkblue.svg)](LICENSE)
[![Docs](https://img.shields.io/badge/Docs-latest-blue.svg)](https://irl2.github.io/nanover-docs) 
[![DOI](https://joss.theoj.org/papers/10.21105/joss.08118/status.svg)](https://doi.org/10.21105/joss.08118) 
[![Build Status](https://github.com/IRL2/nanover-server-py/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/IRL2/nanover-server-py/actions/workflows/tests.yml?query=branch%3Amain)

Repository containing the python based implementations
of servers for NanoVer, providing a framework for developing interactive molecular dynamics simulations.

This software is designed to be used with **NanoVer VR clients**, 
e.g. [NanoVer iMD-VR](https://github.com/IRL2/nanover-imd-vr).

This repository is maintained by the Intangible Realities Laboratory, University of Santiago de Compostela,
and is distributed under the [MIT](LICENSE) license.
See [the list of contributors](CONTRIBUTORS.md) for the individual authors of the project. If you would like to contribute to NanoVer, please see our [contributing file](./CONTRIBUTING.md) for guidelines on how to do this.

For more information please take a look at [the project's documentation](https://irl2.github.io/nanover-docs).

## Contents

1. [Getting started](#Getting-started)
2. [User installation](#User-installation)
3. [Developer installation](#Developer-installation)
4. [Running the tests](#Running-the-tests)
5. [Running the tutorials](#Running-the-tutorials)
6. [Troubleshooting](#Troubleshooting)
7. [Citation and external libraries](#Citation-and-external-libraries)

----

## Getting started

Here are some quick notes to get you started with NanoVer! 
If you haven't installed NanoVer yet, please go to [User installation](#User-installation) or 
[Developer installation](#Developer-installation).

### Running a server via the command line

`nanover.omni` provides a command line interface for running OpenMM simulations. For example, from the `nanover-server-py` directory:

    nanover-omni --omm examples/ase/openmm_files/nanotube.xml

Learn more about running a NanoVer server 
[here in our documentation](https://irl2.github.io/nanover-docs/tutorials/basics.html#running-a-server).

### Tutorials

The [examples](tutorials) folder contains [Jupyter notebooks](https://jupyter.org/) that demostrate how to get started NanoVer. 
Please head to the [Tutorials page](https://irl2.github.io/nanover-docs/tutorials/tutorials.html) of the 
[project's documentation](https://irl2.github.io/nanover-docs) for more information!

### Exploring the code  

The `python-libraries` folder contains the library to write NanoVer clients and
servers in python, as well as the services implemented in python. The
`python-libraries/prototypes` directory contains examples and (sometimes
unmaintained) prototypes using the python libraries.

## User installation

Check out the [Installation & Getting Started](https://irl2.github.io/nanover-docs/installation) 
page in our documentation for detailed instructions on installing NanoVer.

### Updating the conda package

* Run `conda list nanover-server` to determine the currently installed version
* Run `conda update nanover-server` to attempt to update to latest version
* If you can't seem to update to the latest version, run `python --version` to check your python version is at least 
  as recent as in the installation instructions. 
  If it isn't you will need to create a new conda environment with a newer version of python.

## Developer installation

### Windows

* Install Anaconda
* Clone the nanover-server-py repository
* In the "Anaconda Powershell Prompt":
    * Create a conda environment (here we call the environment "nanover-dev") with the required depencies: `conda create -n nanover-dev -c conda-forge "python>3.11" openmm MDAnalysis MDAnalysisTests ase`
    * Activate the conda environment: `conda activate nanover-dev`
    * Install the NanoVer libraries in your conda environment: `./win_compile.ps1`.  If you do not plan on modifying the python packages, run `./win_compile.ps1 -noedit` instead. Otherwise, by default, the nanover packages will be installed in edit mode (`pip install -e`) meaning that changes in the `nanover-server-py` directory will be directly reflected in your python environment.

### Mac and Linux

* Install Anaconda
* Clone the nanover-server-py repository
* In a terminal, in the repository root:
    * Create a conda environment (here we call the environment "nanover-dev") with the required depencies: `conda create -n nanover-dev -c conda-forge "python>3.11" openmm MDAnalysis MDAnalysisTests ase`
    * Activate the conda environment: `conda activate nanover-dev`
    * Install the NanoVer python libraries in your conda environment: `./compile.sh`.  If you do not plan on modifying the python packages, you may run `./compile.sh --no-edit` instead. Otherwise, by default, the NanoVer packages will be installed in edit mode (`pip install -e`) meaning that changes in the `nanover-server-py` directory will be directly reflected in your python environment.

## Running the tests

All code changes have to pass a series of automatic tests ("the CI") that attempt to verify code quality and
continued functionality of the project. You can run these locally to verify your changes in advance.

### Unit Tests

The unit tests check code functionality of the python libraries. To run them:

    python -m pytest python-libraries

Optionally, you can run most of the tests in parallel with pytest-xdist:

    python -m pip install pytest-xdist
    python -m pytest python-libraries -n auto -m 'not serial'
    python -m pytest python-libraries -n0 -m 'serial'

### Formatting & Linting Tests

The formatting and linting tests check code style, and require ruff and black:

    python -m pip install ruff
    python -m pip install black
    python -m ruff check python-libraries
    python -m black --diff --check python-libraries

black can automatically reformat the files for you:

    python -m black python-libraries

### Type Checks

The type checks look at the type hints in the code to make sure they are consistent and help find potential errors:

    python -m pip install mypy
    mypy --ignore-missing-imports --namespace-packages --check-untyped-defs --allow-redefinition nanover-server
## Running the tutorials

The [tutorials](tutorials) folder contains [Jupyter notebooks](https://jupyter.org/) for examples of how to use NanoVer. 
Learn about these [Tutorials](https://irl2.github.io/nanover-docs/tutorials/tutorials.html) or
[how to run a NanoVer server](https://irl2.github.io/nanover-docs/tutorials/basics.html#running-a-server) in this
[project's documentation](https://irl2.github.io/nanover-docs).


### OpenMM IMD Simulations

`nanover.omni` provides a command line interface for running serialised OpenMM simulations. For example, from the 
`nanover-server-py` directory:

    nanover-omni --omm examples/ase/openmm_files/nanotube.xml

### ASE IMD Simulations Jupyter Notebooks

The [`examples/ase`](tutorials/ase) folder contains several Jupyter notebooks that demonstrate visualisation and interaction 
from a notebook.

### MD Analysis Trajectories

`nanover.mdanalysis` provides a server for the trajectory service that infinitely loops over the frames of an example
trajectory. To serve the frames on port 54321, from the `nanover-server-py` directory, run

    python ./examples/mdanalysis/example.py

## Troubleshooting

### Autoconnect

If you are having trouble autoconnecting to servers, you can run `nanover-essd-list` to verify which local network servers are visible to your machine.

## Old recordings

If you try to use the older .traj/.state recordings, you will find NanoVer complaining that they are not zip files.
In this case you can use the [recording converter](https://github.com/IRL2/nanover-recording-converter) in conjunction
with NanoVer to convert them to the new format.

## Citation and external libraries

Any work that uses NanoVer should cite the following publications:

> Stroud, H. J., Wonnacott, M. D., Barnoud, J., Roebuck Williams, R., Dhouioui, M., McSloy, A., Aisa, L., Toledo, L. E., Bates, P., Mulholland, A. J., & Glowacki, D. R. (2025). NanoVer Server: A Python Package for Serving Real-Time Multi-User Interactive Molecular Dynamics in Virtual Reality. *Journal of Open Source Software*, *10* (110), 8118. https://doi.org/10.21105/joss.08118

> Jamieson-Binnie, A. D., O’Connor, M. B., Barnoud, J., Wonnacott, M. D., Bennie, S. J., & Glowacki, D. R. (2020, August 17). Narupa iMD: A VR-Enabled Multiplayer Framework for Streaming Interactive Molecular Simulations. ACM SIGGRAPH 2020 Immersive Pavilion. SIGGRAPH ’20: Special Interest Group on Computer Graphics and Interactive Techniques Conference. https://doi.org/10.1145/3388536.3407891

> O’Connor, M., Bennie, S. J., Deeks, H. M., Jamieson-Binnie, A., Jones, A. J., Shannon, R. J., Walters, R., Mitchell, T., Mulholland, A. J., & Glowacki, D. R. (2019). Interactive molecular dynamics from quantum chemistry to drug binding: an open-source multi-person virtual reality framework, *The Journal of Chemical Physics*, *150* (22), 224703. https://doi.org/10.1021/acs.jcim.0c01030

This project has been made possible by the following open source projects. We gratefully thank them for their efforts, and suggest that you use and cite them:

* [ASE](https://wiki.fysik.dtu.dk/ase/) (LGPLv3): Atomic simulation environment used for running simulations ([citation](https://iopscience.iop.org/article/10.1088/1361-648X/aa680e)).
* [OpenMM](http://openmm.org/) (MIT, LGPLv3): GPU accelerated molecular mechanics library ([citation](https://simtk.org/plugins/publications/index.php/?group_id=161)).
* [MDAnalysis](https://www.mdanalysis.org/) (GPLv2): Molecular dynamics analysis library ([citations](https://www.mdanalysis.org/pages/citations/)).
* [NGLView](https://nglviewer.org/#nglview) (MIT): IPython/Jupyter widget to interactively view structures and trajectories ([citations](http://nglviewer.org/nglview/latest/#cite)).
* [python-osc](https://pypi.org/project/python-osc/) (Public domain) - Open sound control library.
* [Numpy](https://numpy.org/) (BSD) - Numerical computation library.
* [Netifaces](https://pypi.org/project/netifaces/) (MIT) - Portable library for accessing network interface information.
* [Pytest](https://docs.pytest.org/en/latest/) (MIT) - Python testing framework
* [Hypothesis](https://hypothesis.readthedocs.io/en/latest/) ([Mozilla Public License 2.0](https://github.com/HypothesisWorks/hypothesis/blob/master/hypothesis-python/LICENSE.txt)) - Python testing framework.
