[nanover-docs]: https://irl2.github.io/nanover-docs
[nanover-server-pypi]: https://pypi.org/project/nanover-server/
[nanover-server-conda]: https://anaconda.org/channels/irl/packages/nanover-server/overview
[imd-xr-app]: https://www.meta.com/en-gb/experiences/nanover-imd-xr/33606061302340842/
[imd-xr-repo]: https://github.com/IRL2/nanover-imd-xr

# NanoVer Python Server

[![License: MIT](https://img.shields.io/badge/License-MIT-darkblue.svg)](LICENSE)
[![Docs](https://img.shields.io/badge/Docs-latest-blue.svg)](https://irl2.github.io/nanover-docs)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.08118/status.svg)](https://doi.org/10.21105/joss.08118)
[![Build Status](https://github.com/IRL2/nanover-server-py/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/IRL2/nanover-server-py/actions/workflows/tests.yml?query=branch%3Amain)

Repository containing the python based implementations of servers for NanoVer, providing a framework for working with
interactive molecular dynamics simulations for computational chemistry research.

This software is designed to be used with NanoVer's iMD-XR client: [Meta Quest 3 app](imd-xr-app) ([source repo][imd-xr-repo]).

This repository is maintained by the Intangible Realities Laboratory, University of Santiago de Compostela, and is
distributed under the [MIT](LICENSE) license. See [the list of contributors](CONTRIBUTORS.md) for the individual authors
of the project. If you would like to contribute to NanoVer, please see our [contributing file](./CONTRIBUTING.md) for
guidelines on how to do this.

For more information please take a look at [the project's documentation](https://irl2.github.io/nanover-docs).

## Contents

1. [Getting started](#Getting-started)
2. [Installing the nanover-server package](#installing-the-nanover-server-package)
3. [Running the command line servers](#running-the-command-line-servers)
4. [Tutorials](#tutorials)
5. [Troubleshooting](#troubleshooting)
6. [Developer installation](#Developer-installation)
7. [Running the tests](#Running-the-tests)
8. [Citation and external libraries](#Citation-and-external-libraries)

----

## Getting started

The quickest way to see NanoVer in action is to install the NanoVer iMD-XR client app on a Meta Quest 3 headset, install
the `nanover-server` package in a python environment, and then run one of the package's command line servers with one of
the example files in this repo:

1. [Install the Nanover iMD-XR app][imd-xr-app]
2. [Install the nanover-server package](#installing-the-nanover-server-package)
3. [Run a command line server](#running-the-command-line-servers)

## Installing the nanover-server package

NanoVer is available as a [PyPI package][nanover-server-pypi] and as an [Anaconda package][nanover-server-conda]. The
two packages are identical but when using conda you should prefer the Anaconda package to ensure that dependencies are
installed in the manner most compatible with your conda environment.

Check out the [Installation & Getting Started](https://irl2.github.io/nanover-docs/installation) page in our
documentation for more complete instructions on installing NanoVer for those not familiar with conda or python.

### Installing the PyPI package

```shell
# with pip
pip install nanover-server
```

```shell
# with uv
uv pip install nanover-server
```

### Installing the Anaconda package

```shell
# in a fresh environment
conda create -n nanover -c conda-forge irl::nanover-server
```

```shell
# in an existing environment
conda install -c conda-forge irl::nanover-server
```

## Running the command line servers

`nanover` provides a command line interface for running OpenMM simulations. For example, from the `nanover-server-py`
directory:

```shell
# live openmm simulation
nanover-server --omm ./notebooks/systems/nanotube.xml
```

```shell
# static pdb structure
nanover-server --mda ./notebooks/systems/17-ala.pdb
```

```shell
# nanover format recording
nanover-server --playback ./notebooks/systems/nanotube-example-recording.nanover.zip
```

Learn more about running a NanoVer server
[here in our documentation](https://irl2.github.io/nanover-docs/tutorials/basics.html#running-a-server).

## Tutorials

The [tutorials](tutorials) folder contains [Jupyter notebooks](https://jupyter.org/) that demostrate how to get started
NanoVer. Please head to the [Tutorials page](https://irl2.github.io/nanover-docs/tutorials/tutorials.html) of the
[project's documentation][nanover-docs] for more information!

## Troubleshooting

### Connecting to servers

If you are having trouble browsing or auto-connecting to servers, you can run `nanover-essd-list` to verify which local
network servers are visible to your machine. Even if you can't see your server in the server listing, you may be able
to connect to it directly by finding the machine's ip in the local network and the server's port. Connectivity problems
are likely when using bare eduroam, and it is recommended instead to use a private local network.

## Developer installation

To install for development, first follow the [user installation instructions](#installing-the-nanover-server-package)
for either the PyPI or the Anaconda package, and then these additional steps to reinstall the package from the source
code in [development mode](https://setuptools.pypa.io/en/latest/userguide/development_mode.html) and install additional
development dependencies.

### Clone the nanover-server-py repository

```shell
git clone https://github.com/IRL2/nanover-server-py.git
```

### Enter the nanover-server-py directory

```shell
cd nanover-server-py
```

### Install the package from source code

```shell
# with pip
pip install -e .[dev]
```

```shell
# with uv
uv pip install -e .[dev]
```

**NOTE:** You should repeat this step whenever the dependencies or available command line apps are changed. If in doubt
after pulling updates, rerun this installation step.

## Running the tests locally

All code changes have to pass [a series of automatic tests](./.github/workflows/tests.yml) that attempt to verify code
quality and continued functionality of the project. You can run these locally to verify your changes in advance.

### Running the unit tests

The unit tests check code functionality of the python libraries. To run them:

```shell
pytest tests
```

Optionally, you can run most of the tests in parallel with `pytest-xdist`:

```shell
pip install pytest-xdist
pytest tests -n auto -m 'not serial'
pytest tests -n0 -m 'serial'
```

### Running the linting and reformatting

The linting checks code style using ruff:

```shell
ruff check src
ruff format --check src
```

Ruff can also automatically fix and reformat the files:

```shell
ruff check --fix src
ruff format src
```

### Running the type checker

The type checker looks at the type hints in the code to make sure they are consistent and help find potential errors:

```shell
mypy src
```

## Citation and external libraries

Any work that uses NanoVer should cite the following publications:

> Stroud, H. J., Wonnacott, M. D., Barnoud, J., Roebuck Williams, R., Dhouioui, M., McSloy, A., Aisa, L., Toledo, L. E.,
> Bates, P., Mulholland, A. J., & Glowacki, D. R. (2025). NanoVer Server: A Python Package for Serving Real-Time
> Multi-User Interactive Molecular Dynamics in Virtual Reality. *Journal of Open Source Software*, *10* (110), 8118. https://doi.org/10.21105/joss.08118

> Jamieson-Binnie, A. D., O’Connor, M. B., Barnoud, J., Wonnacott, M. D., Bennie, S. J., & Glowacki, D. R. (2020, August
> 17). Narupa iMD: A VR-Enabled Multiplayer Framework for Streaming Interactive Molecular Simulations. ACM SIGGRAPH 2020
> Immersive Pavilion. SIGGRAPH ’20: Special Interest Group on Computer Graphics and Interactive Techniques
> Conference. https://doi.org/10.1145/3388536.3407891

> O’Connor, M., Bennie, S. J., Deeks, H. M., Jamieson-Binnie, A., Jones, A. J., Shannon, R. J., Walters, R., Mitchell,
> T., Mulholland, A. J., & Glowacki, D. R. (2019). Interactive molecular dynamics from quantum chemistry to drug
> binding: an open-source multi-person virtual reality framework, *The Journal of Chemical Physics*, *150* (22), 224703. https://doi.org/10.1021/acs.jcim.0c01030

This project has been made possible by the following open source projects. We gratefully thank them for their efforts,
and suggest that you use and cite them:

* [ASE](https://wiki.fysik.dtu.dk/ase/) (LGPLv3): Atomic simulation environment used for running simulations
  ([citation](https://iopscience.iop.org/article/10.1088/1361-648X/aa680e)).
* [OpenMM](http://openmm.org/) (MIT, LGPLv3): GPU accelerated molecular mechanics library
  ([citation](https://simtk.org/plugins/publications/index.php/?group_id=161)).
* [MDAnalysis](https://www.mdanalysis.org/) (GPLv2): Molecular dynamics analysis library
  ([citations](https://www.mdanalysis.org/pages/citations/)).
* [NGLView](https://nglviewer.org/#nglview) (MIT): IPython/Jupyter widget to interactively view structures and
  trajectories ([citations](http://nglviewer.org/nglview/latest/#cite)).
* [NumPy](https://numpy.org/) (BSD) - Numerical computation library.
* [psutil](https://pypi.org/project/netifaces/) (BSD) - Cross-platform lib for process and system monitoring in Python.
* [pytest](https://docs.pytest.org/en/latest/) (MIT) - Python testing framework
* [Hypothesis](https://hypothesis.readthedocs.io/en/latest/)
  ([Mozilla Public License 2.0](https://github.com/HypothesisWorks/hypothesis/blob/master/hypothesis-python/LICENSE.txt)) -
  Python testing framework.
