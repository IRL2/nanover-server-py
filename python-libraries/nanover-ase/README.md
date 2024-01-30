ASE server for NanoVer
========================

This server implements interactive molecular dynamics (IMD) for an ASE molecular dynamics simulation. 

Running an ASE OpenMM server from the command line
-----------------------------------------------

When `nanover-ase` is installed, it provides the `nanover-omm-ase`
command in the command line. When provided with the description of an
OpenMM simulation as an XML file serialised as described in the [NanoVer OpenMM documentation](../nanover-openmm/README.md) 
, `nanover-omm-ase` runs an interactive simulation. 
The host address and port can be set with
the `--address` and the `--port` option, respectively.


Running a server from python
----------------------------

The `nanover-ase` module provides the
`nanover.ase.NanoVerASEDynamics` class. Given an ASE simulation set up with an 
[ASE molecular dynamics runner](https://wiki.fysik.dtu.dk/ase/ase/md.html), this class will 
attach interactive molecular dynamics functionality and frame serving to the dynamics. 
An example is given below, assuming an ASE Atoms object has been set up, named `atoms`:

```python
from ase import units
from ase.md import Langevin
from nanover.ase.imd import NanoVerASEDynamics

# Given some ASE atoms object appropriately set up, set up dynamics.
dyn = Langevin(atoms, 1 * units.fs, temperature_K=300, fraction=0.1)

# Set up a basic NanoVer server to run the interactive dynamics.
with NanoVerASEDynamics.basic_imd(dyn) as imd:
    while True:
        imd.run(100)
```

Full examples are given in the [ASE examples](../../examples/ase) folder, which additionally
contains several Jupyter notebooks that explore how NanoVer can be used with OpenMM:

* `nanover_ase_client_server`: A notebook showing how one can run the server for an OpenMM simulation, 
connect a client to it, and render a simple visualisation. 
* `nanover_ase_interactive_md`: A notebook that runs a simulation of a carbon nanotube, then applies
interactive forces to it from the notebook.
* `nanover_interactive_visualiser`: A notebook that assumes a server is already running, and visualises it
with [NGLView](https://github.com/arose/nglview). To run this notebook, ensure NGLView is installed with:

```bash
conda install nglview -c conda-forge
# might need: jupyter-nbextension enable nglview --py --sys-prefix

# if you already installed nglview, you can `upgrade`
conda upgrade nglview --force
# might need: jupyter-nbextension enable nglview --py --sys-prefix
```


