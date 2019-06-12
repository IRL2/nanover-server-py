ASE server for Narupa
========================

This server implements interactive molecular dynamics (IMD) for an ASE molecular dynamics simulation. 

Running an OpenMM server from the command line
-----------------------------------------------

When `narupa-ase` is installed, it provides the `narupa-omm-ase`
command in the command line. When provided with the description of an
OpenMM simulation as an XML file, `narupa-omm-ase` runs an interactive simulation. 
The host address and port can be set with
the `--address` and the `--port` option, respectively.


Running a server from python
----------------------------

The `narupa-ase` module provides the
`narupa.ase.ASEImdServer` class. Given an ASE simulation set up with an 
[ASE molecular dynamics runner](https://wiki.fysik.dtu.dk/ase/ase/md.html), this class will 
attach interactive molecular dynamics functionality and frame serving to the dynamics. 

```python
from ase import units
from ase.md import Langevin
from narupa.ase.imd_server import ASEImdServer

# Given some ASE atoms object appropriately set up, set up dynamics.
dyn = Langevin(atoms, 1 * units.fs, 300, 0.1)

# Attach the IMD calculator and server to the dynamics object. 
imd = ASEImdServer(dyn)
while True:
    imd.run(100)
```

Full examples are given in the [examples](./examples) folder, which contains several
Jupyter notebooks that explore how Narupa can be used with OpenMM.

