ASE server for Narupa
========================

This server implements interactive molecular dynamics (IMD) for an ASE molecular dynamics simulation. 

Running a server from python
----------------------------

The `narupa-ase` module provides the
`narupa.ase.ASEImdServer` class. Given an ASE simulation set up with an 
[ASE molecular dynamics runner](https://wiki.fysik.dtu.dk/ase/ase/md.html), this class will 
attach interactive molecular dynamics functionality and frame serving to the dynamics. 

```python
from ase import units
from ase.md import Langevin
from narupa.ase.imd_server import IMDServer

# Given some ASE atoms object appropriately set up, set up dynamics.
dyn = Langevin(atoms, 1 * units.fs, 300, 0.1)

# Attach the IMD calculator and server to the dynamics object. 
imd = IMDServer(dyn)
while True:
    imd.run(100)
```

Full examples are given in the [examples](./examples) folder.

