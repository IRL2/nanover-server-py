OpenMM server for NanoVer
========================

This server sends over the network the frames of a running OpenMM simulation. The server can
be started from two different interfaces: the python interpreter, or the command line interface.
The python interface is the most potent; it gives direct access to OpenMM and allows to use third
party tools around the simulation engine to setup or run simulations. The command line interface
provides a convenient and portable way to run and serve a simulation that has previously been set up.

Running a server from python
----------------------------

When used as a python library, `nanover-openmm` provides the
`nanover.openmm.Server` class. This class can be used as follow:

```python
import openmm as mm
from openmm import app

from nanover.openmm import Server

# Build an OpenMM simulation
simulation = app.Simulation(...)

# Setup the NanoVer server
server = Server(simulation, address='localhost', port=8000)
# Make the server verbose. The server is quiet by default.
# When verbose, the server prints the number of steps run,
# the potential energy of the system in kJ/mol, and the
# performances in ns/day on the standard output.
server.verbose = True
# Run the server for an infinite number of steps
server.run()
# Press Ctrl+C to interrupt the server.

# Close the network connection.
server.close()
```

The [OpenMM documentation](http://docs.openmm.org/latest/userguide/application.html#running-simulations)
describes different ways of setting up a simulation, including from Amber or Gromacs input files.

The NanoVer OpenMM package provides a way of serializing an entire simulation as a single XML file. 
This allows one to decouple the preparation and the
running of a simulation. See the "[Writing a XML simulation file](#writing-a-xml-simulation-file)"
section to learn more about the file format to describe a simulation, and
how to produce such a file.

A server can be created from the XML description of a simulation:

```python
server = Server.from_xml_input('simulation.xml')
```

In any case, the OpenMM simulation object is exposed by the server object through the
`simulation` attribute of the server object so every operation usually performed on OpenMM simulation
object (such as adding a reporter, or accessing the context) can be achieved.
Refer to the OpenMM [documentation](http://openmm.org/documentation.html) for more details.

A NanoVer server can also be included in an existing OpenMM workflow by adding
a `NanoverReporter` to an existing simulation object:

```python
import openmm as mm
from openmm import app

from nanover.trajectory import FrameServer
from nanover.openmm import NanoverReporter

# Create an OpenMM simulation
simulation = app.Simulation(...)

# Start a NanoVer frame server
frame_server = FrameServer(address='localhost', port=8000)
# Setup a reporter to link OpenMM and NanoVer
nanover_reporter = NanoverReporter(report_interval=1, frame_server=frame_server)
# Add the NanoVer reporter to the simulation
simulation.reporters.append(nanover_reporter)

# Run the simulation for 1000 steps
simulation.step(1000)

# Close the network connection.
frame_server.close()
```

Running a server from the command line
--------------------------------------

When `nanover-opemm` is installed, it provides the `nanover-omm-server`
command in the command line. When provided with the description of an
OpenMM simulation as an XML file, `nanover-omm-server` runs the simulation
and serves the frame for NanoVer. The host address and port can be set with
the `--address` and the `--port` option, respectively.

See the "[Writing a XML simulation file](#writing-a-xml-simulation-file)"
section to learn more about the file format to describe a simulation, and
how to produce such a file.

Run `nanover-omm-server --help` in a terminal to read more about the options.

The XML file for the neuraminidase demo can be found in the `examples` directory.

Writing a XML simulation file
-----------------------------

`nanover-openmm` can use an XML file to describe a simulation. Such an XML file is
the concatenation of a PDB file (used to get the initial coordinates and the
topology of the system), an XML serialized OpenMM system (used to get the forces
to apply), and an XML serialized OpenMM integrator. The file looks like:

```xml
<OpenMMSimulation>
    <pdb>
        // pasted content of the PDB file
    </pdb>
    <System ...>
        // XML content of the OpenMM serialized system
    </System>
    <Integrator ...>
        // XML content of the OpenMM serialized integrator
    </Integrator>
</OpenMMSimulation>
```

XML description of simulations can be written with the
`nanover.openmm.serializer.serialize_simulation` function.
