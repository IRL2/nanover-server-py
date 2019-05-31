"""
Serialize and deserialize OpenMM simulations to and from XML files.

A simulation is described as the concatenation of a starting structure as a PDB
file, an OpenMM serialized system, an OpenMM serialized integrator, and,
optionally, an OpenMM serialized state. The resulting XML file looks like:

::

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

The ``System`` and ``Integrator`` tags are the roots of the serialized system
and integrator, respectively.

This module provides a function :func:`serialize_simulation` that generates an
XML file from an existing instance of :class:`simtk.openmm.app.Simulation`, and
a function :func:`deserialize_simulation` that creates an instance of simulation
from an XML file.
"""

from io import StringIO
from tempfile import TemporaryDirectory
import os
from xml.dom.minidom import getDOMImplementation, parseString

from simtk.openmm import app, XmlSerializer

ROOT_TAG = 'OpenMMSimulation'


def serialize_simulation(simulation: app.Simulation) -> str:
    """
    Generate an XML string from a simulation.

    :param simulation: The simulation to serialize.
    :return: A string with the content of an XML file describing the simulation.
    """
    implementation = getDOMImplementation()
    document = implementation.createDocument(None, ROOT_TAG, None)

    # Extract the PDB
    positions = simulation.context.getState(getPositions=True).getPositions()
    pdb_content = StringIO()
    app.PDBFile.writeFile(simulation.topology, positions, pdb_content)
    pdb_node = document.createElement('pdb')
    pdb_node.appendChild(document.createTextNode(pdb_content.getvalue()))

    # Extract the system
    system_xml_str = XmlSerializer.serialize(simulation.system)
    system_document = parseString(system_xml_str)

    # Extract the integrator
    integrator_xml_str = XmlSerializer.serialize(simulation.integrator)
    integrator_document = parseString(integrator_xml_str)

    # Combine the element in a single
    root = document.documentElement
    root.appendChild(pdb_node)
    root.appendChild(system_document.documentElement)
    root.appendChild(integrator_document.documentElement)

    return root.toprettyxml()


def deserialize_simulation(xml_content: str) -> app.Simulation:
    """
    Create an OpenMM simulation from XML.

    :param xml_content: The content of an XML file as a string.
    :return: An instance of the simulation.
    """
    document = parseString(xml_content)

    pdb_node = _get_node_and_raise_if_more_than_one(document, 'pdb')
    pdb_content = StringIO(pdb_node.firstChild.nodeValue)
    with TemporaryDirectory() as tmp_dir:
        pdb_path = os.path.join(tmp_dir, 'configuration.pdb')
        with open(str(pdb_path), 'w') as outfile:
            outfile.write(pdb_content.getvalue())
        pdb = app.PDBFile(str(pdb_path))

    system_node = _get_node_and_raise_if_more_than_one(document, 'System')
    system_content = system_node.toprettyxml()
    system = XmlSerializer.deserialize(system_content)

    integrator_node = _get_node_and_raise_if_more_than_one(document, 'Integrator')
    integrator_content = integrator_node.toprettyxml()
    integrator = XmlSerializer.deserialize(integrator_content)

    simulation = app.Simulation(
        topology=pdb.topology,
        system=system,
        integrator=integrator,
    )
    simulation.context.setPositions(pdb.positions)
    return simulation


def _get_node_and_raise_if_more_than_one(document, tag_name: str):
    nodes = document.getElementsByTagName(tag_name)
    if not nodes:
        raise IOError('No {} tag defined in the XML.'.format(tag_name))
    if len(nodes) != 1:
        raise IOError('More than one {} tag defined in the XML.'.format(tag_name))
    return nodes[0]
