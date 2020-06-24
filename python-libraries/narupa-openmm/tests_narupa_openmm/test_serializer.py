# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Tests for :mod:`narupa.openmm.serializer`.
"""

# Pylint does not recognize pytest fixtures which creates fake warnings.
# pylint: disable=redefined-outer-name,unused-import
from xml.dom.minidom import parseString
import pytest

import simtk.openmm as mm

from narupa.openmm.serializer import (
    serialize_simulation,
    deserialize_simulation,
    ROOT_TAG,
)
from narupa.openmm.imd import create_imd_force

from .simulation_utils import basic_simulation, empty_imd_force


@pytest.fixture
def basic_simulation_xml(basic_simulation):
    """
    Generate a XML serialized simulation from the basic test simulation.
    """
    return serialize_simulation(basic_simulation)


def remove_xml_tag(simulation_xml: str, tag_to_remove: str) -> str:
    """
    Remove an branch of an XML file.

    This is meant to break valid XML serialized simulation by removing a given
    required tag.

    :param simulation_xml: XML content as a string.
    :param tag_to_remove: The name of a tag to remove from the XML. All
        occurrences of that tag will be removed, so will be the children of
        these tags.
    :return: The new XML content as a string.
    """
    document = parseString(simulation_xml)
    for node in document.getElementsByTagName(tag_to_remove):
        node.parentNode.removeChild(node)
    return document.toprettyxml()


def add_extra_xml_tag(simulation_xml: str, extra_tag: str) -> str:
    """
    Add an empty tag to a XML file under the root.
    """
    document = parseString(simulation_xml)
    extra_node = document.createElement(extra_tag)
    root = document.getElementsByTagName(ROOT_TAG)[0]
    root.appendChild(extra_node)
    return document.toprettyxml()


@pytest.mark.parametrize('imd_force', (None, create_imd_force()))
def test_serialize_runs(basic_simulation_xml, imd_force):
    """
    Test that a serialized simulation can be deserialized and run.
    """
    simulation = deserialize_simulation(basic_simulation_xml, imd_force=imd_force)
    assert simulation.currentStep == 0
    simulation.step(5)
    assert simulation.currentStep == 5


@pytest.mark.parametrize('section_missing', ('pdb', 'System', 'Integrator'))
def test_missing_section(basic_simulation_xml, section_missing):
    """
    Make sure that the deserialization raise an exception when a base section
    is missing.
    """
    broken_xml = remove_xml_tag(basic_simulation_xml, section_missing)
    with pytest.raises(IOError):
        deserialize_simulation(broken_xml)


@pytest.mark.parametrize('section_to_duplicate', ('pdb', 'System', 'Integrator'))
def test_duplicate_section(basic_simulation_xml, section_to_duplicate):
    """
    Make sure de deserialization fails if a base section appears more than once.
    """
    broken_xml = add_extra_xml_tag(basic_simulation_xml, section_to_duplicate)
    with pytest.raises(IOError):
        deserialize_simulation(broken_xml)


def test_imd_force(basic_simulation_xml, empty_imd_force):
    """
    When deserializing a simulation and passing an imd force, the force is
    added to the system.
    """
    simulation = deserialize_simulation(basic_simulation_xml, empty_imd_force)

    system = simulation.system
    all_forces = system.getForces()
    # The wrapper for the C++ force is recreated at some point, so what
    # system.getForces returns is a different python object as what we
    # inputted initially, even though it references the same force under the
    # hood. As I do not know how to access the common reference between the
    # wrapper objects for the same force, I have to do some gymnastic.
    putative_imd_forces = [
        force for force in all_forces if (
            isinstance(force, mm.CustomExternalForce)
            and [
                force.getPerParticleParameterName(i)
                for i in range(force.getNumPerParticleParameters())
            ] == ['fx', 'fy', 'fz']
        )
    ]
    # In principle, there could be more than one force that matches the above
    # description. In practice, though, it is not the case for now.
    assert len(putative_imd_forces) == 1
    force_obtained = putative_imd_forces[0]
    force_added = empty_imd_force
    # The forces are the same if by modifying one we also modify the other.
    force_added.setParticleParameters(0, 0, (1.0, 2.0, 3.0))
    parameters = force_obtained.getParticleParameters(0)
    assert parameters == [0, (1.0, 2.0, 3.0)]


