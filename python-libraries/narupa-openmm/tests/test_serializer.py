# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Tests for :mod:`narupa.openmm.serializer`.
"""

# Pylint does not recognize pytest fixtures which creates fake warnings.
# pylint: disable=redefined-outer-name,unused-import
from xml.dom.minidom import parseString
import pytest

from narupa.openmm.serializer import (
    serialize_simulation,
    deserialize_simulation,
    ROOT_TAG,
)

from .simulation_utils import basic_simulation


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


def test_serialize_runs(basic_simulation_xml):
    """
    Test that a serialized simulation can be deserialized and run.
    """
    simulation = deserialize_simulation(basic_simulation_xml)
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
