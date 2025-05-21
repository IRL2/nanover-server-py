"""
Tests for :mod:`nanover.openmm.serializer`.
"""

# Pylint does not recognize pytest fixtures which creates fake warnings.
# pylint: disable=redefined-outer-name,unused-import
from xml.dom.minidom import parseString
import pytest
from io import StringIO
from nanover.openmm.imd import create_imd_force
from nanover.openmm.serializer import (
    serialize_simulation,
    deserialize_simulation,
    ROOT_TAG,
)

from simulation_utils import (
    basic_simulation,
    basic_simulation_xml,
    build_basic_simulation,
    empty_imd_force,
)
from nanover.omni.openmm import OpenMMSimulation
from nanover.app import (
    NanoverImdClient,
    NanoverApplicationServer,
    NanoverImdApplication,
)


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


@pytest.mark.parametrize("imd_force", (None, create_imd_force()))
def test_serialize_runs(basic_simulation_xml, imd_force):
    """
    Test that a serialized simulation can be deserialized and run.
    """
    simulation = deserialize_simulation(basic_simulation_xml, imd_force=imd_force)
    assert simulation.currentStep == 0
    simulation.step(5)
    assert simulation.currentStep == 5


@pytest.mark.parametrize("section_missing", ("System", "Integrator"))
def test_missing_section(basic_simulation_xml, section_missing):
    """
    Make sure that the deserialization raise an exception when a base section
    is missing.
    """
    broken_xml = remove_xml_tag(basic_simulation_xml, section_missing)
    with pytest.raises(IOError):
        deserialize_simulation(broken_xml)


@pytest.mark.parametrize(
    "section_to_duplicate", ("pdbx", "pdb", "System", "Integrator")
)
def test_duplicate_section(basic_simulation_xml, section_to_duplicate):
    """
    Make sure de deserialization fails if a base section appears more than once.
    """
    broken_xml = add_extra_xml_tag(basic_simulation_xml, section_to_duplicate)
    with pytest.raises(IOError):
        deserialize_simulation(broken_xml)


def test_missing_structure(basic_simulation_xml):
    broken_xml = remove_xml_tag(basic_simulation_xml, "pdbx")
    broken_xml = remove_xml_tag(broken_xml, "pdb")
    with pytest.raises(IOError):
        deserialize_simulation(broken_xml)


def test_duplicate_structure(basic_simulation_xml):
    broken_xml = add_extra_xml_tag(basic_simulation_xml, "pdb")
    with pytest.raises(IOError):
        deserialize_simulation(broken_xml)


def test_imd_force(basic_simulation_xml, empty_imd_force):
    """
    When deserializing a simulation and passing an imd force, the force is
    added to the system.
    """
    simulation = deserialize_simulation(basic_simulation_xml, empty_imd_force)
    # The provided force is appended as the last step, so will be the last listed force
    force_obtained = simulation.system.getForces()[-1]
    force_added = empty_imd_force
    # The forces are the same if by modifying one we also modify the other.
    force_added.setParticleParameters(0, 0, (1.0, 2.0, 3.0))
    parameters = force_obtained.getParticleParameters(0)
    assert parameters == [0, (1.0, 2.0, 3.0)]


@pytest.mark.parametrize("platform", ("Reference", "CPU"))
def test_platform(basic_simulation_xml, platform):
    """
    We can choose the platform when deserialising a simulation.

    Only the "Reference" and "CPU" platforms can be expected to be available so
    we only test these ones.
    """
    simulation = deserialize_simulation(basic_simulation_xml, platform_name=platform)
    effective_platform_name = simulation.context.getPlatform().getName()
    assert effective_platform_name == platform


def test_serializer_pbc():
    """
    Check if the periodic boundary conditions are correctly serialized.
    The test deserializes two simulations with different setting of pbc and check if the positions are correct.
    """
    omm_sim = build_basic_simulation()
    UNIT_SIMULATION_BOX_VECTORS = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    omm_sim.context.setPeriodicBoxVectors(*UNIT_SIMULATION_BOX_VECTORS)

    def out_of_bounds(coord):
        return coord < 0 or coord > 1

    def get_sim_position_coords(sim):
        for position in sim.make_regular_frame().particle_positions:
            for coord in position:
                yield coord

    with StringIO() as xml_pbc:
        xml_pbc.write(serialize_simulation(omm_sim, pbc_wrapping=True))
        xml_pbc.seek(0)
        sim_pbc_obj = deserialize_simulation(xml_pbc.read())

    openmm_sim = OpenMMSimulation.from_simulation(sim_pbc_obj)
    openmm_sim.load()
    with NanoverImdApplication.basic_server(port=0) as app_server:
        openmm_sim.reset(app_server)

    assert not any(
        out_of_bounds(coord) for coord in get_sim_position_coords(openmm_sim)
    )

    with StringIO() as xml_no_pbc:
        xml_no_pbc.write(serialize_simulation(omm_sim, pbc_wrapping=False))
        xml_no_pbc.seek(0)
        sim_no_pbc_obj = deserialize_simulation(xml_no_pbc.read())

    openmm_sim = OpenMMSimulation.from_simulation(sim_no_pbc_obj)
    openmm_sim.load()
    with NanoverImdApplication.basic_server(port=0) as app_server:
        openmm_sim.reset(app_server)

    assert any(out_of_bounds(coord) for coord in get_sim_position_coords(openmm_sim))
