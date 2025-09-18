"""
Fixtures and utilities for tests that requires OpenMM simulations.
"""

# Pylint does not recognize pytest fixtures, which causes some false warnings.
# pylint: disable=unused-argument,redefined-outer-name
import pytest
import numpy as np

import openmm as mm
from openmm import app

# Prefixed units in `openmm.unit` are added programmatically and are not
# recognized by pylint and PyCharm.
from openmm.unit import (
    kelvin,
    picosecond,
    femtosecond,
    nanometer,
)  # pylint: disable=no-name-in-module

from nanover.openmm import serializer
import nanover.openmm.imd


BASIC_SIMULATION_BOX_VECTORS = [[50, 0, 0], [0, 50, 0], [0, 0, 50]]
BASIC_SIMULATION_POSITIONS = [
    # First residue
    [0, 0, 0],  # C
    [5.288, 1.610, 9.359],  # H
    [2.051, 8.240, -6.786],  # H
    [-10.685, -0.537, 1.921],  # H
    # Second residue, copied from the first but shifted
    # by 5 nm along the Z axis
    [0, 0, 5],  # C
    [5.288, 1.610, 14.359],  # H
    [2.051, 8.240, -1.786],  # H
    [-10.685, -0.537, 6.921],  # H
]
ARGON_SIMULATION_POSITION = [[0.0, 0.0, 0.0]]


def build_basic_system():
    periodic_box_vector = BASIC_SIMULATION_BOX_VECTORS
    system = mm.System()
    system.setDefaultPeriodicBoxVectors(*periodic_box_vector)
    system.addParticle(mass=12)
    system.addParticle(mass=1)
    system.addParticle(mass=1)
    system.addParticle(mass=1)
    system.addParticle(mass=12)
    system.addParticle(mass=1)
    system.addParticle(mass=1)
    system.addParticle(mass=1)
    return system


def build_basic_topology() -> app.Topology:
    topology = app.Topology()
    carbon = app.Element.getBySymbol("C")
    hydrogen = app.Element.getBySymbol("H")
    chain = topology.addChain()
    residue = topology.addResidue(name="METH1", chain=chain)
    atom_c1 = topology.addAtom(element=carbon, name="C1", residue=residue)
    atom_h2 = topology.addAtom(element=hydrogen, name="H2", residue=residue)
    atom_h3 = topology.addAtom(element=hydrogen, name="H3", residue=residue)
    atom_h4 = topology.addAtom(element=hydrogen, name="H4", residue=residue)
    topology.addBond(atom_c1, atom_h2)
    topology.addBond(atom_c1, atom_h3)
    topology.addBond(atom_c1, atom_h4)
    chain = topology.addChain()
    residue = topology.addResidue(name="METH2", chain=chain)
    atom_c1 = topology.addAtom(element=carbon, name="C1", residue=residue)
    atom_h2 = topology.addAtom(element=hydrogen, name="H2", residue=residue)
    atom_h3 = topology.addAtom(element=hydrogen, name="H3", residue=residue)
    atom_h4 = topology.addAtom(element=hydrogen, name="H4", residue=residue)
    topology.addBond(atom_c1, atom_h2)
    topology.addBond(atom_c1, atom_h3)
    topology.addBond(atom_c1, atom_h4)
    return topology


def build_basic_simulation():
    """
    Setup a minimal OpenMM simulation with two methane molecules.
    """
    # In ths function, we define matrices and we want to align the column.
    # We disable the pylint warning about bad spacing for the scope of the
    # function.
    # pylint: disable=bad-whitespace
    periodic_box_vector = BASIC_SIMULATION_BOX_VECTORS
    positions = np.array(BASIC_SIMULATION_POSITIONS, dtype=np.float32)

    topology = build_basic_topology()
    system = build_basic_system()

    force = mm.NonbondedForce()
    force.setNonbondedMethod(force.NoCutoff)
    # These non-bonded parameters are completely wrong, but it does not matter
    # for the tests as long as we do not start testing the dynamic and
    # thermodynamics properties of methane.
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    system.addForce(force)

    integrator = mm.LangevinIntegrator(300 * kelvin, 1 / picosecond, 2 * femtosecond)

    platform = mm.Platform.getPlatformByName("CPU")
    simulation = app.Simulation(topology, system, integrator, platform=platform)
    simulation.context.setPeriodicBoxVectors(*periodic_box_vector)
    simulation.context.setPositions(positions * nanometer)

    return simulation


@pytest.fixture
def basic_system():
    return build_basic_system()


@pytest.fixture
def basic_simulation():
    return build_basic_simulation()


@pytest.fixture
def basic_simulation_xml():
    """
    Generate an XML serialized simulation from the basic test simulation.
    """
    xml_string = serializer.serialize_simulation(
        build_basic_simulation(), save_state=True
    )
    return xml_string


@pytest.fixture
def empty_imd_force():
    return nanover.openmm.imd.create_imd_force()
