# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Fixtures and utilities for tests that requires OpenMM simulations.
"""
# Pylint does not recognize pytest fixtures, which causes some false warnings.
# pylint: disable=unused-argument,redefined-outer-name
import pytest
from typing import Optional
import numpy as np

import simtk.openmm as mm
from simtk.openmm import app
# Prefixed units in `simtk.unit` are added programmatically and are not
# recognized by pylint and PyCharm.
from simtk.unit import kelvin, picosecond, femtosecond, nanometer  # pylint: disable=no-name-in-module

from narupa.openmm import serializer
import narupa.openmm.imd


BASIC_SIMULATION_BOX_VECTORS = [
    [50,  0,  0],
    [ 0, 50,  0],
    [ 0,  0, 50]
]
BASIC_SIMULATION_POSITIONS = [
    # First residue
    [ 0,       0,      0],  # C
    [ 5.288,   1.610,  9.359],  # H
    [ 2.051,   8.240, -6.786],  # H
    [-10.685, -0.537,  1.921],  # H
    # Second residue, copied from the first but shifted
    # by 5 nm along the Z axis
    [  0,      0,      5],  # C
    [  5.288,  1.610, 14.359],  # H
    [  2.051,  8.240, -1.786],  # H
    [-10.685, -0.537,  6.921],  # H
]


class DoNothingReporter:
    """
    OpenMM reporter that does nothing.

    The reporter does nothing but is valid. It is meant to populate the list of
    reporters of an OpenMM simulation.
    """
    # The name of the method is part of the OpenMM API. It cannot be made to
    # conform PEP8.
    def describeNextReport(self, simulation):  # pylint: disable=invalid-name,no-self-use
        """
        Activate the reporting every step, but collect no data.
        """
        return 1, False, False, False, False

    def report(self, simulation, state):
        """
        Do not report anything.
        """
        pass


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
    carbon = app.Element.getBySymbol('C')
    hydrogen = app.Element.getBySymbol('H')
    chain = topology.addChain()
    residue = topology.addResidue(name='METH1', chain=chain)
    atom_c1 = topology.addAtom(element=carbon, name='C1', residue=residue)
    atom_h2 = topology.addAtom(element=hydrogen, name='H2', residue=residue)
    atom_h3 = topology.addAtom(element=hydrogen, name='H3', residue=residue)
    atom_h4 = topology.addAtom(element=hydrogen, name='H4', residue=residue)
    topology.addBond(atom_c1, atom_h2)
    topology.addBond(atom_c1, atom_h3)
    topology.addBond(atom_c1, atom_h4)
    chain = topology.addChain()
    residue = topology.addResidue(name='METH2', chain=chain)
    atom_c1 = topology.addAtom(element=carbon, name='C1', residue=residue)
    atom_h2 = topology.addAtom(element=hydrogen, name='H2', residue=residue)
    atom_h3 = topology.addAtom(element=hydrogen, name='H3', residue=residue)
    atom_h4 = topology.addAtom(element=hydrogen, name='H4', residue=residue)
    topology.addBond(atom_c1, atom_h2)
    topology.addBond(atom_c1, atom_h3)
    topology.addBond(atom_c1, atom_h4)
    return topology


def build_basic_simulation(
        imd_force: Optional[mm.CustomExternalForce] = None) -> app.Simulation:
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
    if imd_force is not None:
        narupa.openmm.imd.populate_imd_force(imd_force, system)
        system.addForce(imd_force)

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

    platform = mm.Platform.getPlatformByName('CPU')
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
def basic_simulation_with_imd_force():
    imd_force = narupa.openmm.imd.create_imd_force()
    return build_basic_simulation(imd_force), imd_force


@pytest.fixture
def serialized_simulation_path(basic_simulation, tmp_path):
    """
    Setup an XML serialized simulation as a temporary files.
    """
    serialized_simulation = serializer.serialize_simulation(basic_simulation)
    xml_path = tmp_path / "system.xml"
    with open(str(xml_path), 'w') as outfile:
        outfile.write(serialized_simulation)
    return xml_path


@pytest.fixture
def basic_simulation_xml(basic_simulation):
    """
    Generate a XML serialized simulation from the basic test simulation.
    """
    xml_string = serializer.serialize_simulation(basic_simulation)
    return xml_string


@pytest.fixture
def empty_imd_force():
    return narupa.openmm.imd.create_imd_force()


def assert_basic_simulation_topology(frame):
    """
    Fails with an :exc:`AssertError` if the topology of the given frame does
    not match the expectation for the basic simulation.
    """
    assert frame.residue_names == ['METH1', 'METH2']
    assert frame.residue_chains == [0, 1]
    assert frame.particle_names == ['C1', 'H2', 'H3', 'H4'] * 2
    assert frame.particle_elements == [6, 1, 1, 1] * 2
    assert frame.particle_residues == [0] * 4 + [1] * 4
    assert frame.bond_pairs == [
        [0, 1], [0, 2], [0, 3],  # First residue
        [4, 5], [4, 6], [4, 7],  # Second residue
    ]