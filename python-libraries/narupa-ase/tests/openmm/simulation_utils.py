# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Fixtures and utilities for tests that requires OpenMM simulations.
"""
# TODO: This is a duplicated file from narupa-openmm. See issue #60.

# Pylint does not recognize pytest fixtures, which causes some false warnings.
# pylint: disable=unused-argument,redefined-outer-name
import pytest

import numpy as np

import simtk.openmm as mm
from simtk.openmm import app
# Prefixed units in `simtk.unit` are added programmatically and are not
# recognized by pylint and PyCharm.
from simtk.unit import kelvin, picosecond, femtosecond, nanometer  # pylint: disable=no-name-in-module

from narupa.openmm import serializer


def build_basic_simulation():
    # In ths function, we define matrices and we want to align the column.
    # We disable the pylint warning about bad spacing for the scope of the
    # function.
    # pylint: disable=bad-whitespace
    periodic_box_vector = [
        [15,  0,  0],
        [ 0, 15,  0],
        [ 0,  0, 15]
    ]
    positions = np.array([
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
    ], dtype=np.float32)

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

    force = mm.NonbondedForce()
    force.setNonbondedMethod(force.CutoffPeriodic)
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
def basic_simulation():
    """
    Setup a minimal OpenMM simulation with two methane molecules.
    """
    return build_basic_simulation()


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
