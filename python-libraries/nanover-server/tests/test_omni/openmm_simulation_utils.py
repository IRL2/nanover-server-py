"""
Fixtures and utilities for tests that requires OpenMM simulations.
"""

# Pylint does not recognize pytest fixtures, which causes some false warnings.
# pylint: disable=unused-argument,redefined-outer-name
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


def build_basic_simulation_periodic():
    """
    Setup a minimal OpenMM simulation with two methane molecules.
    """
    periodic_box_vector = BASIC_SIMULATION_BOX_VECTORS
    positions = np.array(BASIC_SIMULATION_POSITIONS, dtype=np.float32)

    topology = build_basic_topology()
    system = build_basic_system()

    force = mm.NonbondedForce()
    force.setNonbondedMethod(force.PME)
    force.setCutoffDistance(1.2 * nanometer)
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


def assert_basic_simulation_topology(frame):
    """
    Fails with an :exc:`AssertError` if the topology of the given frame does
    not match the expectation for the basic simulation.
    """
    assert frame.residue_names == ["METH1", "METH2"]
    assert frame.residue_chains == [0, 1]
    assert frame.particle_names == ["C1", "H2", "H3", "H4"] * 2
    assert frame.particle_elements == [6, 1, 1, 1] * 2
    assert frame.particle_residues == [0] * 4 + [1] * 4
    assert frame.bond_pairs == [
        [0, 1],
        [0, 2],
        [0, 3],  # First residue
        [4, 5],
        [4, 6],
        [4, 7],  # Second residue
    ]


def build_single_atom_system():
    box_vector = BASIC_SIMULATION_BOX_VECTORS
    system = mm.System()
    system.setDefaultPeriodicBoxVectors(*box_vector)
    system.addParticle(mass=40)

    return system


def build_single_atom_topology() -> app.Topology:
    topology = app.Topology()
    argon = app.Element.getBySymbol("Ar")
    chain = topology.addChain()
    residue = topology.addResidue(name="ARGON", chain=chain)
    topology.addAtom(element=argon, name="AR1", residue=residue)
    return topology


def build_single_atom_simulation():
    periodic_box_vector = BASIC_SIMULATION_BOX_VECTORS
    positions = np.array(ARGON_SIMULATION_POSITION, dtype=np.float32)
    topology = build_single_atom_topology()
    system = build_single_atom_system()

    # As we are only dealing with a single atom system, it is unnecessary to
    # add a non-bonded force.

    # Use a Verlet integrator to make the dynamics predictable, avoiding the
    # random kicks introduced by Langevin.
    integrator = mm.VerletIntegrator(2 * femtosecond)

    platform = mm.Platform.getPlatformByName("CPU")
    simulation = app.Simulation(topology, system, integrator, platform=platform)
    simulation.context.setPeriodicBoxVectors(*periodic_box_vector)
    simulation.context.setPositions(positions * nanometer)

    return simulation


def assert_single_atom_simulation_topology(frame):
    """
    Fails with an :exc:`AssertError` if the topology of the given frame does
    not match the expectation for the basic simulation.
    """
    assert frame.residue_names == ["ARGON"]
    assert frame.residue_chains == [0]
    assert frame.particle_names == ["AR1"]
    assert frame.particle_elements == [40]
