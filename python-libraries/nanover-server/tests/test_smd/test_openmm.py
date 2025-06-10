"""
Tests for :mod:`nanover.smd.openmm`.

Things to test:
- An SMD simulation can be created either from an existing OpenMM simulation or
  a NanoVer OpenMM XML file
- OpenMMSMDSimulation returns OpenMMSMDSimulationAtom or OpenMMSMDSimulationCOM
  as appropriate
- The PBCs of the loaded simulation are respected by the SMD force, and the SMD
  force shares this periodicity
- The SMD force attaches to the correct atom (dictated by the index/indices passed
  to the class upon creation)
- The simulation can be reset correctly, with all attributes returning to the same
  state as immediately after the creation of the class itself
- SMD force is correctly added to the system
- SMD force can be correctly removed from the system
- SMD force position is correctly updated
- Running the equilibration with the initial restraint throws an error correctly
  when run with the SMD force not located at the initial position
- The SMD simulation can be saved correctly, with or without the SMD force
- If the SMD simulation is being loaded from a NanoVer OpenMM XML file created
  via one of the SMDSimulation classes and contains an SMD force already, this
  SMD force is correctly loaded and matches the expected force constant specified
  when creating the class
- The class can generate the correct number of starting structures in the specified
  time interval, and that these are saved to the correct location
- Running an SMD simulation produces reasonable results for the cumulative work done
  (may need to think about a specific test case for this...)
- _calculate_forces works as expected
- _calculate_work_done works as expected
- Simulation data is saved in the correct format to the correct location, and can be
  subsequently loaded back into python correctly
- General SMD data is saved in the correct format to the correct location, and can be
  subsequently loaded back into python correctly
- For OpenMMSMDSimulationCOM, the COM of the specified atoms is correctly calculated
- For OpenMMSMDSimulationCOM, the trajectory of the COM of the specified atoms is
  correctly calculated
- smd_com_force works as expected
- smd_single_atom_force works as expected
"""
import pytest
from contextlib import contextmanager

import numpy as np
import openmm as mm
from openmm import app
from openmm.unit import (
    kelvin,
    picosecond,
    femtosecond,
    nanometer,
)

from nanover.openmm import serializer
from nanover.smd.openmm import *

# Very basic thing to test entire class as it would be used: tutorial notebook that can be tested

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

TEST_SMD_SINGLE_INDEX = np.array(0)
TEST_SMD_MULTIPLE_INDICES = np.array([0, 1, 2, 3])
TEST_SMD_PATH = np.array([np.linspace(0.0,1.0,101), np.zeros(101), np.zeros(101)]).transpose()
TEST_SMD_FORCE_CONSTANT = 3011.0

TEST_SMD_INDEX_OPTIONS = [TEST_SMD_SINGLE_INDEX, TEST_SMD_MULTIPLE_INDICES]
EXPECTED_SMD_SIMULATION_TYPES = [OpenMMSMDSimulationAtom, OpenMMSMDSimulationCOM]

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


def build_basic_simulation() -> app.Simulation:
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
def make_basic_simulation_xml(tmp_path):
    serialized_simulation = serializer.serialize_simulation(build_basic_simulation(), save_state=True)
    xml_path = tmp_path / "basic_simulation.xml"
    with open(str(xml_path), "w") as xml_file:
        xml_file.write(serialized_simulation)
    return xml_path

@pytest.fixture
def make_smd_simulation(simulation, indices):
    """
    Fixture that creates an OpenMMSMDSimulation from
    """
    smd_sim = OpenMMSMDSimulation.from_simulation(simulation, indices, TEST_SMD_PATH, TEST_SMD_FORCE_CONSTANT)
    return smd_sim


def test_load_smd_sim_from_simulation():
    """
    Test that an OpenMMSMDSimulation can be correctly loaded from an OpenMM simulation.
    """
    smd_sim = OpenMMSMDSimulation.from_simulation(build_basic_simulation(), TEST_SMD_SINGLE_INDEX, TEST_SMD_PATH, TEST_SMD_FORCE_CONSTANT)
    assert smd_sim
    assert smd_sim.simulation
    assert (smd_sim.smd_path == TEST_SMD_PATH).all()
    assert (smd_sim.smd_atom_indices == TEST_SMD_SINGLE_INDEX).all()
    assert smd_sim.smd_force_constant == TEST_SMD_FORCE_CONSTANT


def test_load_smd_sim_from_xml_path(make_basic_simulation_xml):
    """
    Test that an OpenMMSMDSimulation can be correctly loaded from a NanoVer OpenMM XML file.
    """
    smd_sim = OpenMMSMDSimulation.from_xml_path(make_basic_simulation_xml, TEST_SMD_SINGLE_INDEX, TEST_SMD_PATH, TEST_SMD_FORCE_CONSTANT)
    assert smd_sim
    assert smd_sim.xml_path == make_basic_simulation_xml
    assert smd_sim.simulation
    assert (smd_sim.smd_path == TEST_SMD_PATH).all()
    assert (smd_sim.smd_atom_indices == TEST_SMD_SINGLE_INDEX).all()
    assert smd_sim.smd_force_constant == TEST_SMD_FORCE_CONSTANT
