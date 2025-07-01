"""
Tests for :mod:`nanover.smd.openmm`.

Things to test:
- An SMD simulation can be created either from an existing OpenMM simulation or
  a NanoVer OpenMM XML file [√]
- OpenMMSMDSimulation returns OpenMMSMDSimulationAtom or OpenMMSMDSimulationCOM
  as appropriate [√]
- The PBCs of the loaded simulation are respected by the SMD force, and the SMD
  force shares this periodicity [√]
- The SMD force attaches to the correct atom (dictated by the index/indices passed
  to the class upon creation) [√]
- The simulation can be reset correctly, with all attributes returning to the same
  state as immediately after the creation of the class itself [√]
- SMD force is correctly added to the system [√]
- SMD force can be correctly removed from the system [√]
- SMD force position is correctly updated [√]
- Running the equilibration with the initial restraint throws an error correctly
  when run with the SMD force not located at the initial position [√]
- The SMD simulation can be saved correctly, with or without the SMD force [√]
- If the SMD simulation is being loaded from a NanoVer OpenMM XML file created
  via one of the SMDSimulation classes and contains an SMD force already, this
  SMD force is correctly loaded and matches the expected force constant specified
  when creating the class [√]
- The class can generate the correct number of starting structures in the specified
  time interval, and that these are saved to the correct location [√]
- Running an SMD simulation produces reasonable results for the cumulative work done [√]
- _calculate_smd_forces works as expected [√]
- _calculate_work_done works as expected [√]
- Simulation data is saved in the correct format to the correct location, and can be
  subsequently loaded back into python correctly [√]
- General SMD data is saved in the correct format to the correct location, and can be
  subsequently loaded back into python correctly [√]
- The COM of a specified group of atoms is correctly calculated [√]
- For OpenMMSMDSimulationCOM, the COM of the specified atoms is correctly calculated [√]
- For OpenMMSMDSimulationCOM, the trajectory of the COM of the specified atoms is
  correctly calculated [√]
- smd_com_force works as expected [√]
- smd_single_atom_force works as expected [√]
- OpenMMSMDSimulation correctly loads the state of a simulation [√]
"""

import tempfile
from io import StringIO

import numpy as np
import pytest
from contextlib import redirect_stdout

import openmm as mm
from openmm import app
from openmm.unit import (
    kelvin,
    picosecond,
    femtosecond,
    nanometer,
)

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

# Test parameters for OpenMMSMDSimulation
TEST_SMD_SINGLE_INDEX = np.array(0)
TEST_SMD_MULTIPLE_INDICES = np.array([0, 1, 2, 3])
TEST_SMD_PATH = np.array(
    [np.linspace(0.05, 1.05, 101), np.zeros(101), np.zeros(101)]
).transpose()
TEST_SMD_FORCE_CONSTANT = 3011.0

TEST_SMD_ARGON_INDEX = np.array(0)
TEST_SMD_ARGON_PATH = np.array(
    [np.linspace(0.00, 0.02, 3), np.zeros(3), np.zeros(3)]
).transpose()
TEST_SMD_ARGON_FORCE_CONSTANT = 100.0

# Test systems for COM calculations, formatted as (positions, masses, expected COM)
TEST_COM_TWO_ATOMS = (
    np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
    np.array([1.0, 1.0]),
    np.array([0.0, 0.0, 0.0]),
)
TEST_COM_METHANE = (
    np.array(
        [
            [1.0, 1.0, 0.0],
            [0.0, 1.0, -1.0 / np.sqrt(2)],
            [2.0, 1.0, -1.0 / np.sqrt(2)],
            [1.0, 0.0, 1.0 / np.sqrt(2)],
            [1.0, 2.0, 1.0 / np.sqrt(2)],
        ]
    ),
    np.array([12.01, 1.00, 1.00, 1.00, 1.00]),
    np.array([1.0, 1.0, 0.0]),
)
TEST_COM_CIRCLE = (
    np.array(
        [[np.cos(i * (np.pi / 4)), 0.0, np.sin(i * (np.pi / 4))] for i in range(8)]
    ),
    np.array([2.0, 1.05, 2.0, 1.05, 2.0, 1.05, 2.0, 1.05]),
    np.array([0.0, 0.0, 0.0]),
)
TEST_COM_CUBE = (
    np.array(
        [
            [0.0, 1.0, 2.0],
            [0.0, 1.0, 4.0],
            [0.0, 3.0, 2.0],
            [0.0, 3.0, 4.0],
            [2.0, 1.0, 2.0],
            [2.0, 1.0, 4.0],
            [2.0, 3.0, 2.0],
            [2.0, 3.0, 4.0],
        ]
    ),
    np.array([6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0]),
    np.array([1.0, 2.0, 3.0]),
)


def build_com_system(parameters: tuple):
    positions, masses, com = parameters
    box_vector = BASIC_SIMULATION_BOX_VECTORS
    system = mm.System()
    system.setDefaultPeriodicBoxVectors(*box_vector)
    for atom in range(masses.size):
        system.addParticle(mass=masses[atom])

    return system


def build_com_topology(parameters: tuple) -> app.Topology:
    positions, masses, com = parameters
    topology = app.Topology()
    for atom in range(masses.size):
        element = app.Element.getByMass(masses.size)
        chain = topology.addChain()
        residue = topology.addResidue(name=f"atom_{atom}", chain=chain)
        topology.addAtom(element=element, name=f"AT{atom}", residue=residue)
    return topology


def build_com_simulation(parameters: tuple) -> app.Simulation:
    positions, masses, com = parameters
    periodic_box_vector = BASIC_SIMULATION_BOX_VECTORS
    topology = build_com_topology(parameters)
    system = build_com_system(parameters)

    # No forces added to system, non-interacting particles
    integrator = mm.LangevinIntegrator(300 * kelvin, 1 / picosecond, 2 * femtosecond)

    platform = mm.Platform.getPlatformByName("CPU")
    simulation = app.Simulation(topology, system, integrator, platform=platform)
    simulation.context.setPeriodicBoxVectors(*periodic_box_vector)
    simulation.context.setPositions(positions * nanometer)

    return simulation


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


def build_basic_simulation(pbcs: bool = False) -> app.Simulation:
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
    if pbcs:
        force.setNonbondedMethod(force.CutoffPeriodic)
    else:
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


def build_single_atom_simulation(pbcs: bool = False):
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


@pytest.fixture
def make_basic_simulation_xml(tmp_path):
    serialized_simulation = serializer.serialize_simulation(
        build_basic_simulation(), save_state=True
    )
    xml_path = tmp_path / "basic_simulation.xml"
    with open(str(xml_path), "w") as xml_file:
        xml_file.write(serialized_simulation)
    return xml_path


@pytest.fixture
def make_basic_smd_simulation_with_atom_smd_force_xml(tmp_path):
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        TEST_SMD_SINGLE_INDEX,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    xml_path = tmp_path / "basic_smd_simulation.xml"
    smd_sim.save_simulation(xml_path, save_state=True, save_smd_force=True)
    return xml_path


@pytest.fixture
def make_basic_smd_simulation_with_com_smd_force_xml(tmp_path):
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        TEST_SMD_MULTIPLE_INDICES,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    xml_path = tmp_path / "basic_smd_simulation.xml"
    smd_sim.save_simulation(xml_path, save_state=True, save_smd_force=True)
    return xml_path


@pytest.fixture
def make_basic_smd_simulation_without_atom_smd_force_xml(tmp_path):
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        TEST_SMD_SINGLE_INDEX,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    xml_path = tmp_path / "basic_smd_simulation.xml"
    smd_sim.save_simulation(xml_path, save_state=True, save_smd_force=False)
    return xml_path


@pytest.fixture
def make_basic_smd_simulation_without_com_smd_force_xml(tmp_path):
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        TEST_SMD_MULTIPLE_INDICES,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    xml_path = tmp_path / "basic_smd_simulation.xml"
    smd_sim.save_simulation(xml_path, save_state=True, save_smd_force=False)
    return xml_path


def test_load_smd_sim_from_simulation():
    """
    Test that an OpenMMSMDSimulation can be correctly loaded from an OpenMM simulation.
    """
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        TEST_SMD_SINGLE_INDEX,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    assert smd_sim
    assert smd_sim.simulation
    assert np.array_equal(smd_sim.smd_path, TEST_SMD_PATH)
    assert np.array_equal(smd_sim.smd_atom_indices, TEST_SMD_SINGLE_INDEX)
    assert smd_sim.smd_force_constant == TEST_SMD_FORCE_CONSTANT


def test_load_smd_sim_from_xml_path(make_basic_simulation_xml):
    """
    Test that an OpenMMSMDSimulation can be correctly loaded from a NanoVer OpenMM XML file.
    """
    smd_sim = OpenMMSMDSimulation.from_xml_path(
        make_basic_simulation_xml,
        TEST_SMD_SINGLE_INDEX,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    assert smd_sim
    assert smd_sim.xml_path == make_basic_simulation_xml
    assert smd_sim.simulation
    assert np.array_equal(smd_sim.smd_path, TEST_SMD_PATH)
    assert np.array_equal(smd_sim.smd_atom_indices, TEST_SMD_SINGLE_INDEX)
    assert smd_sim.smd_force_constant == TEST_SMD_FORCE_CONSTANT


def test_load_smd_simulation_with_atom_smd_force_from_xml_path(
    make_basic_smd_simulation_with_atom_smd_force_xml,
):
    """
    Check that when an input file containing a single atom SMD force is passed to the
    OpenMMSMDSimulation class, the SMD force is loaded correctly from the file using
    check_for_existing_smd_force(), and that the parameters for the SMD force match
    those that are passed via the file.
    """
    with redirect_stdout(StringIO()) as _:
        smd_sim = OpenMMSMDSimulation.from_xml_path(
            make_basic_smd_simulation_with_atom_smd_force_xml,
            TEST_SMD_SINGLE_INDEX,
            TEST_SMD_PATH,
            TEST_SMD_FORCE_CONSTANT,
        )
        assert smd_sim.loaded_smd_force_from_sim


def test_load_smd_simulation_with_com_smd_force_from_xml_path(
    make_basic_smd_simulation_with_com_smd_force_xml,
):
    """
    Check that when an input file containing a COM SMD force is passed to the
    OpenMMSMDSimulation class, the SMD force is loaded correctly from the file using
    check_for_existing_smd_force(), and that the parameters for the SMD force match
    those that are passed via the file.
    """
    with redirect_stdout(StringIO()) as _:
        smd_sim = OpenMMSMDSimulation.from_xml_path(
            make_basic_smd_simulation_with_com_smd_force_xml,
            TEST_SMD_MULTIPLE_INDICES,
            TEST_SMD_PATH,
            TEST_SMD_FORCE_CONSTANT,
        )
        assert smd_sim.loaded_smd_force_from_sim


def test_load_smd_simulation_without_atom_smd_force_from_xml_path(
    make_basic_smd_simulation_without_atom_smd_force_xml,
):
    """
    Check that when an xml input file is saved from an OpenMMSMDSimulationAtom class
    without the SMD force, the SMD force is not loaded from the file.
    """
    smd_sim = OpenMMSMDSimulation.from_xml_path(
        make_basic_smd_simulation_without_atom_smd_force_xml,
        TEST_SMD_SINGLE_INDEX,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    assert not smd_sim.loaded_smd_force_from_sim


def test_load_smd_simulation_without_com_smd_force_from_xml_path(
    make_basic_smd_simulation_without_com_smd_force_xml,
):
    """
    Check that when an xml input file is saved from an OpenMMSMDSimulationCOM class
    without the SMD force, the SMD force is not loaded from the file.
    """
    smd_sim = OpenMMSMDSimulation.from_xml_path(
        make_basic_smd_simulation_without_com_smd_force_xml,
        TEST_SMD_MULTIPLE_INDICES,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    assert not smd_sim.loaded_smd_force_from_sim


@pytest.mark.parametrize(
    "indices, sim_type",
    [
        (TEST_SMD_SINGLE_INDEX, OpenMMSMDSimulationAtom),
        (TEST_SMD_MULTIPLE_INDICES, OpenMMSMDSimulationCOM),
    ],
)
def test_return_correct_smd_sim_type(indices, sim_type):
    """
    Check that the OpenMMSMDSimulation class returns the correct subclass depending on the
    number of indices that are passed to it (one for OpenMMSMDSimulationAtom, more than one
    for OpenMMSMDSimulationCOM).

    :param indices: Indices of atoms to apply the SMD force to (should at least
      test one single index and one set of indices)
    :param sim_type: Type of simulation to expect for the indices given
    """
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    assert type(smd_sim) == sim_type


@pytest.mark.parametrize("apply_pbcs", [True, False])
@pytest.mark.parametrize("indices", [TEST_SMD_SINGLE_INDEX, TEST_SMD_MULTIPLE_INDICES])
def test_simulation_pbcs_are_respected(apply_pbcs, indices):
    """
    Check that the periodic boundary conditions of the OpenMMSimulation passed to the
    OpenMMSMDSimulation class are respected (i.e. the PBCs of the SMD simulation match
    those of the OpenMM simulation), and that the PBCs of the SMD force match the PBCs
    of the simulation.

    :param apply_pbcs: Boolean value indicating whether to apply PBCs to the simulation
    :param indices: Indices of atoms to apply the SMD force to (should at least
      test one single index and one set of indices)
    """
    sim = build_basic_simulation(pbcs=apply_pbcs)
    uses_pbcs = sim.system.usesPeriodicBoundaryConditions()
    smd_sim = OpenMMSMDSimulation.from_simulation(
        sim,
        indices,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    assert smd_sim.smd_force.usesPeriodicBoundaryConditions() == uses_pbcs
    assert smd_sim.simulation.system.usesPeriodicBoundaryConditions() == uses_pbcs


@pytest.mark.parametrize(
    "index", [np.array(0), np.array(1), np.array(4), np.array(5), np.array(7)]
)
def test_smd_force_attaches_to_correct_atom(index):
    """
    Check that the SMD force is attached to the correct atom when a single index is passed.
    Should use the OpenMMSMDSimulationAtom class, with only one CustomExternalForce.

    :param index: Indices of atoms to apply the SMD force to (should be arrays containing
      a single index)
    """
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        index,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    # Attaches force to single atom, so index of atom within force is zero
    p_index, p_params = smd_sim.smd_force.getParticleParameters(0)
    assert p_index == index
    assert np.array_equal(np.array(p_params), TEST_SMD_PATH[0])


@pytest.mark.parametrize(
    "indices",
    [
        np.array([0, 1, 2, 3]),
        np.array([1, 2, 3, 4]),
        np.array([0, 1, 4, 5]),
        np.array([1, 3, 4, 7]),
    ],
)
def test_smd_force_attaches_to_correct_atoms(indices):
    """
    Check that the SMD force attaches to the correct atoms when an array of indices is passed.
    Should use the OpenMMSMDSimulationCOM class, with only one CustomCentroidBondForce.

    :param indices: Indices of atoms to apply the SMD force to (should be arrays of multiple
      indices)
    """
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    # Only one centroid force added, index of force is zero
    p_indices, _ = smd_sim.smd_force.getGroupParameters(0)
    assert np.array_equal(np.array(p_indices), indices)


@pytest.mark.parametrize("indices", [TEST_SMD_SINGLE_INDEX, TEST_SMD_MULTIPLE_INDICES])
def test_reset(indices):
    """
    Check that all the attributes of the OpenMMSMDSimulation subclasses are reset to their initial
    state by the .reset() function of the class.

    :param indices: Indices of atoms to apply the SMD force to (should at least
      test one single index and one set of indices)
    """
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )

    with redirect_stdout(StringIO()) as _:
        smd_sim.run_smd()
        smd_sim.reset()

    # Create a fresh copy of the SMD simulation
    smd_sim_copy = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )

    # Check that the attributes that were created during the SMD simulation
    # are no longer present in the class
    try:
        assert smd_sim.smd_simulation_forces or smd_sim_copy.smd_simulation_work_done
    except AttributeError:
        pass

    # Check relevant observables from the simulation against the fresh simulation
    smd_sim_state = smd_sim.simulation.context.getState(
        getPositions=True, getVelocities=True, getForces=True, getEnergy=True
    )
    smd_sim_copy_state = smd_sim_copy.simulation.context.getState(
        getPositions=True, getVelocities=True, getForces=True, getEnergy=True
    )
    assert np.array_equal(
        smd_sim_state.getPositions(asNumpy=True),
        smd_sim_copy_state.getPositions(asNumpy=True),
    )
    assert np.array_equal(
        smd_sim_state.getVelocities(asNumpy=True),
        smd_sim_copy_state.getVelocities(asNumpy=True),
    )
    assert np.array_equal(
        smd_sim_state.getForces(asNumpy=True),
        smd_sim_copy_state.getForces(asNumpy=True),
    )
    assert smd_sim_state.getKineticEnergy() == smd_sim_copy_state.getKineticEnergy()
    assert smd_sim_state.getPotentialEnergy() == smd_sim_copy_state.getPotentialEnergy()
    assert np.array_equal(
        smd_sim.smd_simulation_atom_positions,
        smd_sim_copy.smd_simulation_atom_positions,
    )

    # Check the arguments passed to the OpenMMSMDSimulation class are unchanged by the reset
    assert np.array_equal(smd_sim.smd_atom_indices, smd_sim_copy.smd_atom_indices)
    assert np.array_equal(smd_sim.smd_path, smd_sim_copy.smd_path)
    assert np.array_equal(smd_sim.smd_force_constant, smd_sim_copy.smd_force_constant)

    # Check that the SMD force attached to the simulation is correctly reset
    assert np.array_equal(
        smd_sim.current_smd_force_position, smd_sim_copy.current_smd_force_position
    )
    assert (
        smd_sim.current_smd_force_position_index
        == smd_sim_copy.current_smd_force_position_index
    )
    assert (
        smd_sim.smd_force.getEnergyFunction()
        == smd_sim_copy.smd_force.getEnergyFunction()
    )
    # Class-specific checks
    if type(smd_sim) == OpenMMSMDSimulationAtom:
        assert smd_sim.smd_force.getParticleParameters(
            0
        ) == smd_sim_copy.smd_force.getParticleParameters(0)
        assert smd_sim.smd_force.getNumParticles() == 1
    elif type(smd_sim) == OpenMMSMDSimulationCOM:
        assert smd_sim.smd_force.getGroupParameters(
            0
        ) == smd_sim_copy.smd_force.getGroupParameters(0)
        assert smd_sim.smd_force.getNumBonds() == 1

    # Check other relevant properties of the OpenMMSMDSimulation class match those of the
    # fresh copy after the reset
    assert np.array_equal(
        smd_sim.smd_simulation_atom_positions,
        smd_sim_copy.smd_simulation_atom_positions,
    )


@pytest.mark.parametrize("indices", [TEST_SMD_SINGLE_INDEX, TEST_SMD_MULTIPLE_INDICES])
def test_smd_force_added_to_system(indices):
    """
    Check that the last force to be added to the OpenMM simulation is the SMD force added during
    initialisation of the OpenMMSMDSimulation class, with force group 31.

    :param indices: Indices of atoms to apply the SMD force to (should at least
      test one single index and one set of indices)
    """
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    last_force = smd_sim.simulation.system.getForces()[-1]
    assert type(last_force) == type(smd_sim.smd_force)
    assert last_force.getEnergyFunction() == smd_sim.smd_force.getEnergyFunction()
    assert last_force.getForceGroup() == 31
    # Subclass-specific force type check
    if type(smd_sim) == OpenMMSMDSimulationAtom:
        assert type(smd_sim.smd_force) == CustomExternalForce
    elif type(smd_sim) == OpenMMSMDSimulationCOM:
        assert type(smd_sim.smd_force) == CustomCentroidBondForce


@pytest.mark.parametrize("indices", [TEST_SMD_SINGLE_INDEX, TEST_SMD_MULTIPLE_INDICES])
def test_smd_force_removed_from_system(indices):
    """
    Check that the SMD force is correctly removed from the OpenMM simulation upon calling
    remove_smd_force_from_system().

    :param indices: Indices of atoms to apply the SMD force to (should at least
      test one single index and one set of indices)
    """
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    # Add arbitrary force to system (to test scenario when extra forces added after
    # creation of the SMD class)
    arb_force = CustomExternalForce("0.5 * k * (x)^2")
    arb_force.addGlobalParameter("k", 100.0)
    arb_force.addPerParticleParameter("x")
    smd_sim.simulation.system.addForce(arb_force)

    # Check that the number of forces before and after removal of the SMD
    # force make sense (that only a single SMD force is removed)
    n_forces_before_removal = smd_sim.simulation.system.getNumForces()
    smd_sim.remove_smd_force_from_system()
    n_forces_after_removal = smd_sim.simulation.system.getNumForces()
    assert n_forces_before_removal == n_forces_after_removal + 1

    # Check that none of the energy functions of the remaining system forces
    # match that of the SMD force removed from the system
    system_forces = smd_sim.simulation.system.getForces()
    for force in system_forces:
        try:
            assert force.getEnergyFunction() != smd_sim.smd_force.getEnergyFunction()
        except AttributeError:
            pass
    # TODO: Figure out if it's possible to remove the force constant
    #  associated with the force from global parameters (doesn't seem
    #  to be implemented in OpenMM right now)
    assert smd_sim.simulation.context.getParameter("smd_k") == TEST_SMD_FORCE_CONSTANT


@pytest.mark.parametrize("indices", [TEST_SMD_SINGLE_INDEX, TEST_SMD_MULTIPLE_INDICES])
def test_smd_force_updates_correctly(indices):
    """
    Check that the position of the SMD force is correctly updated upon calling
    update_smd_force_position().

    :param indices: Indices of atoms to apply the SMD force to (should at least
      test one single index and one set of indices)
    """
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    # Choose next force position to be the final position defined by
    # the SMD path
    new_force_position_index = TEST_SMD_PATH.shape[0] - 1
    new_force_position = TEST_SMD_PATH[new_force_position_index]

    # Update the force position and check the relevant class parameters
    # update accordingly
    smd_sim.current_smd_force_position_index = new_force_position_index
    smd_sim.update_smd_force_position()
    assert np.array_equal(smd_sim.current_smd_force_position, new_force_position)

    # Check the subclass-specific force parameters in both the class and the system
    # which should be identical
    n_system_forces = smd_sim.simulation.system.getNumForces()
    if type(smd_sim.smd_force) == CustomExternalForce:

        # OpenMMSMDSimulationAtom force parameters
        index, position = smd_sim.smd_force.getParticleParameters(0)
        assert index == indices
        assert np.array_equal(np.array(position), new_force_position)

        # Force parameters from system
        sys_index, sys_position = smd_sim.simulation.system.getForce(
            n_system_forces - 1
        ).getParticleParameters(0)
        assert sys_index == indices
        assert np.array_equal(np.array(sys_position), new_force_position)

    elif type(smd_sim.smd_force) == CustomCentroidBondForce:

        # OpenMMSMDSimulationCOM force parameters
        _, bond_params = smd_sim.smd_force.getBondParameters(0)
        assert np.array_equal(np.array(bond_params), new_force_position)

        # Force parameters from system
        _, sys_bond_params = smd_sim.simulation.system.getForce(
            n_system_forces - 1
        ).getBondParameters(0)
        assert np.array_equal(np.array(sys_bond_params), new_force_position)


def test_error_for_non_initial_restraint_during_equilibration():
    """
    Check that the SMD simulation throws an error if the user attempts to perform an
    equilibration after updating the position of the SMD force.
    """
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        TEST_SMD_SINGLE_INDEX,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    smd_sim.current_smd_force_position_index = 1
    smd_sim.update_smd_force_position()
    try:
        smd_sim.run_equilibration_with_initial_restraint(n_steps=10)
    except AssertionError:
        pass


@pytest.mark.parametrize("n_structures", [10, 100, 328, 1000])
@pytest.mark.parametrize("interval_ps", [10.0, 100.0])
def test_generate_starting_structures(n_structures, interval_ps):
    """
    Check that the SMD simulation class generates the correct number of starting
    structures in a given time interval, saves them to the correct path, and check
    that the generated files aren't empty.

    :param n_structures: Number of structures to generate
    :param interval_ps: Simulation time interval (in picoseconds) in which to
      generate the structures
    """
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        TEST_SMD_SINGLE_INDEX,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )

    structure_file_prefix = "starting_structure"
    with redirect_stdout(StringIO()) as _:
        with tempfile.TemporaryDirectory() as tmpdir:

            output_path = Path(tmpdir)

            smd_sim.generate_starting_structures(
                interval_ps=interval_ps,
                n_structures=n_structures,
                output_directory=output_path,
                filename_prefix=structure_file_prefix,
                save_smd_force=False,
            )

            # Check that correct number of files are generated
            generated_files = list(output_path.glob(f"{structure_file_prefix}_*.xml"))
            assert len(generated_files) == n_structures

            # Check that filenames are as expected
            expected_filenames = sorted(
                [f"{structure_file_prefix}_{i + 1}.xml" for i in range(n_structures)]
            )
            actual_filenames = sorted(f.name for f in generated_files)
            assert actual_filenames == expected_filenames

            # Check that files aren't empty
            for file in generated_files:
                size = file.stat().st_size
                assert size > 0, f"File {file.name} is unexpectedly empty."


@pytest.mark.parametrize(
    "position_shifts",
    [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
        np.array([2.0, 0.0, 0.0]),
        np.array([1.75, -3.0, 5.263]),
    ],
)
@pytest.mark.parametrize("indices", [TEST_SMD_SINGLE_INDEX, TEST_SMD_MULTIPLE_INDICES])
def test_calculate_smd_forces(position_shifts, indices):
    """
    Test that the function _calculate_smd_forces correctly calculates the SMD forces
    for a given set of positions that is passed to it. As the SMD force is harmonic,
    we expect the force to take the form

    F = - k * (position - smd_force_position)

    This is tested below using the SMD force path given to the simulation, which is
    shifted by some defined by the position_shifts parameter, meaning that we expect
    the forces calculated to take the form

    F = - k * position_shift

    :param position_shifts: Array defining the offset for the positions defined by
      the positions from the test SMD path
    :param indices: Indices of atoms to apply the SMD force to (should at least
      test one single index and one set of indices)
    """
    test_positions = TEST_SMD_PATH + position_shifts
    expected_forces = (
        np.zeros(TEST_SMD_PATH.shape) - TEST_SMD_FORCE_CONSTANT * position_shifts
    )
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    smd_sim._calculate_smd_forces(test_positions)
    assert np.allclose(smd_sim.smd_simulation_forces, expected_forces, atol=1e-16)


@pytest.mark.parametrize(
    "position_shifts",
    [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
        np.array([2.0, 0.0, 0.0]),
        np.array([1.75, -3.0, 5.263]),
    ],
)
@pytest.mark.parametrize("indices", [TEST_SMD_SINGLE_INDEX, TEST_SMD_MULTIPLE_INDICES])
def test_calculate_work_done(position_shifts, indices):
    """
    Check that the work done by the SMD force on the system along the reaction
    coordinate defined by the SMD path is correctly calculated in the function
    _calculate_work_done. Uses the same logic as test_calculate_smd_forces.

    :param position_shifts: Array defining the offset for the positions defined by
      the positions from the test SMD path
    :param indices: Indices of atoms to apply the SMD force to (should at least
      test one single index and one set of indices)
    """
    # TODO: Generalise to curved paths?
    test_positions = TEST_SMD_PATH + position_shifts

    # Calculate displacements of force along test SMD path and
    # check they are all approximately equal
    smd_force_displacements = np.diff(TEST_SMD_PATH, axis=0)
    diff = smd_force_displacements[0]
    assert np.allclose(
        smd_force_displacements,
        np.full(smd_force_displacements.shape, diff),
        atol=1e-16,
    )

    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    smd_sim._calculate_smd_forces(test_positions)

    # Calculate expected work done as a function of time, based on forces and
    # the vector between successive points defining the SMD coordinate. Zeroth
    # value corresponds to work done at t=0 (i.e. zero), so non-zero values
    # start at index 1. SMD paths are straight lines in the current examples,
    # so work done between each step is the same.
    smd_force = smd_sim.smd_simulation_forces[0]
    work_per_step = np.dot(diff, smd_force)
    expected_work_done = np.array(
        [i * work_per_step for i in range(test_positions.shape[0])]
    )

    # Calculate work done using function and check that the values match the
    # expected values
    smd_sim._calculate_work_done()
    assert np.allclose(smd_sim.smd_simulation_work_done, expected_work_done, atol=1e-16)


@pytest.mark.parametrize("indices", [TEST_SMD_SINGLE_INDEX, TEST_SMD_MULTIPLE_INDICES])
def test_save_smd_simulation_data(indices):
    """
    Check that the function save_smd_simulation_data correctly saves the
    data from the specific SMD simulation in the correct format to the
    correct location, and that the data can be subsequently loaded into
    Python, giving the same results as before saving

    :param indices: Indices of atoms to apply the SMD force to (should at least
      test one single index and one set of indices)
    """

    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    with redirect_stdout(StringIO()) as _:
        smd_sim.run_smd()

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir)
        filename = "test_simulation_data.npy"
        file_path = output_path.joinpath(filename)
        smd_sim.save_smd_simulation_data(file_path)
        assert file_path.exists()

        with open(file_path, "rb") as infile:
            loaded_smd_simulation_atom_positions = np.load(infile)
            loaded_smd_simulation_work_done = np.load(infile)

            assert np.array_equal(
                smd_sim.smd_simulation_atom_positions,
                loaded_smd_simulation_atom_positions,
            )
            assert np.array_equal(
                smd_sim.smd_simulation_work_done, loaded_smd_simulation_work_done
            )


@pytest.mark.parametrize("indices", [TEST_SMD_SINGLE_INDEX, TEST_SMD_MULTIPLE_INDICES])
def test_save_general_smd_data(indices):
    """
    Check that the function save_general_smd_data correctly saves the
    data from the SMD simulation in the correct format to the
    correct location, and that the data can be subsequently loaded into
    Python, giving the same results as before saving

    :param indices: Indices of atoms to apply the SMD force to (should at least
      test one single index and one set of indices)
    """
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    with redirect_stdout(StringIO()) as _:
        smd_sim.run_smd()

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir)
        filename = "test_simulation_data.npy"
        file_path = output_path.joinpath(filename)
        smd_sim.save_general_smd_data(file_path)
        assert file_path.exists()

        with open(file_path, "rb") as infile:
            loaded_smd_atom_indices = np.load(infile)
            loaded_smd_path = np.load(infile)
            loaded_smd_force_constant = np.load(infile)
            loaded_temperature = np.load(infile)
            loaded_timestep_ps = np.load(infile)

            assert np.array_equal(smd_sim.smd_atom_indices, loaded_smd_atom_indices)
            assert np.array_equal(smd_sim.smd_path, loaded_smd_path)
            assert np.array_equal(smd_sim.smd_force_constant, loaded_smd_force_constant)
            assert (
                smd_sim.simulation.integrator.getTemperature()._value
                == loaded_temperature
            )
            assert (
                smd_sim.simulation.integrator.getStepSize()._value == loaded_timestep_ps
            )


@pytest.mark.parametrize(
    "positions, masses, com",
    [TEST_COM_TWO_ATOMS, TEST_COM_METHANE, TEST_COM_CIRCLE, TEST_COM_CUBE],
)
def test_calculate_com(positions, masses, com):
    """
    Check that the function calculate_com correctly calculates
    the centre of mass of a set of atoms, given their positions
    and masses.
    """
    calculated_com = calculate_com(positions, masses)
    expected_com = com
    assert np.allclose(calculated_com, expected_com, atol=1e-16)


@pytest.mark.parametrize(
    "positions, masses, com",
    [TEST_COM_TWO_ATOMS, TEST_COM_METHANE, TEST_COM_CIRCLE, TEST_COM_CUBE],
)
def test_calculate_com_smd_simulation_class(positions, masses, com):
    """
    Check that the function _calculate_com correctly calculates
    the centre of mass of the atoms to which the SMD force is
    applied in the OpenMMSMDSimulationCOM class, given their
    positions and masses.
    """
    # Create the simulation and retrieve indices for all atoms
    simulation = build_com_simulation((positions, masses, com))
    indices = np.array([i for i in range(masses.size)])

    # Create the SMD simulation
    smd_sim = OpenMMSMDSimulation.from_simulation(
        simulation,
        indices,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )

    # Calculate the COM using the class function to caand check it against the expected COM
    calculated_com = smd_sim._calculate_com(positions, masses)
    assert np.allclose(calculated_com, com, atol=1e-16)


@pytest.mark.parametrize(
    "positions, masses, com",
    [TEST_COM_TWO_ATOMS, TEST_COM_METHANE, TEST_COM_CIRCLE, TEST_COM_CUBE],
)
def test_calculate_com_trajectory_smd_simulation_class(positions, masses, com):
    """
    Check that the function _calculate_com_trajectory correctly calculates
    the trajectory of the centre of mass of the atoms to which the SMD force is
    applied in the OpenMMSMDSimulationCOM class.
    """
    # Create the simulation and retrieve indices for all atoms
    simulation = build_com_simulation((positions, masses, com))
    indices = np.array([i for i in range(masses.size)])

    # Create the SMD simulation
    smd_sim = OpenMMSMDSimulation.from_simulation(
        simulation,
        indices,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )

    # Manually set the trajectory of atom positions and calculate the
    # corresponding trajectory of the COM
    atom_positions = np.zeros((TEST_SMD_PATH.shape[0], *positions.shape))
    expected_com_array = np.zeros(TEST_SMD_PATH.shape)
    for i in range(TEST_SMD_PATH.shape[0]):
        atom_positions[i] = positions + np.array(
            [TEST_SMD_PATH[i] for j in range(indices.size)]
        )
        expected_com_array[i] = com + TEST_SMD_PATH[i]

    # Set SMD atom positions equal to the trajectory of calculated atom positions
    smd_sim.smd_simulation_atom_positions = atom_positions

    # Calculate COM trajectory using internal function and check the calculated
    # COMs match the predicted COMs
    smd_sim._calculate_com_trajectory()
    assert np.allclose(smd_sim.com_positions, expected_com_array, atol=1e-16)


@pytest.mark.parametrize("pbcs", [True, False])
def test_smd_com_force(pbcs):
    """
    Check that the force produced by the function smd_com_force returns
    a force with the correct properties.
    """
    smd_force = smd_com_force(TEST_SMD_FORCE_CONSTANT, uses_pbcs=pbcs)
    assert type(smd_force) == CustomCentroidBondForce
    assert smd_force.usesPeriodicBoundaryConditions() == pbcs
    assert smd_force.getEnergyFunction() == SMD_FORCE_EXPRESSION_COM
    assert smd_force.getGlobalParameterName(0) == SMD_FORCE_CONSTANT_PARAMETER_NAME
    assert smd_force.getGlobalParameterDefaultValue(0) == TEST_SMD_FORCE_CONSTANT
    assert smd_force.getNumPerBondParameters() == 3
    assert smd_force.getPerBondParameterName(0) == "x0"
    assert smd_force.getPerBondParameterName(1) == "y0"
    assert smd_force.getPerBondParameterName(2) == "z0"
    assert smd_force.getForceGroup() == 31


@pytest.mark.parametrize("pbcs", [True, False])
def test_smd_single_atom_force(pbcs):
    """
    Check that the force produced by the function smd_single_atom_force
    returns a force with the correct properties.
    """
    smd_force = smd_single_atom_force(TEST_SMD_FORCE_CONSTANT, uses_pbcs=pbcs)
    assert type(smd_force) == CustomExternalForce
    assert smd_force.usesPeriodicBoundaryConditions() == pbcs
    if pbcs:
        assert smd_force.getEnergyFunction() == SMD_FORCE_EXPRESSION_ATOM_PERIODIC
    else:
        assert smd_force.getEnergyFunction() == SMD_FORCE_EXPRESSION_ATOM_NONPERIODIC
    assert smd_force.getGlobalParameterName(0) == SMD_FORCE_CONSTANT_PARAMETER_NAME
    assert smd_force.getGlobalParameterDefaultValue(0) == TEST_SMD_FORCE_CONSTANT
    assert smd_force.getNumPerParticleParameters() == 3
    assert smd_force.getPerParticleParameterName(0) == "x0"
    assert smd_force.getPerParticleParameterName(1) == "y0"
    assert smd_force.getPerParticleParameterName(2) == "z0"
    assert smd_force.getForceGroup() == 31


# TODO: Tests single atom case only!
@pytest.mark.parametrize("fc_multiplier", (1.0, 2.02, 0.00750, 4002.8, 5.0003))
def test_calculate_cumulative_work_done(fc_multiplier):
    r"""
    Check that the work done along the reaction coordinate is correctly
    calculated. Use a single atom argon simulation with a Verlet integrator
    to test this to eliminate random kicks induced by a thermostat and
    guarantee reproducibility.

    The expression for work done being tested is

    .. math::
        W(n) = \sum_{i=0}^{n-1} F_i \cdot v_{i} \Delta t

    Where v_{i} is the velocity of the restraint and F_{i} is the force
    applied by the restraint at the ith step. Defining velocity in terms
    of position using the formal definition of the derivative, the
    expression above can be written as

    .. math::
        W(n) = \sum_{i=0}^{n-1} F_i \cdot (x_{i+1} - x_{i})

    where x_{i+1} is the position of the restraint at the (i + 1)th step
    and x_{i} is the position of the restraint at the ith step.

    TEST CASE: a single Argon atom (Ar) that starts at the origin and a
    3 point path defining the positions of the SMD force, starting at the
    origin and increasing along the x-axis in increments of 0.01 nm, with
    a force constant of 100 kJ mol-1 nm-2. The force can be calculated
    using

    .. math::
        F = - k (X_{Ar, i} - x_{i})

    where X_{Ar, i} is the position of the argon atom on the ith step.
    According to the scheme above, the total work done by the time the
    restraint reaches the final point should be

    .. math::
        W(2) = (F_{0} \cdot (x_{1} - x_{0})) + (F_{1} \cdot (x_{2} - x_{1}))
             = (F_{0} \cdot 0.01 nm) + (F_{1} \cdot 0.01 nm)

    Given

    .. math::
        F_{0} = -100.0 kJ mol-1 nm-2 * (0.00 - 0.00) nm = 0.0 kJ mol-1 nm-1

    and

    .. math::
        F_{1} = -100.0 kJ mol-1 nm-2 * (0.00 - 0.01) nm = 1.0 kJ mol-1 nm-1

    Then we expect the work done by step 2 to be

    .. math::
        W(2) = (0.00 kJ mol-1 nm-1 * 0.01 nm) + (1.00 kJ mol-1 nm-1 * 0.01 nm)
             = 0.01 kJ mol-1

    Thus, the work done by the SMD force applied in this simulation should
    be 0.01 kJ mol-1.
    """
    assert TEST_SMD_ARGON_PATH.shape == (3, 3)
    assert TEST_SMD_ARGON_FORCE_CONSTANT == 100.0
    # Create the SMD simulation
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_single_atom_simulation(),
        TEST_SMD_ARGON_INDEX,
        TEST_SMD_ARGON_PATH,
        fc_multiplier * TEST_SMD_ARGON_FORCE_CONSTANT,
    )

    # Run SMD procedure
    with redirect_stdout(StringIO()) as _:
        smd_sim.run_smd()

    # Check values of work done
    assert smd_sim.smd_simulation_work_done[-1] == fc_multiplier * 0.01
    assert np.array_equal(
        smd_sim.smd_simulation_work_done, fc_multiplier * np.array([0.0, 0.0, 0.01])
    )


@pytest.mark.parametrize("apply_pbcs", [True, False])
@pytest.mark.parametrize("indices", [TEST_SMD_SINGLE_INDEX, TEST_SMD_MULTIPLE_INDICES])
def test_load_openmm_state(apply_pbcs, indices):
    """
    Check that the OpenMMSMDSimulation correctly loads the state
    of the system by checking that the velocities loaded are
    correct.
    """
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(apply_pbcs),
        indices,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )

    # Run simulation for a few steps
    smd_sim.run_equilibration_with_initial_restraint(n_steps=1000)
    original_velocities = smd_sim.simulation.context.getState(
        getVelocities=True
    ).getVelocities(asNumpy=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        # Save simulation to file
        output_path = Path(tmpdir)
        filename = "test_velocities.xml"
        file_path = output_path.joinpath(filename)
        smd_sim.save_simulation(output_filepath=file_path, save_state=True)
        assert file_path.exists()

        # Load saved simulation
        loaded_smd_sim = OpenMMSMDSimulation.from_xml_path(
            file_path,
            indices,
            TEST_SMD_PATH,
            TEST_SMD_FORCE_CONSTANT,
        )

        # Retrieve and compare velocities
        loaded_velocities = loaded_smd_sim.simulation.context.getState(
            getVelocities=True
        ).getVelocities(asNumpy=True)
        assert np.allclose(original_velocities, loaded_velocities, rtol=1e-7)
