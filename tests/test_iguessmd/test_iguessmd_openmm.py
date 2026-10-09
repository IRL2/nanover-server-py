"""
Tests for :mod:`nanover.iguessmd.openmm`.

Things to test for class performing iGUESSMD calculations:
- An iGUESSMD simulation can be created either from an existing OpenMM simulation or
  a NanoVer OpenMM XML file [√]
- OMMiGUESSMDSimulation returns OMMiGUESSMDSimulationAtom or OMMiGUESSMDSimulationCOM
  as appropriate [√]
- The PBCs of the loaded simulation are respected by the iGUESSMD force, and the iGUESSMD
  force shares this periodicity [√]
- The iGUESSMD force attaches to the correct atom (dictated by the index/indices passed
  to the class upon creation) [√]
- The positions retrieved from the simulation are correctly unwrapped [√]
- The simulation can be reset correctly, with all attributes returning to the same
  state as immediately after the creation of the class itself [√]
- iGUESSMD force is correctly added to the system [√]
- iGUESSMD force can be correctly removed from the system [√]
- iGUESSMD force position is correctly updated [√]
- Running the equilibration with the initial restraint throws an error correctly
  when run with the iGUESSMD force not located at the initial position [√]
- The iGUESSMD simulation can be saved correctly, with or without the iGUESSMD force [√]
- If the iGUESSMD simulation is being loaded from a NanoVer OpenMM XML file created
  via one of the OMMiGUESSMDSimulation classes and contains an iGUESSMD force already, this
  iGUESSMD force is correctly loaded and matches the expected force constant specified
  when creating the class [√]
- The class can generate the correct number of starting structures in the specified
  time interval, and that these are saved to the correct location [√]
- Running an iGUESSMD simulation produces reasonable results for the cumulative work done [√]
- _calculate_iguessmd_forces works as expected [√]
- _calculate_work_done works as expected [√]
- Simulation data is saved in the correct format to the correct location, and can be
  subsequently loaded back into python correctly [√]
- General iGUESSMD data is saved in the correct format to the correct location, and can be
  subsequently loaded back into python correctly [√]
- The COM of a specified group of atoms is correctly calculated [√]
- For OMMiGUESSMDSimulationCOM, the COM of the specified atoms is correctly calculated [√]
- For OMMiGUESSMDSimulationCOM, the trajectory of the COM of the specified atoms is
  correctly calculated [√]
- iguessmd_com_force works as expected [√]
- iguessmd_single_atom_force works as expected [√]
- OMMiGUESSMDSimulation correctly loads the state of a simulation [√]
- OMMiGUESSMDSimulation correctly calculates the unit tangents of the iGUESSMD path [√]
- Different parallel and perpendicular force constants can be given and correctly applied
  to the simulation [√]
- Different parallel and perpendicular force constants can be correctly saved and loaded
  from the general iGUESSMD data file [√]
- Different parallel and perpendicular force constants can be correctly saved and loaded
  from an iGUESSMD simulation saved to a NanoVer OpenMM XML file [√]
"""

import tempfile
from io import StringIO
from itertools import product

import pytest
from contextlib import redirect_stdout

import openmm as mm
from mypy.main import fail
from openmm import app
from openmm.unit import (
    kelvin,
    picosecond,
    femtosecond,
    nanometer,
)

from nanover.iguessmd.openmm import *

from .iguessmd_test_utilities import define_circular_path

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

# Test parameters for OMMiGUESSMDSimulation
TEST_iGUESSMD_SINGLE_INDEX = np.array(0)
TEST_iGUESSMD_MULTIPLE_INDICES = np.array([0, 1, 2, 3])
TEST_iGUESSMD_INDICES = [TEST_iGUESSMD_SINGLE_INDEX, TEST_iGUESSMD_MULTIPLE_INDICES]
TEST_iGUESSMD_LINEAR_PATH = np.array(
    [np.linspace(0.05, 1.05, 101), np.zeros(101), np.zeros(101)]
).transpose()
TEST_iGUESSMD_LINEAR_PATH_TANGENTS = np.array(
    [np.ones(101), np.zeros(101), np.zeros(101)]
).transpose()
TEST_iGUESSMD_CIRCULAR_PATH, TEST_iGUESSMD_CIRCULAR_PATH_TANGENTS = define_circular_path(10)
TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL = 3011.0
TEST_iGUESSMD_FORCE_CONSTANT_PAR_PERP = np.array([3011.0, 301.1])
TEST_iGUESSMD_FORCE_CONSTANTS = [
    TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
    TEST_iGUESSMD_FORCE_CONSTANT_PAR_PERP,
]

TEST_iGUESSMD_ARGON_INDEX = np.array(0)
TEST_iGUESSMD_ARGON_PATH = np.array(
    [np.linspace(0.00, 0.02, 3), np.zeros(3), np.zeros(3)]
).transpose()
TEST_iGUESSMD_ARGON_FORCE_CONSTANT = 100.0

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

TEST_BOOLS = [True, False]


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
def make_basic_iguessmd_simulation_with_iguessmd_force_xml(
    tmp_path,
    request,
):
    atom_indices, force_constant = request.param
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        atom_indices,
        TEST_iGUESSMD_LINEAR_PATH,
        force_constant,
    )
    xml_path = tmp_path / "basic_iguessmd_simulation.xml"
    iguessmd_sim.save_simulation(xml_path, save_state=True, save_iguessmd_force=True)
    return xml_path, atom_indices, force_constant


@pytest.fixture
def make_basic_iguessmd_simulation_without_iguessmd_force_xml(
    tmp_path,
    request,
):
    atom_indices, force_constant = request.param
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        atom_indices,
        TEST_iGUESSMD_LINEAR_PATH,
        force_constant,
    )
    xml_path = tmp_path / "basic_iguessmd_simulation.xml"
    iguessmd_sim.save_simulation(xml_path, save_state=True, save_iguessmd_force=False)
    return xml_path, atom_indices, force_constant


@pytest.mark.parametrize(
    "force_constant, atom_indices",
    product(TEST_iGUESSMD_FORCE_CONSTANTS, TEST_iGUESSMD_INDICES),
)
def test_load_iguessmd_sim_from_simulation(force_constant, atom_indices):
    """
    Test that an OMMiGUESSMDSimulation can be correctly loaded from an OpenMM simulation.
    """
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        atom_indices,
        TEST_iGUESSMD_LINEAR_PATH,
        force_constant,
    )
    assert iguessmd_sim
    assert iguessmd_sim.simulation
    assert np.array_equal(iguessmd_sim.iguessmd_path, TEST_iGUESSMD_LINEAR_PATH)
    assert np.array_equal(iguessmd_sim.iguessmd_atom_indices, atom_indices)
    assert np.array_equal(iguessmd_sim.iguessmd_force_constant, force_constant)


@pytest.mark.parametrize(
    "force_constant, atom_indices",
    product(TEST_iGUESSMD_FORCE_CONSTANTS, TEST_iGUESSMD_INDICES),
)
def test_load_iguessmd_sim_from_xml_path(
    make_basic_simulation_xml, force_constant, atom_indices
):
    """
    Test that an OMMiGUESSMDSimulation can be correctly loaded from a NanoVer OpenMM XML file.
    """
    iguessmd_sim = OMMiGUESSMDSimulation.from_xml_path(
        make_basic_simulation_xml,
        atom_indices,
        TEST_iGUESSMD_LINEAR_PATH,
        force_constant,
    )
    assert iguessmd_sim
    assert iguessmd_sim.xml_path == make_basic_simulation_xml
    assert iguessmd_sim.simulation
    assert np.array_equal(iguessmd_sim.iguessmd_path, TEST_iGUESSMD_LINEAR_PATH)
    assert np.array_equal(iguessmd_sim.iguessmd_atom_indices, atom_indices)
    assert np.array_equal(iguessmd_sim.iguessmd_force_constant, force_constant)


@pytest.mark.parametrize(
    "make_basic_iguessmd_simulation_with_iguessmd_force_xml",
    product(TEST_iGUESSMD_INDICES, TEST_iGUESSMD_FORCE_CONSTANTS),
    indirect=True,
)
def test_load_iguessmd_simulation_with_iguessmd_force_from_xml_path(
    make_basic_iguessmd_simulation_with_iguessmd_force_xml,
):
    """
    Check that when an input file containing an iGUESSMD force is passed to the
    OMMiGUESSMDSimulation class, the iGUESSMD force is loaded correctly from the file using
    check_for_existing_iguessmd_force(), and that the parameters for the iGUESSMD force match
    those that are passed via the file.
    """
    xml_path, atom_indices, force_constant = (
        make_basic_iguessmd_simulation_with_iguessmd_force_xml
    )
    with redirect_stdout(StringIO()) as _:
        iguessmd_sim = OMMiGUESSMDSimulation.from_xml_path(
            xml_path,
            atom_indices,
            TEST_iGUESSMD_LINEAR_PATH,
            force_constant,
        )
        assert iguessmd_sim.loaded_iguessmd_force_from_sim
        assert np.array_equal(iguessmd_sim.iguessmd_force_constant, force_constant)


@pytest.mark.parametrize(
    "make_basic_iguessmd_simulation_without_iguessmd_force_xml",
    product(TEST_iGUESSMD_INDICES, TEST_iGUESSMD_FORCE_CONSTANTS),
    indirect=True,
)
def test_load_iguessmd_simulation_without_iguessmd_force_from_xml_path(
    make_basic_iguessmd_simulation_without_iguessmd_force_xml,
):
    """
    Check that when an xml input file is saved from an OMMiGUESSMDSimulationAtom class
    without the iGUESSMD force, the iGUESSMD force is not loaded from the file.
    """
    xml_path, atom_indices, force_constant = (
        make_basic_iguessmd_simulation_without_iguessmd_force_xml
    )
    iguessmd_sim = OMMiGUESSMDSimulation.from_xml_path(
        xml_path,
        atom_indices,
        TEST_iGUESSMD_LINEAR_PATH,
        force_constant,
    )
    assert not iguessmd_sim.loaded_iguessmd_force_from_sim


@pytest.mark.parametrize(
    "indices, sim_type",
    [
        (TEST_iGUESSMD_SINGLE_INDEX, OMMiGUESSMDSimulationAtom),
        (TEST_iGUESSMD_MULTIPLE_INDICES, OMMiGUESSMDSimulationCOM),
    ],
)
def test_return_correct_iguessmd_sim_type(indices, sim_type):
    """
    Check that the OMMiGUESSMDSimulation class returns the correct subclass depending on the
    number of indices that are passed to it (one for OMMiGUESSMDSimulationAtom, more than one
    for OMMiGUESSMDSimulationCOM).

    :param indices: Indices of atoms to apply the iGUESSMD force to (should at least
      test one single index and one set of indices)
    :param sim_type: Type of simulation to expect for the indices given
    """
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_iGUESSMD_LINEAR_PATH,
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
    )
    assert type(iguessmd_sim) == sim_type


@pytest.mark.parametrize(
    "apply_pbcs, indices", product(TEST_BOOLS, TEST_iGUESSMD_INDICES)
)
def test_simulation_pbcs_are_respected(apply_pbcs, indices):
    """
    Check that the periodic boundary conditions of the OpenMMSimulation passed to the
    OMMiGUESSMDSimulation class are respected (i.e. the PBCs of the iGUESSMD simulation match
    those of the OpenMM simulation).

    :param apply_pbcs: Boolean value indicating whether to apply PBCs to the simulation
    :param indices: Indices of atoms to apply the iGUESSMD force to (should at least
      test one single index and one set of indices)
    """
    sim = build_basic_simulation(pbcs=apply_pbcs)
    uses_pbcs = sim.system.usesPeriodicBoundaryConditions()
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        sim,
        indices,
        TEST_iGUESSMD_LINEAR_PATH,
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
    )

    # Check PBCs are read correctly by iGUESSMD simulation
    assert iguessmd_sim._sim_uses_pbcs == uses_pbcs
    assert iguessmd_sim.simulation.system.usesPeriodicBoundaryConditions() == uses_pbcs


@pytest.mark.parametrize(
    "index, force_constant",
    product(
        [np.array(0), np.array(1), np.array(4), np.array(5), np.array(7)],
        TEST_iGUESSMD_FORCE_CONSTANTS,
    ),
)
def test_iguessmd_force_attaches_to_correct_atom(index, force_constant):
    """
    Check that the iGUESSMD force is attached to the correct atom when a single index is passed.
    Should use the OMMiGUESSMDSimulationAtom class, with only one CustomExternalForce.

    :param index: Indices of atoms to apply the iGUESSMD force to (should be arrays containing
      a single index)
    """
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        index,
        TEST_iGUESSMD_LINEAR_PATH,
        force_constant,
    )
    # Attaches force to single atom, so index of atom within force is zero
    p_index, p_params = iguessmd_sim.iguessmd_force.getParticleParameters(0)
    assert p_index == index
    assert np.array_equal(
        np.array(p_params),
        np.array([*TEST_iGUESSMD_LINEAR_PATH[0], *TEST_iGUESSMD_LINEAR_PATH_TANGENTS[0]]),
    )


@pytest.mark.parametrize(
    "indices, force_constant",
    product(
        [
            np.array([0, 1, 2, 3]),
            np.array([1, 2, 3, 4]),
            np.array([0, 1, 4, 5]),
            np.array([1, 3, 4, 7]),
        ],
        TEST_iGUESSMD_FORCE_CONSTANTS,
    ),
)
def test_iguessmd_force_attaches_to_correct_atoms(indices, force_constant):
    """
    Check that the iGUESSMD force attaches to the correct atoms when an array of indices is passed.
    Should use the OMMiGUESSMDSimulationCOM class, with only one CustomCentroidBondForce.

    :param indices: Indices of atoms to apply the iGUESSMD force to (should be arrays of multiple
      indices)
    """
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_iGUESSMD_LINEAR_PATH,
        force_constant,
    )
    # Only one centroid force added, index of force is zero
    p_indices, _ = iguessmd_sim.iguessmd_force.getGroupParameters(0)
    assert np.array_equal(np.array(p_indices), indices)


@pytest.mark.parametrize("indices", TEST_iGUESSMD_INDICES)
def test_reset(indices):
    """
    Check that all the attributes of the OMMiGUESSMDSimulation subclasses are reset to their initial
    state by the .reset() function of the class.

    :param indices: Indices of atoms to apply the iGUESSMD force to (should at least
      test one single index and one set of indices)
    """
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_iGUESSMD_LINEAR_PATH,
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
    )

    with redirect_stdout(StringIO()) as _:
        iguessmd_sim.run_iguessmd()
        iguessmd_sim.reset()

    # Create a fresh copy of the iGUESSMD simulation
    iguessmd_sim_copy = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_iGUESSMD_LINEAR_PATH,
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
    )

    # Check that the attributes that were created during the iGUESSMD simulation
    # are no longer present in the class
    try:
        assert (
            iguessmd_sim.iguessmd_simulation_forces
            or iguessmd_sim_copy.iguessmd_simulation_work_done
        )
    except AttributeError:
        pass

    # Check relevant observables from the simulation against the fresh simulation
    iguessmd_sim_state = iguessmd_sim.simulation.context.getState(
        getPositions=True, getVelocities=True, getForces=True, getEnergy=True
    )
    iguessmd_sim_copy_state = iguessmd_sim_copy.simulation.context.getState(
        getPositions=True, getVelocities=True, getForces=True, getEnergy=True
    )
    assert np.array_equal(
        iguessmd_sim_state.getPositions(asNumpy=True),
        iguessmd_sim_copy_state.getPositions(asNumpy=True),
    )
    assert np.array_equal(
        iguessmd_sim_state.getVelocities(asNumpy=True),
        iguessmd_sim_copy_state.getVelocities(asNumpy=True),
    )
    assert np.array_equal(
        iguessmd_sim_state.getForces(asNumpy=True),
        iguessmd_sim_copy_state.getForces(asNumpy=True),
    )
    assert (
        iguessmd_sim_state.getKineticEnergy()
        == iguessmd_sim_copy_state.getKineticEnergy()
    )
    assert (
        iguessmd_sim_state.getPotentialEnergy()
        == iguessmd_sim_copy_state.getPotentialEnergy()
    )
    assert np.array_equal(
        iguessmd_sim.iguessmd_simulation_atom_positions,
        iguessmd_sim_copy.iguessmd_simulation_atom_positions,
    )

    # Check the arguments passed to the OMMiGUESSMDSimulation class are unchanged by the reset
    assert np.array_equal(
        iguessmd_sim.iguessmd_atom_indices, iguessmd_sim_copy.iguessmd_atom_indices
    )
    assert np.array_equal(iguessmd_sim.iguessmd_path, iguessmd_sim_copy.iguessmd_path)
    assert np.array_equal(
        iguessmd_sim.iguessmd_force_constant, iguessmd_sim_copy.iguessmd_force_constant
    )

    # Check that the iGUESSMD force attached to the simulation is correctly reset
    assert np.array_equal(
        iguessmd_sim.current_iguessmd_force_position,
        iguessmd_sim_copy.current_iguessmd_force_position,
    )
    assert (
        iguessmd_sim.current_iguessmd_force_position_index
        == iguessmd_sim_copy.current_iguessmd_force_position_index
    )
    assert (
        iguessmd_sim.iguessmd_force.getEnergyFunction()
        == iguessmd_sim_copy.iguessmd_force.getEnergyFunction()
    )
    # Class-specific checks
    if type(iguessmd_sim) == OMMiGUESSMDSimulationAtom:
        assert iguessmd_sim.iguessmd_force.getParticleParameters(
            0
        ) == iguessmd_sim_copy.iguessmd_force.getParticleParameters(0)
        assert iguessmd_sim.iguessmd_force.getNumParticles() == 1
    elif type(iguessmd_sim) == OMMiGUESSMDSimulationCOM:
        assert iguessmd_sim.iguessmd_force.getGroupParameters(
            0
        ) == iguessmd_sim_copy.iguessmd_force.getGroupParameters(0)
        assert iguessmd_sim.iguessmd_force.getNumBonds() == 1

    # Check other relevant properties of the OMMiGUESSMDSimulation class match those of the
    # fresh copy after the reset
    assert np.array_equal(
        iguessmd_sim.iguessmd_simulation_atom_positions,
        iguessmd_sim_copy.iguessmd_simulation_atom_positions,
    )


@pytest.mark.parametrize("indices", TEST_iGUESSMD_INDICES)
def test_iguessmd_force_added_to_system(indices):
    """
    Check that the last force to be added to the OpenMM simulation is the iGUESSMD force added during
    initialisation of the OMMiGUESSMDSimulation class, with force group 31.

    :param indices: Indices of atoms to apply the iGUESSMD force to (should at least
      test one single index and one set of indices)
    """
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_iGUESSMD_LINEAR_PATH,
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
    )
    last_force = iguessmd_sim.simulation.system.getForces()[-1]
    assert type(last_force) == type(iguessmd_sim.iguessmd_force)
    assert (
        last_force.getEnergyFunction()
        == iguessmd_sim.iguessmd_force.getEnergyFunction()
    )
    assert last_force.getForceGroup() == 31
    # Subclass-specific force type check
    if type(iguessmd_sim) == OMMiGUESSMDSimulationAtom:
        assert type(iguessmd_sim.iguessmd_force) == CustomExternalForce
    elif type(iguessmd_sim) == OMMiGUESSMDSimulationCOM:
        assert type(iguessmd_sim.iguessmd_force) == CustomCentroidBondForce


@pytest.mark.parametrize(
    "indices, force_constant",
    product(TEST_iGUESSMD_INDICES, TEST_iGUESSMD_FORCE_CONSTANTS),
)
def test_iguessmd_force_removed_from_system(indices, force_constant):
    """
    Check that the iGUESSMD force is correctly removed from the OpenMM simulation upon calling
    remove_iguessmd_force_from_system().

    :param indices: Indices of atoms to apply the iGUESSMD force to (should at least
      test one single index and one set of indices)
    """
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_iGUESSMD_LINEAR_PATH,
        force_constant,
    )
    # Add arbitrary force to system (to test scenario when extra forces added after
    # creation of the iGUESSMD class)
    arb_force = CustomExternalForce("0.5 * k * (x)^2")
    arb_force.addGlobalParameter("k", 100.0)
    arb_force.addPerParticleParameter("x")
    iguessmd_sim.simulation.system.addForce(arb_force)

    # Check that the number of forces before and after removal of the iGUESSMD
    # force make sense (that only a single iGUESSMD force is removed)
    n_forces_before_removal = iguessmd_sim.simulation.system.getNumForces()
    iguessmd_sim.remove_iguessmd_force_from_system()
    n_forces_after_removal = iguessmd_sim.simulation.system.getNumForces()
    assert n_forces_before_removal == n_forces_after_removal + 1

    # Check that none of the energy functions of the remaining system forces
    # match that of the iGUESSMD force removed from the system
    system_forces = iguessmd_sim.simulation.system.getForces()
    for force in system_forces:
        try:
            assert (
                force.getEnergyFunction()
                != iguessmd_sim.iguessmd_force.getEnergyFunction()
            )
        except AttributeError:
            pass

    # Cannot currently fully remove global parameters from OpenMM
    # simulation context, so check these parameters and their values
    if isinstance(force_constant, np.ndarray):
        assert (
            iguessmd_sim.simulation.context.getParameter("smd_k_par")
            == force_constant[0]
        )
        assert (
            iguessmd_sim.simulation.context.getParameter("smd_k_perp")
            == force_constant[1]
        )
    else:
        assert (
            iguessmd_sim.simulation.context.getParameter("smd_k_par") == force_constant
        )


@pytest.mark.parametrize("indices", TEST_iGUESSMD_INDICES)
def test_iguessmd_force_updates_correctly(indices):
    """
    Check that the position of the iGUESSMD force is correctly updated upon calling
    update_iguessmd_force_position().

    :param indices: Indices of atoms to apply the iGUESSMD force to (should at least
      test one single index and one set of indices)
    """
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_iGUESSMD_LINEAR_PATH,
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
    )
    # Choose next force position to be the final position defined by
    # the iGUESSMD path
    new_force_position_index = TEST_iGUESSMD_LINEAR_PATH.shape[0] - 1
    new_force_position = TEST_iGUESSMD_LINEAR_PATH[new_force_position_index]
    new_force_tangent = TEST_iGUESSMD_LINEAR_PATH_TANGENTS[new_force_position_index]

    # Update the force position and check the relevant class parameters
    # update accordingly
    iguessmd_sim.current_iguessmd_force_position_index = new_force_position_index
    iguessmd_sim.update_iguessmd_force_position()
    assert np.array_equal(
        iguessmd_sim.current_iguessmd_force_position, new_force_position
    )
    assert np.array_equal(
        iguessmd_sim.current_iguessmd_force_tangent, new_force_tangent
    )

    # Check the subclass-specific force parameters in both the class and the system
    # which should be identical
    n_system_forces = iguessmd_sim.simulation.system.getNumForces()
    if type(iguessmd_sim.iguessmd_force) == CustomExternalForce:
        # OMMiGUESSMDSimulationAtom force parameters
        index, position = iguessmd_sim.iguessmd_force.getParticleParameters(0)
        assert index == indices
        assert np.array_equal(
            np.array(position), np.array([*new_force_position, *new_force_tangent])
        )

        # Force parameters from system
        sys_index, sys_position = iguessmd_sim.simulation.system.getForce(
            n_system_forces - 1
        ).getParticleParameters(0)
        assert sys_index == indices
        assert np.array_equal(
            np.array(sys_position), np.array([*new_force_position, *new_force_tangent])
        )

    elif type(iguessmd_sim.iguessmd_force) == CustomCentroidBondForce:
        # OMMiGUESSMDSimulationCOM force parameters
        _, bond_params = iguessmd_sim.iguessmd_force.getBondParameters(0)
        assert np.array_equal(
            np.array(bond_params), np.array([*new_force_position, *new_force_tangent])
        )

        # Force parameters from system
        _, sys_bond_params = iguessmd_sim.simulation.system.getForce(
            n_system_forces - 1
        ).getBondParameters(0)
        assert np.array_equal(
            np.array(sys_bond_params),
            np.array([*new_force_position, *new_force_tangent]),
        )


def test_error_for_non_initial_restraint_during_equilibration():
    """
    Check that the iGUESSMD simulation throws an error if the user attempts to perform an
    equilibration after updating the position of the iGUESSMD force.
    """
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        TEST_iGUESSMD_SINGLE_INDEX,
        TEST_iGUESSMD_LINEAR_PATH,
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
    )
    iguessmd_sim.current_iguessmd_force_position_index = 1
    iguessmd_sim.update_iguessmd_force_position()
    with pytest.raises(AssertionError):
        iguessmd_sim.run_equilibration_with_initial_restraint(n_steps=10)


@pytest.mark.parametrize("n_structures", [10, 100, 328, 1000])
@pytest.mark.parametrize("interval_ps", [10.0, 100.0])
def test_generate_starting_structures(n_structures, interval_ps):
    """
    Check that the iGUESSMD simulation class generates the correct number of starting
    structures in a given time interval, saves them to the correct path, and check
    that the generated files aren't empty.

    :param n_structures: Number of structures to generate
    :param interval_ps: Simulation time interval (in picoseconds) in which to
      generate the structures
    """
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        TEST_iGUESSMD_SINGLE_INDEX,
        TEST_iGUESSMD_LINEAR_PATH,
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
    )

    structure_file_prefix = "starting_structure"
    with redirect_stdout(StringIO()) as _:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir)

            iguessmd_sim.generate_starting_structures(
                interval_ps=interval_ps,
                n_structures=n_structures,
                output_directory=output_path,
                filename_prefix=structure_file_prefix,
                save_iguessmd_force=False,
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
    "name, path, tangents",
    [
        ("linear", TEST_iGUESSMD_LINEAR_PATH, TEST_iGUESSMD_LINEAR_PATH_TANGENTS),
        ("circular", TEST_iGUESSMD_CIRCULAR_PATH, TEST_iGUESSMD_CIRCULAR_PATH_TANGENTS),
    ],
)
def test_calculate_iguessmd_path_tangents(name, path, tangents):
    """
    Test that checks whether the normalised tangent vectors of an
    iGUESSMD path are calculated correctly, using the forward difference
    approximation and with the final tangent vector equal to the
    penultimate tangent vector.
    """
    # Create iGUESSMD simulation
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        TEST_iGUESSMD_SINGLE_INDEX,
        path,
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
    )

    # Calculate tangent vectors
    iguessmd_sim.calculate_iguessmd_path_tangents()

    # Check that calculated tangents are as expected
    assert iguessmd_sim.iguessmd_path_tangents.shape == tangents.shape
    assert np.allclose(
        iguessmd_sim.iguessmd_path_tangents, tangents, rtol=1e-8
    )

    # For circular path, check calculated tangent vectors obey symmetry
    # of circle (excluding final tangent)
    if name == "circular":
        n_points = path.shape[0]
        assert np.allclose(
            iguessmd_sim.iguessmd_path_tangents[: int((n_points - 1) / 2)],
            -iguessmd_sim.iguessmd_path_tangents[int((n_points - 1) / 2): -1],
            rtol=1e-8,
        )


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
@pytest.mark.parametrize(
    "indices, force_constant",
    product(TEST_iGUESSMD_INDICES, TEST_iGUESSMD_FORCE_CONSTANTS),
)
def test_calculate_iguessmd_forces(position_shifts, indices, force_constant):
    """
    Test that the function _calculate_iguessmd_forces correctly calculates the iGUESSMD forces
    for a given set of positions that is passed to it. As the iGUESSMD force is harmonic,
    we expect the force to take the general form

    F = - k_par * (position - iguessmd_force_position)_par
          - k_perp * (position - iguessmd_force_position)_perp

    which reduces to

    F = - k * (position - iguessmd_force_position)

    in the case that the potential is spherically symmetric (the parallel and perpendicular force
    constants are equal)

    This is tested below using the iGUESSMD force path given to the simulation, which is
    shifted by some defined by the position_shifts parameter, meaning that we expect
    the forces calculated to take the form

    F = - k_par * position_shift_par - k_perp * position_shift_perp

    :param position_shifts: Array defining the offset for the positions defined by
      the positions from the test iGUESSMD path
    :param indices: Indices of atoms to apply the iGUESSMD force to (should at least
      test one single index and one set of indices)
    """
    # Test atom positions and atom-interaction centre displacement vectors
    test_positions = TEST_iGUESSMD_LINEAR_PATH + position_shifts
    displacements = test_positions - TEST_iGUESSMD_LINEAR_PATH

    # Check for different parallel and perpendicular force constants
    if isinstance(force_constant, np.ndarray):
        fc_par = force_constant[0]
        fc_perp = force_constant[1]
    else:
        fc_par = fc_perp = force_constant

    # Calculate normalised tangent vectors along RC
    tangents = np.diff(TEST_iGUESSMD_LINEAR_PATH, axis=0)
    tangents = np.array([*tangents, tangents[-1]])
    tangents /= np.linalg.norm(tangents, axis=1, keepdims=True)

    # Calculate parallel and perpendicular components of displacement vectors
    displacements_par = (
        np.linalg.vecdot(displacements, tangents, axis=1)[:, None] * tangents
    )
    displacements_perp = displacements - displacements_par

    # Calculate force components and total force
    expected_forces_par = -fc_par * displacements_par
    expected_forces_perp = -fc_perp * displacements_perp
    expected_forces = expected_forces_par + expected_forces_perp

    # Create iGUESSMD simulation and pass test positions to
    # function to calculate forces
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_iGUESSMD_LINEAR_PATH,
        force_constant,
    )
    iguessmd_sim._calculate_iguessmd_forces(test_positions)

    # Check calculated forces are equal to expected forces
    assert np.allclose(
        iguessmd_sim.iguessmd_simulation_forces, expected_forces, atol=1e-16
    )


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
@pytest.mark.parametrize("indices", TEST_iGUESSMD_INDICES)
def test_calculate_work_done(position_shifts, indices):
    """
    Check that the work done by the iGUESSMD force on the system along the reaction
    coordinate defined by the iGUESSMD path is correctly calculated in the function
    _calculate_work_done. Uses the same logic as test_calculate_iguessmd_forces.

    :param position_shifts: Array defining the offset for the positions defined by
      the positions from the test iGUESSMD path
    :param indices: Indices of atoms to apply the iGUESSMD force to (should at least
      test one single index and one set of indices)
    """
    # TODO: Generalise to curved paths?
    test_positions = TEST_iGUESSMD_LINEAR_PATH + position_shifts

    # Calculate displacements of force along test iGUESSMD path and
    # check they are all approximately equal
    iguessmd_force_displacements = np.diff(TEST_iGUESSMD_LINEAR_PATH, axis=0)
    diff = iguessmd_force_displacements[0]
    assert np.allclose(
        iguessmd_force_displacements,
        np.full(iguessmd_force_displacements.shape, diff),
        atol=1e-16,
    )

    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_iGUESSMD_LINEAR_PATH,
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
    )
    iguessmd_sim._calculate_iguessmd_forces(test_positions)

    # Calculate expected work done as a function of time, based on forces and
    # the vector between successive points defining the iGUESSMD coordinate. Zeroth
    # value corresponds to work done at t=0 (i.e. zero), so non-zero values
    # start at index 1. iGUESSMD paths are straight lines in the current examples,
    # so work done between each step is the same.
    iguessmd_force = iguessmd_sim.iguessmd_simulation_forces[0]
    work_per_step = np.dot(diff, iguessmd_force)
    expected_work_done = np.array(
        [i * work_per_step for i in range(test_positions.shape[0])]
    )

    # Calculate work done using function and check that the values match the
    # expected values
    iguessmd_sim._calculate_work_done()
    assert np.allclose(
        iguessmd_sim.iguessmd_simulation_work_done, expected_work_done, atol=1e-16
    )


@pytest.mark.parametrize(
    "indices, force_constants, every_nth, dtypes",
    product(
        TEST_iGUESSMD_INDICES,
        TEST_iGUESSMD_FORCE_CONSTANTS,
        [1, 2, 5, 23],
        [np.float16, np.float32, np.float64],
    ),
)
def test_save_iguessmd_simulation_data(indices, force_constants, every_nth, dtypes):
    """
    Check that the iGUESSMD simulation data can be saved to a file with the
    expected behaviour:

        - the function save_iguessmd_simulation_data correctly saves the
          data from the specific iGUESSMD simulation in the correct format
          to the correct location
        - the data can be subsequently loaded into Python
        - the data is saved as the desired datatype
        - the saved data contain the same results as contained in
          the simulation instance before saving

    :param indices: Indices of atoms to apply the iGUESSMD force to (should at least
      test one single index and one set of indices)
    """
    # Create the simulation
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_iGUESSMD_LINEAR_PATH,
        force_constants,
    )

    # Run iGUESSMD simulation to generate example data
    with redirect_stdout(StringIO()) as _:
        iguessmd_sim.run_iguessmd()

    with tempfile.TemporaryDirectory() as tmpdir:

        # Define output file
        output_path = Path(tmpdir)
        filename = "test_simulation_data.npz"
        file_path = output_path.joinpath(filename)

        # Save as specified dtype (currently arrays are dtype float64 internally)
        if indices.size > 1:
            iguessmd_sim.save_iguessmd_simulation_data(
                file_path, every_nth=every_nth, atom_positions_dtype=dtypes, work_done_dtype=dtypes, com_positions_dtype=dtypes
            )
        else:
            iguessmd_sim.save_iguessmd_simulation_data(
                file_path, every_nth=every_nth, atom_positions_dtype=dtypes, work_done_dtype=dtypes
            )

        # Check the outfile exists
        assert file_path.exists()

        # Load saved data
        iguessmd_simulation_data = dict(np.load(file_path))

        # Check that the data timestep saved correctly
        assert iguessmd_simulation_data["data_timestep_ps"] == every_nth * iguessmd_sim.simulation.integrator.getStepSize()._value


        # Check that arrays save correctly (including dtype)
        assert iguessmd_simulation_data["iguessmd_simulation_atom_positions"].dtype == dtypes
        assert np.array_equal(
            get_every_nth(
                iguessmd_sim.iguessmd_simulation_atom_positions,
                axis=0,
                every_nth=every_nth,
                include_end=False
            ).astype(dtypes),
            iguessmd_simulation_data["iguessmd_simulation_atom_positions"],
        )
        assert iguessmd_simulation_data["iguessmd_simulation_work_done"].dtype == dtypes
        assert np.array_equal(
            get_every_nth(
                iguessmd_sim.iguessmd_simulation_work_done,
                axis=0,
                every_nth=every_nth,
                include_end=False
            ).astype(dtypes),
            iguessmd_simulation_data["iguessmd_simulation_work_done"],
        )
        if indices.size > 1:
            assert iguessmd_simulation_data["iguessmd_com_positions"].dtype == dtypes
            assert np.array_equal(
                get_every_nth(
                    iguessmd_sim.iguessmd_com_positions,
                    axis=0,
                    every_nth=every_nth,
                    include_end=False
                ).astype(dtypes),
                iguessmd_simulation_data["iguessmd_com_positions"],
            )


@pytest.mark.parametrize(
    "indices, force_constants",
    product(TEST_iGUESSMD_INDICES, TEST_iGUESSMD_FORCE_CONSTANTS),
)
def test_save_general_iguessmd_data(indices, force_constants):
    """
    Check that the function save_general_iguessmd_data correctly saves the
    data from the iGUESSMD simulation in the correct format to the
    correct location, and that the data can be subsequently loaded into
    Python, giving the same results as before saving

    :param indices: Indices of atoms to apply the iGUESSMD force to (should at least
      test one single index and one set of indices)
    """
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_iGUESSMD_LINEAR_PATH,
        force_constants,
    )
    with redirect_stdout(StringIO()) as _:
        iguessmd_sim.run_iguessmd()

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir)
        filename = "test_simulation_data.npz"
        file_path = output_path.joinpath(filename)
        iguessmd_sim.save_general_iguessmd_data(file_path)
        assert file_path.exists()

        # Load saved data
        general_iguessmd_data = dict(np.load(file_path))

        assert np.array_equal(
            iguessmd_sim.iguessmd_atom_indices, general_iguessmd_data["iguessmd_atom_indices"]
        )
        assert np.array_equal(iguessmd_sim.iguessmd_path, general_iguessmd_data["iguessmd_path"])
        assert np.array_equal(
            iguessmd_sim.iguessmd_force_constant, general_iguessmd_data["iguessmd_force_constant"]
        )
        assert (
                iguessmd_sim.simulation.integrator.getTemperature()._value
                == general_iguessmd_data["temperature_K"]
        )
        assert (
                iguessmd_sim.simulation.integrator.getStepSize()._value
                == general_iguessmd_data["timestep_ps"]
        )


@pytest.mark.parametrize(
    "positions, masses, com",
    [TEST_COM_TWO_ATOMS, TEST_COM_METHANE, TEST_COM_CIRCLE, TEST_COM_CUBE],
)
def test_calculate_com_iguessmd_simulation_class(positions, masses, com):
    """
    Check that the function _calculate_com correctly calculates
    the centre of mass of the atoms to which the iGUESSMD force is
    applied in the OMMiGUESSMDSimulationCOM class, given their
    positions and masses.
    """
    # Create the simulation and retrieve indices for all atoms
    simulation = build_com_simulation((positions, masses, com))
    indices = np.array([i for i in range(masses.size)])

    # Create the iGUESSMD simulation
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        simulation,
        indices,
        TEST_iGUESSMD_LINEAR_PATH,
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
    )

    # Calculate the COM using the class function to caand check it against the expected COM
    calculated_com = iguessmd_sim._calculate_com(positions, masses)
    assert np.allclose(calculated_com, com, atol=1e-16)


@pytest.mark.parametrize(
    "positions, masses, com",
    [TEST_COM_TWO_ATOMS, TEST_COM_METHANE, TEST_COM_CIRCLE, TEST_COM_CUBE],
)
def test_calculate_com_trajectory_iguessmd_simulation_class(positions, masses, com):
    """
    Check that the function _calculate_com_trajectory correctly calculates
    the trajectory of the centre of mass of the atoms to which the iGUESSMD force is
    applied in the OMMiGUESSMDSimulationCOM class.
    """
    # Create the simulation and retrieve indices for all atoms
    simulation = build_com_simulation((positions, masses, com))
    indices = np.array([i for i in range(masses.size)])

    # Create the iGUESSMD simulation
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        simulation,
        indices,
        TEST_iGUESSMD_LINEAR_PATH,
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
    )

    # Manually set the trajectory of atom positions and calculate the
    # corresponding trajectory of the COM
    atom_positions = np.zeros((TEST_iGUESSMD_LINEAR_PATH.shape[0], *positions.shape))
    expected_com_array = np.zeros(TEST_iGUESSMD_LINEAR_PATH.shape)
    for i in range(TEST_iGUESSMD_LINEAR_PATH.shape[0]):
        atom_positions[i] = positions + np.array(
            [TEST_iGUESSMD_LINEAR_PATH[i] for j in range(indices.size)]
        )
        expected_com_array[i] = com + TEST_iGUESSMD_LINEAR_PATH[i]

    # Set iGUESSMD atom positions equal to the trajectory of calculated atom positions
    iguessmd_sim.iguessmd_simulation_atom_positions = atom_positions

    # Calculate COM trajectory using internal function and check the calculated
    # COMs match the predicted COMs
    iguessmd_sim._calculate_com_trajectory()
    assert np.allclose(iguessmd_sim.iguessmd_com_positions, expected_com_array, atol=1e-16)


@pytest.mark.parametrize("pbcs", TEST_BOOLS)
def test_iguessmd_com_force(pbcs):
    """
    Check that the force produced by the function iguessmd_com_force returns
    a force with the correct properties.
    """
    iguessmd_force = iguessmd_com_force(
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
        uses_pbcs=pbcs,
    )
    assert type(iguessmd_force) == CustomCentroidBondForce
    assert iguessmd_force.usesPeriodicBoundaryConditions() == pbcs
    assert (
        iguessmd_force.getEnergyFunction() == iGUESSMD_FORCE_EXPRESSION_COM_NONPERIODIC
        or iguessmd_force.getEnergyFunction() == iGUESSMD_FORCE_EXPRESSION_COM_PERIODIC
    )
    assert (
        iguessmd_force.getGlobalParameterName(0)
        == iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME
    )
    assert (
        iguessmd_force.getGlobalParameterName(1)
        == iGUESSMD_FORCE_CONSTANT_PERPENDICULAR_PARAMETER_NAME
    )
    assert (
        iguessmd_force.getGlobalParameterDefaultValue(0)
        == TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL
    )
    assert iguessmd_force.getNumPerBondParameters() == 6
    assert iguessmd_force.getPerBondParameterName(0) == "x0"
    assert iguessmd_force.getPerBondParameterName(1) == "y0"
    assert iguessmd_force.getPerBondParameterName(2) == "z0"
    assert iguessmd_force.getPerBondParameterName(3) == "tx"
    assert iguessmd_force.getPerBondParameterName(4) == "ty"
    assert iguessmd_force.getPerBondParameterName(5) == "tz"
    assert iguessmd_force.getForceGroup() == 31


@pytest.mark.parametrize("pbcs", TEST_BOOLS)
def test_iguessmd_single_atom_force(pbcs):
    """
    Check that the force produced by the function iguessmd_single_atom_force
    returns a force with the correct properties.
    """
    iguessmd_force = iguessmd_single_atom_force(
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
        uses_pbcs=pbcs,
    )
    assert type(iguessmd_force) == CustomExternalForce
    if pbcs:
        assert (
            iguessmd_force.getEnergyFunction()
            == iGUESSMD_FORCE_EXPRESSION_ATOM_PERIODIC
        )
    else:
        assert (
            iguessmd_force.getEnergyFunction()
            == iGUESSMD_FORCE_EXPRESSION_ATOM_NONPERIODIC
        )
    assert (
        iguessmd_force.getGlobalParameterName(0)
        == iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME
    )
    assert (
        iguessmd_force.getGlobalParameterName(1)
        == iGUESSMD_FORCE_CONSTANT_PERPENDICULAR_PARAMETER_NAME
    )
    assert (
        iguessmd_force.getGlobalParameterDefaultValue(0)
        == TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL
    )
    assert iguessmd_force.getNumPerParticleParameters() == 6
    assert iguessmd_force.getPerParticleParameterName(0) == "x0"
    assert iguessmd_force.getPerParticleParameterName(1) == "y0"
    assert iguessmd_force.getPerParticleParameterName(2) == "z0"
    assert iguessmd_force.getPerParticleParameterName(3) == "tx"
    assert iguessmd_force.getPerParticleParameterName(4) == "ty"
    assert iguessmd_force.getPerParticleParameterName(5) == "tz"
    assert iguessmd_force.getForceGroup() == 31


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
    3 point path defining the positions of the iGUESSMD force, starting at the
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

    Thus, the work done by the iGUESSMD force applied in this simulation should
    be 0.01 kJ mol-1.
    """
    assert TEST_iGUESSMD_ARGON_PATH.shape == (3, 3)
    assert TEST_iGUESSMD_ARGON_FORCE_CONSTANT == 100.0
    # Create the iGUESSMD simulation
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_single_atom_simulation(),
        TEST_iGUESSMD_ARGON_INDEX,
        TEST_iGUESSMD_ARGON_PATH,
        fc_multiplier * TEST_iGUESSMD_ARGON_FORCE_CONSTANT,
    )

    # Run iGUESSMD procedure
    with redirect_stdout(StringIO()) as _:
        iguessmd_sim.run_iguessmd()

    # Check values of work done
    assert iguessmd_sim.iguessmd_simulation_work_done[-1] == fc_multiplier * 0.01
    assert np.array_equal(
        iguessmd_sim.iguessmd_simulation_work_done,
        fc_multiplier * np.array([0.0, 0.0, 0.01]),
    )


@pytest.mark.parametrize(
    "apply_pbcs, save_iguessmd_force, indices",
    product(TEST_BOOLS, TEST_BOOLS, TEST_iGUESSMD_INDICES),
)
def test_load_openmm_state(apply_pbcs, save_iguessmd_force, indices):
    """
    Check that the OMMiGUESSMDSimulation correctly loads the state
    of the system by checking that the velocities loaded are
    correct. This also implicitly tests the save_simulation
    function.
    """
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(apply_pbcs),
        indices,
        TEST_iGUESSMD_LINEAR_PATH,
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
    )

    # Run simulation for a few steps
    iguessmd_sim.run_equilibration_with_initial_restraint(n_steps=1000)
    original_velocities = iguessmd_sim.simulation.context.getState(
        getVelocities=True
    ).getVelocities(asNumpy=True)

    with redirect_stdout(StringIO()) as _:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Save simulation to file
            output_path = Path(tmpdir)
            filename = "test_velocities.xml"
            file_path = output_path.joinpath(filename)
            iguessmd_sim.save_simulation(
                output_filepath=file_path,
                save_state=True,
                save_iguessmd_force=save_iguessmd_force,
            )
            assert file_path.exists()

            # Load saved simulation
            loaded_iguessmd_sim = OMMiGUESSMDSimulation.from_xml_path(
                file_path,
                indices,
                TEST_iGUESSMD_LINEAR_PATH,
                TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
            )
            # Retrieve and compare velocities
            loaded_velocities = loaded_iguessmd_sim.simulation.context.getState(
                getVelocities=True
            ).getVelocities(asNumpy=True)
            assert np.allclose(original_velocities, loaded_velocities, rtol=1e-7)

@pytest.mark.parametrize("indices", TEST_iGUESSMD_INDICES)
def test_get_iguessmd_atom_positions_periodic(indices):
    """
    Test that the atom positions retrieved from a periodic simulation
    that has crossed the boundary are unwrapped correctly.
    """
    # Define test path that crosses PBC in negative x direction
    test_pbc_crossing_path = np.array([np.linspace(1, -10, 1000), np.zeros(1000), np.zeros(1000)]).transpose()

    # Define iGUESSMD simulation
    iguessmd_sim = OMMiGUESSMDSimulation.from_simulation(
        build_basic_simulation(pbcs=True),
        indices,
        test_pbc_crossing_path,
        TEST_iGUESSMD_FORCE_CONSTANT_SPHERICAL,
    )

    # Run iGUESSMD to make molecule cross PBC in x direction
    iguessmd_sim.run_iguessmd()

    # Retrieve wrapped and unwrapped positions
    wrapped_positions = iguessmd_sim.simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(asNumpy=True)[indices]._value
    unwrapped_positions = \
    iguessmd_sim.simulation.context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(asNumpy=True)[
        indices]._value

    # Check that x values of atoms to which iGUESSMD force have been
    # unwrapped correctly
    if isinstance(iguessmd_sim, OMMiGUESSMDSimulationAtom):
        assert (wrapped_positions[0] != unwrapped_positions[0]).all()
        assert np.all(iguessmd_sim.iguessmd_simulation_atom_positions[-1,0] < 0.0)
        assert (unwrapped_positions[0] == iguessmd_sim.iguessmd_simulation_atom_positions[-1,0])
    elif isinstance(iguessmd_sim, OMMiGUESSMDSimulationCOM):
        assert (wrapped_positions[:, 0] != unwrapped_positions[:, 0]).all()
        assert np.all(iguessmd_sim.iguessmd_simulation_atom_positions[-1,:,0] < 0.0)
        assert (unwrapped_positions[:, 0] == iguessmd_sim.iguessmd_simulation_atom_positions[-1,:,0]).all()

