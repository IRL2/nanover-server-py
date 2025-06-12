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
  when run with the SMD force not located at the initial position []
- The SMD simulation can be saved correctly, with or without the SMD force []
- If the SMD simulation is being loaded from a NanoVer OpenMM XML file created
  via one of the SMDSimulation classes and contains an SMD force already, this
  SMD force is correctly loaded and matches the expected force constant specified
  when creating the class []
- The class can generate the correct number of starting structures in the specified
  time interval, and that these are saved to the correct location []
- Running an SMD simulation produces reasonable results for the cumulative work done
  (may need to think about a specific test case for this...) []
- _calculate_forces works as expected []
- _calculate_work_done works as expected []
- Simulation data is saved in the correct format to the correct location, and can be
  subsequently loaded back into python correctly []
- General SMD data is saved in the correct format to the correct location, and can be
  subsequently loaded back into python correctly []
- For OpenMMSMDSimulationCOM, the COM of the specified atoms is correctly calculated []
- For OpenMMSMDSimulationCOM, the trajectory of the COM of the specified atoms is
  correctly calculated []
- smd_com_force works as expected []
- smd_single_atom_force works as expected []
"""
import copy

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

from nanover.omni.openmm import OpenMMSimulation
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
TEST_SMD_PATH = np.array(
    [np.linspace(0.05, 1.05, 101), np.zeros(101), np.zeros(101)]
).transpose()
TEST_SMD_FORCE_CONSTANT = 3011.0


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


@pytest.fixture
def make_basic_simulation_xml(tmp_path):
    serialized_simulation = serializer.serialize_simulation(
        build_basic_simulation(), save_state=True
    )
    xml_path = tmp_path / "basic_simulation.xml"
    with open(str(xml_path), "w") as xml_file:
        xml_file.write(serialized_simulation)
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
    """
    smd_sim = OpenMMSMDSimulation.from_simulation(
        build_basic_simulation(),
        indices,
        TEST_SMD_PATH,
        TEST_SMD_FORCE_CONSTANT,
    )
    # TODO: May wish to consider suppressing the output of the run_smd()
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
    smd_sim_state = smd_sim.simulation.context.getState(getPositions=True, getVelocities=True, getForces=True,
                                                        getEnergy=True)
    smd_sim_copy_state = smd_sim_copy.simulation.context.getState(getPositions=True, getVelocities=True, getForces=True,
                                                                  getEnergy=True)
    assert np.array_equal(smd_sim_state.getPositions(asNumpy=True), smd_sim_copy_state.getPositions(asNumpy=True))
    assert np.array_equal(smd_sim_state.getVelocities(asNumpy=True), smd_sim_copy_state.getVelocities(asNumpy=True))
    assert np.array_equal(smd_sim_state.getForces(asNumpy=True), smd_sim_copy_state.getForces(asNumpy=True))
    assert smd_sim_state.getKineticEnergy() == smd_sim_copy_state.getKineticEnergy()
    assert smd_sim_state.getPotentialEnergy() == smd_sim_copy_state.getPotentialEnergy()
    assert np.array_equal(smd_sim.smd_simulation_atom_positions, smd_sim_copy.smd_simulation_atom_positions)

    # Check the arguments passed to the OpenMMSMDSimulation class are unchanged by the reset
    assert np.array_equal(smd_sim.smd_atom_indices, smd_sim_copy.smd_atom_indices)
    assert np.array_equal(smd_sim.smd_path, smd_sim_copy.smd_path)
    assert np.array_equal(smd_sim.smd_force_constant, smd_sim_copy.smd_force_constant)

    # Check that the SMD force attached to the simulation is correctly reset
    assert np.array_equal(smd_sim.current_smd_force_position, smd_sim_copy.current_smd_force_position)
    assert smd_sim.current_smd_force_position_index == smd_sim_copy.current_smd_force_position_index
    assert smd_sim.smd_force.getEnergyFunction() == smd_sim_copy.smd_force.getEnergyFunction()
    # Class-specific checks
    if type(smd_sim) == OpenMMSMDSimulationAtom:
        assert smd_sim.smd_force.getParticleParameters(0) == smd_sim_copy.smd_force.getParticleParameters(0)
        assert smd_sim.smd_force.getNumParticles() == 1
    elif type(smd_sim) == OpenMMSMDSimulationCOM:
        assert smd_sim.smd_force.getGroupParameters(0) == smd_sim_copy.smd_force.getGroupParameters(0)
        assert smd_sim.smd_force.getNumBonds() == 1

    # Check other relevant properties of the OpenMMSMDSimulation class match those of the
    # fresh copy after the reset
    assert np.array_equal(smd_sim.smd_simulation_atom_positions, smd_sim_copy.smd_simulation_atom_positions)


@pytest.mark.parametrize("indices", [TEST_SMD_SINGLE_INDEX, TEST_SMD_MULTIPLE_INDICES])
def test_smd_force_added_to_system(indices):
    """
    Check that the last force to be added to the OpenMM simulation is the SMD force added during
    initialisation of the OpenMMSMDSimulation class, with force group 31.
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


@pytest.mark.parametrize("indices", [TEST_SMD_SINGLE_INDEX, TEST_SMD_MULTIPLE_INDICES])
def test_smd_force_updates_correctly(indices):
    """
    Check that the position of the SMD force is correctly updated upon calling
    update_smd_force_position().
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
        sys_index, sys_position = smd_sim.simulation.system.getForce(n_system_forces-1).getParticleParameters(0)
        assert sys_index == indices
        assert np.array_equal(np.array(sys_position), new_force_position)

    elif type(smd_sim.smd_force) == CustomCentroidBondForce:

        # OpenMMSMDSimulationCOM force parameters
        _, bond_params = smd_sim.smd_force.getBondParameters(0)
        assert np.array_equal(np.array(bond_params), new_force_position)

        # Force parameters from system
        _, sys_bond_params = smd_sim.simulation.system.getForce(n_system_forces-1).getBondParameters(0)
        assert np.array_equal(np.array(sys_bond_params), new_force_position)

