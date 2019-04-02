"""
Tests for :mod:`narupa.openmm.runner`.
"""
import pytest

import numpy as np

import simtk.openmm as mm
from simtk.openmm import app
from simtk.unit import kelvin, picosecond, femtosecond, nanometer

from narupa.openmm import Runner


class DoNothingReporter:
    """
    OpenMM reporter that does nothing.

    The reporter does nothing but is valid. It is meant to populate the list of
    reporters of an OpenMM simulation.
    """
    # The name of the method is part of the OpenMM API. It cannot be made to
    # conform PEP8.
    def describeNextReport(self, simulation):  # pylint: disable=invalid-name
        """
        Activate the reporting every step, but collect no data.
        """
        return 1, False, False, False, False

    def report(self, simulation, state):
        """
        Do not report anything.
        """
        pass


@pytest.fixture
def basic_simulation():
    """
    Setup a minimal OpenMM simulation with two atoms.
    """
    periodic_box_vector = [
        [10,  0,  0],
        [ 0, 10,  0],
        [ 0,  0, 10],
    ]
    positions = np.array([[0, 0, 0], [0, 3, 0]], dtype=float)

    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue(name='RES', chain=chain)
    topology.addAtom(element=0, name='A1', residue=residue)
    topology.addAtom(element=0, name='A2', residue=residue)

    system = mm.System()
    system.setDefaultPeriodicBoxVectors(*periodic_box_vector)
    system.addParticle(mass=72)
    system.addParticle(mass=72)

    force = mm.NonbondedForce()
    force.setNonbondedMethod(force.NoCutoff)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    system.addForce(force)

    integrator = mm.LangevinIntegrator(300 * kelvin, 1 / picosecond, 2 * femtosecond)

    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPeriodicBoxVectors(*periodic_box_vector)
    simulation.context.setPositions(positions * nanometer)

    return simulation


@pytest.fixture
def runner(basic_simulation):
    """
    Setup a :class:`Runner` on a basic simulation.

    The simulation has a reporter attached to it to assure removing a reporter
    only removes only that reporter.
    """
    runner = Runner(basic_simulation)
    runner.simulation.reporters.append(DoNothingReporter())
    return runner


@pytest.fixture
def serialized_system(basic_simulation, tmp_path):
    """
    Setup an XML serialized system and a PDB as a temporary files.
    """
    xml_path = tmp_path / "system.xml"
    system_xml = mm.XmlSerializer.serialize(basic_simulation.system)
    with open(str(xml_path), 'w') as outfile:
        outfile.write(system_xml)

    pdb_path = tmp_path / 'system.pdb'
    positions = basic_simulation.context.getState(getPositions=True).getPositions()
    with open(str(pdb_path), 'w') as outfile:
        app.PDBFile.writeFile(basic_simulation.topology, positions, outfile)
    return xml_path, pdb_path


def test_default_verbosity(runner):
    """
    Test that the verbosity is off by default
    """
    assert not runner.verbose


@pytest.mark.parametrize('initial_value', (True, False))
@pytest.mark.parametrize('set_value_to', (True, False))
def test_set_verbosity_from_property(runner, initial_value, set_value_to):
    """
    Test that the verbosity can be set from the :attr:`Runner.verbose` property.

    The test makes sure that the value can be set from one value to an other,
    and from one value to itself.
    """
    expected_number_of_reporters = {
        True: 2,
        False: 1,
    }
    reporters = runner.simulation.reporters
    runner.verbose = initial_value
    assert runner.verbose == initial_value
    assert len(reporters) == expected_number_of_reporters[initial_value]
    runner.verbose = set_value_to
    assert runner.verbose == set_value_to
    assert len(reporters) == expected_number_of_reporters[set_value_to]


@pytest.mark.parametrize('initial_value', (True, False))
def test_make_verbose(runner, initial_value):
    """
    Test that :meth:`Runner.make_verbose` sets the verbosity on.
    """
    reporters = runner.simulation.reporters
    runner.verbose = initial_value
    assert runner.verbose == initial_value
    runner.make_verbose()
    assert runner.verbose == True
    assert len(reporters) == 2


@pytest.mark.parametrize('initial_value', (True, False))
def test_make_quiet(runner, initial_value):
    """
    Test that :meth:`Runner.make_quiet` sets the verbosity off.
    """
    reporters = runner.simulation.reporters
    runner.verbose = initial_value
    assert runner.verbose == initial_value
    runner.make_quiet()
    assert runner.verbose == False
    assert len(reporters) == 1


def test_run(runner):
    """
    Test that :meth:`Runner.run` runs the simulation.
    """
    assert runner.simulation.currentStep == 0
    runner.run(n_steps=5)
    assert runner.simulation.currentStep == 5


def test_from_xml_input(serialized_system):
    """
    Test that a :class:`Runner` can be built from a serialized system.
    """
    xml_path, pdb_path = serialized_system
    runner = Runner.from_xml_input(xml_path, pdb_path)
    assert runner.simulation.system.getNumParticles() == 2

