"""
Tests for :mod:`narupa.openmm.runner`.
"""
import pytest

import simtk.openmm as mm
from simtk.openmm import app

from narupa.openmm import Runner

from simulation_utils import (
    DoNothingReporter,
    basic_simulation,
)



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

