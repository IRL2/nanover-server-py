# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import pytest
from narupa.ase.converter import KJMOL_TO_EV
from narupa.ase.openmm import OpenMMCalculator
from simtk.unit import kilojoules_per_mole, nanometer, angstroms  # pylint: disable=no-name-in-module
import numpy as np

from .simulation_utils import basic_simulation, serialized_simulation_path


@pytest.fixture
def calculator(basic_simulation):
    return OpenMMCalculator(basic_simulation)


def test_energy(basic_simulation):
    calculator = OpenMMCalculator(basic_simulation)
    atoms = calculator.generate_atoms()
    expected_energy = basic_simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        kilojoules_per_mole)
    calculator.calculate(atoms, properties=['energy'])
    energy = calculator.results['energy'] / KJMOL_TO_EV
    assert expected_energy == pytest.approx(energy)


def test_forces(basic_simulation):
    calculator = OpenMMCalculator(basic_simulation)
    atoms = calculator.generate_atoms()
    expected_forces = basic_simulation.context.getState(getForces=True).getForces().value_in_unit(
        kilojoules_per_mole / nanometer)
    calculator.calculate(atoms, properties=['forces'])
    forces = calculator.results['forces'] / KJMOL_TO_EV * 10
    assert np.allclose(expected_forces, forces)


def test_from_serialized_sim(basic_simulation, serialized_simulation_path):
    calculator = OpenMMCalculator.from_xml(serialized_simulation_path)
    assert basic_simulation.system.getNumParticles() == calculator.simulation.system.getNumParticles()


def test_pbc(basic_simulation):
    calculator = OpenMMCalculator(basic_simulation)
    atoms = calculator.generate_atoms()
    assert np.allclose(atoms.get_cell(), [x.value_in_unit(angstroms) for x in
                                          calculator.simulation.system.getDefaultPeriodicBoxVectors()])


def test_no_atoms(basic_simulation):
    calculator = OpenMMCalculator(basic_simulation)
    with pytest.raises(ValueError):
        calculator.calculate()


def test_self_atoms(basic_simulation):
    calculator = OpenMMCalculator(basic_simulation)
    atoms = calculator.generate_atoms()
    calculator.atoms = atoms
    calculator.calculate(properties=['energy'])
    expected_energy = basic_simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        kilojoules_per_mole)
    energy = calculator.results['energy'] / KJMOL_TO_EV
    assert energy == pytest.approx(expected_energy)


def test_get_topology(basic_simulation):
    calculator = OpenMMCalculator(basic_simulation)
    assert calculator.topology == basic_simulation.topology
