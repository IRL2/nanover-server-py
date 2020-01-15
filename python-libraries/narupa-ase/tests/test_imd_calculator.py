# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import time
from typing import Dict

import pytest
import numpy as np
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.md import Langevin
from ase.md.nvtberendsen import NVTBerendsen
from narupa.imd import ImdClient

from narupa.ase import converter
from narupa.core.timing import delayed_generator

from narupa.imd.imd_server import ImdServer
from narupa.ase.imd_calculator import ImdCalculator, get_periodic_box_lengths, _get_cancelled_interactions, \
    _get_atoms_to_reset
from narupa.imd.particle_interaction import ParticleInteraction
from util import co_atoms

@pytest.fixture
def imd_client():
    with ImdClient.insecure_channel(
            address='localhost', port=54322) as client:
        yield client



@pytest.fixture
def atoms():
    return co_atoms()


@pytest.fixture
def interact_c():
    interaction = ParticleInteraction(position=[1, 0, 0], particles=[0], scale=100., interaction_type='spring')
    return interaction


@pytest.fixture
def imd_calculator_co():
    server = ImdServer(address=None, port=None)
    atoms = co_atoms()
    calculator = LennardJones()
    imd_calculator = ImdCalculator(server.service, calculator, atoms)
    yield imd_calculator, atoms
    server.close()


@pytest.fixture
def imd_calculator_no_atoms():
    server = ImdServer(address=None, port=None)
    calculator = LennardJones()
    imd_calculator = ImdCalculator(server.service, calculator)
    yield imd_calculator
    server.close()

def test_imd_calculator_no_interactions(imd_calculator_co):
    imd_calculator, atoms = imd_calculator_co
    properties = ('energy', 'forces')
    imd_calculator.calculator.calculate(atoms=atoms, properties=properties)
    expected_results = imd_calculator.calculator.results
    imd_calculator.calculate()
    results = imd_calculator.results
    for key in properties:
        assert np.allclose(results[key], expected_results[key])
    assert results['interactive_energy'] == 0
    assert np.all(results['interactive_forces'] == np.zeros((len(atoms), 3)))


def test_imd_calculator_one_dimension_pbc(imd_calculator_co):
    imd_calculator, atoms = imd_calculator_co
    atoms.set_pbc((True, False, False))
    with pytest.raises(NotImplementedError):
        imd_calculator.calculate()


def test_imd_calculator_no_pbc(imd_calculator_co):
    imd_calculator, atoms = imd_calculator_co
    atoms.set_pbc((False, False, False))
    assert get_periodic_box_lengths(atoms) is None


def test_imd_calculator_not_orthorhombic(imd_calculator_co):
    imd_calculator, atoms = imd_calculator_co
    atoms.set_cell([1, 1, 1, 45, 45, 45])
    with pytest.raises(NotImplementedError):
        imd_calculator.calculate()


def test_imd_calculator_late_atoms(imd_calculator_no_atoms, atoms):
    """
    tests that the imd calculator works if atoms supplied after initialisation.
    """
    imd_calculator_no_atoms.calculate(atoms=atoms)


def test_imd_calculator_no_atoms(imd_calculator_no_atoms):
    """
    tests that the imd calculator throws an exception if no atoms are supplied.
    """
    with pytest.raises(ValueError):
        imd_calculator_no_atoms.calculate()


@pytest.mark.parametrize("position, imd_energy, imd_forces",
                         [([1, 0, 0], 1, [2, 0, 0, 0, 0, 0]),
                          ([3, 0, 0], 1, [-2, 0, 0, 0, 0, 0]),
                          ([-1, 0, 0], 1, [-2, 0, 0, 0, 0, 0]),
                          ([0, 1, 0], 1, [0, 2, 0, 0, 0, 0]),
                          ([0, 3, 0], 1, [0, -2, 0, 0, 0, 0]),
                          ([0, -1, 0], 1, [0, -2, 0, 0, 0, 0]),
                          ([0, 0, 1], 1, [0, 0, 2, 0, 0, 0]),
                          ([0, 0, 3], 1, [0, 0, -2, 0, 0, 0]),
                          ([0, 0, -1], 1, [0, 0, -2, 0, 0, 0]),
                          ([5, 0, 0], 1, [2, 0, 0, 0, 0, 0]),
                          ])
def test_one_interaction(position, imd_energy, imd_forces, imd_calculator_co, interact_c, imd_client):
    """
    tests an interaction in several different positions, including periodic boundary positions.
    """
    imd_calculator, atoms = imd_calculator_co
    # reshape the expected imd forces.
    imd_forces = np.array(imd_forces)
    imd_forces = imd_forces.reshape((2, 3))

    # perform the internal energy calculation.
    properties = ('energy', 'forces', 'interactive_energy', 'interactive_forces')
    imd_calculator.calculator.calculate(atoms=atoms, properties=properties)
    internal_energy = imd_calculator.calculator.results['energy']
    internal_forces = imd_calculator.calculator.results['forces']

    # perform the calculation with interaction applied.
    interact_c.position = position
    imd_client.publish_interactions_async(delayed_generator([interact_c] * 20, delay=0.01))
    time.sleep(0.05)
    assert len(imd_calculator.interactions) == 1
    imd_calculator.calculate(properties=properties)

    # set up the expected energy and forces.
    expected_imd_energy_kjmol = interact_c.scale * imd_energy * atoms.get_masses()[0]
    expected_imd_energy = expected_imd_energy_kjmol * converter.KJMOL_TO_EV
    expected_imd_forces = interact_c.scale * atoms.get_masses()[0] * (
            imd_forces * converter.KJMOL_TO_EV / converter.NM_TO_ANG)
    expected_forces = internal_forces + expected_imd_forces
    expected_energy = internal_energy + expected_imd_energy

    results = imd_calculator.results
    forces = results['forces']
    assert np.allclose(results['interactive_energy'], expected_imd_energy)
    assert np.allclose(results['interactive_forces'], expected_imd_forces)
    assert np.allclose(results['energy'], expected_energy)
    assert np.allclose(forces, expected_forces)


