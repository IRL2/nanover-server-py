# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

import functools
import pytest
import numpy as np
from ase import Atoms, units
from ase.cell import Cell
from ase.md import VelocityVerlet
from narupa.ase.null_calculator import NullCalculator
from narupa.ase.openmm.calculator import OpenMMCalculator
from openmm_ase.simulation_utils import basic_simulation
from narupa.ase.wall_constraint import VelocityWallConstraint

@pytest.fixture
def walled_dynamics_and_expectations():
    """
    A MD system with a single atom moving toward the wall and the expectations.

    The expectations are the expected positions and velocity at each step.
    """
    position = [1, 2, 3]
    velocity = [-0.5, 1, 3.1]
    box = [2., 3.5, 6.]
    n_steps = 5

    atoms = Atoms(
        'C',
        positions=[[p * units.Ang for p in position]],
        pbc=(True, True, True),
    )
    atoms.set_masses([1])
    atoms.set_velocities([[v * units.Ang / units.fs for v in velocity]])

    cell_array = np.zeros((3, 3))
    np.fill_diagonal(cell_array, box)
    cell = Cell(cell_array)
    atoms.set_cell(cell)

    atoms.calc = NullCalculator()
    atoms.constraints.append(VelocityWallConstraint())
    dynamics = VelocityVerlet(atoms=atoms, timestep=1 * units.fs)

    expected_positions = np.zeros((n_steps, 3))
    expected_velocities = np.zeros((n_steps, 3))
    expected_positions[0, :] = position
    expected_velocities[0, :] = velocity
    for step in range(1, n_steps):
        step_position = expected_positions[step - 1, :].copy()
        step_velocity = expected_velocities[step - 1, :].copy()
        step_position += step_velocity
        inversion = (
            ((step_position <= 0) & (step_velocity < 0))
            | ((step_position >= box) & (step_velocity > 0))
        )
        step_velocity[inversion] *= -1
        expected_positions[step, :] = step_position
        expected_velocities[step, :] = step_velocity

    return dynamics, expected_positions, expected_velocities


@pytest.fixture
def walled_atoms(basic_simulation):
    calculator = OpenMMCalculator(basic_simulation)
    walled_atoms = calculator.generate_atoms()
    walled_atoms.calc = calculator
    walled_atoms.constraints.append(VelocityWallConstraint())
    return walled_atoms


@pytest.fixture
def openmm_calculator_and_atoms(basic_simulation):
    reference_calculator = OpenMMCalculator(basic_simulation)
    reference_atoms = reference_calculator.generate_atoms()
    reference_atoms.calc = reference_calculator
    return reference_calculator, reference_atoms


def test_velocity_wall(walled_dynamics_and_expectations):
    dynamics, expected_positions, expected_velocities = walled_dynamics_and_expectations
    n_steps = expected_positions.shape[0] - 1
    positions = []
    velocities = []

    def register_coordinates(atoms):
        """
        During the dynamics, store the positions and velocities.
        Positions are expressed in Å and velocities in Å/fs.
        """
        positions.append(atoms.get_positions()[0])
        velocities.append(atoms.get_velocities()[0] * (units.fs / units.Ang))

    dynamics.attach(functools.partial(register_coordinates, dynamics.atoms), interval=1)
    dynamics.run(n_steps)

    assert np.allclose(positions, expected_positions)
    assert np.allclose(velocities, expected_velocities)


@pytest.mark.parametrize('broken_cell', (
    Cell(np.zeros((3, 3), dtype=float)),  # Nul volume
    Cell(np.array([[1., 1., 0.], [0., 1., 1.], [0., 0., 1.]])), # Not orthorhombic
))
def test_validate_box(broken_cell):
    with pytest.raises(ValueError):
        VelocityWallConstraint._validate_box(broken_cell)


def test_chaining_calculators(walled_atoms, openmm_calculator_and_atoms):
    reference_calculator, reference_atoms = openmm_calculator_and_atoms
    velocities = np.ones((len(reference_atoms), 3)) * units.Ang / units.fs
    reference_atoms.set_velocities(velocities)
    reference_dynamics = VelocityVerlet(atoms=reference_atoms, timestep=1 * units.fs)

    walled_atoms = walled_atoms
    walled_atoms.set_velocities(velocities)
    walled_dynamics = VelocityVerlet(atoms=walled_atoms, timestep=1 * units.fs)

    # We make sure we do not use the same atom object for the reference and for
    # the walled system. This should not happen expect in the event of
    # accidentally referencing the same object during construction of the test.
    assert walled_atoms is not reference_atoms

    # The first 7 steps should be identical with and without the wall:
    # * some of the coordinates are negative, but the velocity is pushing toward
    # the box so it is not modified by the wall;
    # * one hydrogen has an initial Z of 14.359 Å with a box size of 20 Å and a
    # velocity of 1 Å/fs.
    reference_dynamics.run(7)
    walled_dynamics.run(7)
    assert walled_atoms.get_potential_energy() == reference_atoms.get_potential_energy()
    assert walled_atoms.get_positions() == pytest.approx(reference_atoms.get_positions())

    # During the 7th step, the further hydrogen should have crossed the box. So,
    # at the 8th step, the wall should put it back in only in the walled system.
    reference_dynamics.run(1)
    walled_dynamics.run(1)
    # Only the positions are different. The wall changed the velocity during the
    # integration step, but the energy will only be affected next step.
    assert walled_atoms.get_potential_energy() == reference_atoms.get_potential_energy()
    assert walled_atoms.get_positions() != pytest.approx(reference_atoms.get_positions())
