# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

import functools
import pytest
import numpy as np
from ase import Atoms, units
from ase.md import VelocityVerlet
from narupa.ase import wall_calculator


@pytest.fixture
def walled_dynamics_and_expectations():
    """
    A MD system with a single atom moving toward the wall and the expectations.

    The expectations are the expected positions and velocity at each step.
    """
    position = [1, 2, 3]
    velocity = [-0.5, 1, 3.1]
    box = [2, 3.5, 6]
    n_steps = 5

    atoms = Atoms(
        'C',
        positions=[[p * units.Ang for p in position]],
        pbc=box,
    )
    atoms.set_velocities([[v * 0.1 * units.Ang / units.fs for v in velocity]])
    atoms.set_masses([1])
    calculator = wall_calculator.VelocityWallCalculator(atoms=atoms)
    atoms.set_calculator(calculator)
    dynamics = VelocityVerlet(atoms=atoms, dt=1 * units.fs)

    expected_positions = np.zeros((n_steps, 3))
    expected_velocities = np.zeros((n_steps, 3))
    expected_positions[0, :] = position
    expected_velocities[0, :] = velocity
    for step in range(1, n_steps):
        step_position = expected_positions[step - 1, :].copy()
        step_velocity = expected_velocities[step - 1, :].copy()
        inversion = (
            ((step_position < 0) & (step_velocity < 0))
            | ((step_position > box) & (step_velocity > 0))
        )
        step_velocity[inversion] *= -1
        step_position += step_velocity
        expected_positions[step, :] = step_position
        expected_velocities[step, :] = step_velocity

    return dynamics, expected_positions, expected_velocities


def test_velocity_wall(walled_dynamics_and_expectations):
    dynamics, expected_positions, expected_velocities = walled_dynamics_and_expectations
    n_steps = expected_positions.shape[0]
    positions = []
    velocities = []

    def register_coordinates(atoms):
        positions.append(atoms.get_positions().copy())
        velocities.append(atoms.get_velocities().copy())

    dynamics.attach(functools.partial(register_coordinates, dynamics.atoms), interval=1)
    dynamics.run(expected_positions.shape[0])

    print(dynamics.atoms.get_masses())
    print(expected_positions)
    print('-' * 10)
    print(positions)
    print('-' * 10)
    print(expected_velocities)
    print('-' * 10)
    print(velocities)

    assert np.allclose(positions, expected_positions)
    assert np.allclose(velocities, expected_velocities)
