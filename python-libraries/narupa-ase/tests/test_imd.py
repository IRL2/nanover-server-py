# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import time
from typing import Iterable, Optional, Union

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes, all_properties
from ase.cell import Cell
from ase.md import VelocityVerlet
from narupa.ase.imd import NarupaASEDynamics
from narupa.core import NarupaClient
from narupa.trajectory.frame_server import PLAY_COMMAND_KEY, PAUSE_COMMAND_KEY, RESET_COMMAND_KEY, STEP_COMMAND_KEY

DUMMY_ATOMS_COUNT = 4
DUMMY_ATOMS_POSITIONS = np.arange(DUMMY_ATOMS_COUNT * 3, dtype=float).reshape((-1, 3))
DUMMY_ATOMS_VELOCITIES = DUMMY_ATOMS_POSITIONS[::-1]
DUMMY_ATOMS_CELL = Cell(np.array([
    [10., 0., 0.],
    [0., 11., 0.],
    [0., 0., 12.],
], dtype=float))
OneOrTwoDArray = Union[Iterable[float], Iterable[Iterable[float]], np.ndarray]


class ArbitraryCalculator(Calculator):
    """
    A calculator that sets the positions, velocities, forces, box, and potential
    energy to arbitrary values.

    The values to set are stored in the attribute named ``request_*`. The
    requested positions, velocities, and forces can be expressed as values for
    all the atoms, or as the value for one atom that will be repeated.
    """
    # We set the defaults as class attributes. If they are overwritten for an
    # instance, then the instance will use its own version instead of the
    # class one.
    request_positions: OneOrTwoDArray = np.array([0, 0, 0])
    request_velocities: OneOrTwoDArray = np.array([0, 0, 0])
    request_forces: OneOrTwoDArray = np.array([0, 0, 0])
    request_box: Cell = Cell(
        np.array([
            [1, 0, 0],
            [0, 2, 0],
            [0, 0, 3],
        ], dtype=float)
    )
    request_potential_energy: float = 0

    def calculate(self, atoms: Optional[Atoms] = None, properties=all_properties,
                  system_changes=all_changes):
        if Atoms is None:
            raise ValueError('No atom provided to the calculator.')

        per_atom_suffixes = ('positions', 'velocities')
        for suffix in per_atom_suffixes:
            attribute_name = f'request_{suffix}'
            attribute = getattr(self, attribute_name)
            attribute = self._repeat_per_atom(atoms, attribute)
            set_method_name = f'set_{suffix}'
            set_method = getattr(atoms, set_method_name)
            set_method(attribute)
        atoms.set_cell(self.request_box)

        forces = self._repeat_per_atom(atoms, self.request_forces)
        self.results['forces'] = forces
        self.results['energy'] = self.request_potential_energy

    @staticmethod
    def _repeat_per_atom(atoms, value: OneOrTwoDArray) -> np.ndarray:
        try:
            _ = value[0][0]
        except (TypeError, IndexError):
            # We get a `TypeError` if value is a list, a `IndexError` if it is
            # a numpy array.
            value = np.tile(value, len(atoms)).reshape((-1, 3))
        return value


def do_nothing_producer(*args, **kwrags):
    """
    Take arguments and return a function that takes and does nothing.

    Useful to provide a dummy callback.
    """
    return lambda: None


@pytest.fixture
def dummy_atoms():
    n_atoms = 4
    atoms = Atoms(positions=DUMMY_ATOMS_POSITIONS.copy())
    atoms.set_velocities(DUMMY_ATOMS_VELOCITIES.copy())
    atoms.set_cell(DUMMY_ATOMS_CELL.copy())
    return atoms


@pytest.fixture
def arbitrary_dynamics(dummy_atoms):
    dynamics = VelocityVerlet(dummy_atoms, timestep=1)
    calculator = ArbitraryCalculator()
    dynamics.atoms.set_calculator(calculator)
    return dynamics


@pytest.fixture
def arbitrary_ase_server(arbitrary_dynamics):
    with NarupaASEDynamics.basic_imd(arbitrary_dynamics,
                                     port=0,
                                     frame_method=do_nothing_producer) as ase_server:
        yield ase_server


@pytest.fixture
def client_server(arbitrary_ase_server):
    with NarupaClient.insecure_channel(address='localhost', port=arbitrary_ase_server.port) as c:
        yield c, arbitrary_ase_server


def test_address(arbitrary_ase_server):
    assert arbitrary_ase_server.address == '[::]'


def test_reset(arbitrary_ase_server):
    arbitrary_ase_server.run(1)
    positions = arbitrary_ase_server.atoms.get_positions()

    # Making sure that running the dynamics changed the system.
    assert not np.allclose(positions, DUMMY_ATOMS_POSITIONS)

    arbitrary_ase_server.reset()
    positions = arbitrary_ase_server.atoms.get_positions()
    velocities = arbitrary_ase_server.atoms.get_velocities()
    cell = arbitrary_ase_server.atoms.get_cell()
    assert np.allclose(positions, DUMMY_ATOMS_POSITIONS)
    assert np.allclose(velocities, DUMMY_ATOMS_VELOCITIES)
    assert np.allclose(cell, DUMMY_ATOMS_CELL)


def test_on_reset_event(arbitrary_ase_server):
    received = []

    def on_reset_0():
        received.append(0)

    def on_reset_1():
        received.append(1)

    arbitrary_ase_server.on_reset_listeners.append(on_reset_0)
    arbitrary_ase_server.on_reset_listeners.append(on_reset_1)
    arbitrary_ase_server.run(1)
    arbitrary_ase_server.reset()

    assert received == [0, 1]


def test_auto_reset(arbitrary_ase_server):
    reset = False

    def on_reset():
        nonlocal reset
        reset = True

    arbitrary_ase_server.on_reset_listeners.append(on_reset)

    arbitrary_ase_server.internal_calculator.requested_potential_energy = 1
    arbitrary_ase_server.run(1, reset_energy=10)
    assert not reset

    arbitrary_ase_server.internal_calculator.request_potential_energy = 1e9
    arbitrary_ase_server.run(1, reset_energy=10)
    assert reset


def test_run_blocking(arbitrary_ase_server):
    frame_count = 0

    def count_frames(*args, **kwargs):
        nonlocal frame_count
        frame_count += 1

    arbitrary_ase_server.dynamics.attach(count_frames, interval=1)
    arbitrary_ase_server.run(100)
    # ASE calls the observer callbacks once before the first step.
    assert frame_count == 101


def test_run_non_blocking(arbitrary_ase_server):
    frame_count = 0

    def count_frames(*args, **kwargs):
        nonlocal frame_count
        frame_count += 1

    arbitrary_ase_server.dynamics.attach(count_frames, interval=1)
    arbitrary_ase_server.run(100, block=False)
    # Here we count on the context switching to the assertions before finishing
    # the run.
    assert frame_count < 100
    assert arbitrary_ase_server._run_task is not None
    assert not arbitrary_ase_server._cancelled


def test_cancel_run(arbitrary_ase_server):
    arbitrary_ase_server.run(block=False)
    assert arbitrary_ase_server.is_running
    arbitrary_ase_server.cancel_run(wait=True)
    assert arbitrary_ase_server.is_running is False


def test_cancel_never_running(arbitrary_ase_server):
    # Cancelling a non running simulation should not raise an exception.
    arbitrary_ase_server.cancel_run()


def test_run_twice(arbitrary_ase_server):
    arbitrary_ase_server.run(block=False)
    with pytest.raises(RuntimeError):
        arbitrary_ase_server.run(block=False)


def test_step(arbitrary_ase_server):
    step_count = arbitrary_ase_server.dynamics.get_number_of_steps()
    arbitrary_ase_server.step()
    assert arbitrary_ase_server.dynamics.get_number_of_steps() == step_count + 1


def test_multiple_steps(arbitrary_ase_server):
    num_steps = 10
    arbitrary_ase_server.run(block=False)
    arbitrary_ase_server.cancel_run(wait=True)
    step_count = arbitrary_ase_server.dynamics.get_number_of_steps()
    for i in range(num_steps):
        arbitrary_ase_server.step()
    assert arbitrary_ase_server.dynamics.get_number_of_steps() == step_count + num_steps


def test_pause(arbitrary_ase_server):
    arbitrary_ase_server.run(block=False)
    arbitrary_ase_server.pause()
    step_count = arbitrary_ase_server.dynamics.get_number_of_steps()
    assert arbitrary_ase_server._run_task.done() is True
    time.sleep(0.1)
    assert arbitrary_ase_server.dynamics.get_number_of_steps() == step_count


def test_play(arbitrary_ase_server):
    arbitrary_ase_server.run(block=False)
    assert arbitrary_ase_server.is_running
    arbitrary_ase_server.pause()
    assert arbitrary_ase_server.is_running is False
    arbitrary_ase_server.play()
    assert arbitrary_ase_server.is_running


def test_play_twice(arbitrary_ase_server):
    arbitrary_ase_server.run(block=False)
    arbitrary_ase_server.pause()
    assert arbitrary_ase_server._run_task.done() is True
    arbitrary_ase_server.play()
    assert arbitrary_ase_server.is_running
    arbitrary_ase_server.play()
    assert arbitrary_ase_server.is_running


@pytest.mark.timeout(1)
def test_play_command(client_server):
    client, server = client_server
    assert not server.is_running
    client.run_command(PLAY_COMMAND_KEY)
    while True:
        if server.is_running:
            return


@pytest.mark.timeout(1)
def test_pause_command(client_server):
    client, server = client_server
    server.run()
    client.run_command(PAUSE_COMMAND_KEY)
    while server.is_running:
        continue
    step_count = server.dynamics.get_number_of_steps()
    time.sleep(0.1)
    assert server.dynamics.get_number_of_steps() == step_count


@pytest.mark.timeout(1)
def test_reset_command(client_server):
    client, server = client_server
    server.run()
    reset = False

    def on_reset():
        nonlocal reset
        reset = True

    server.on_reset_listeners.append(on_reset)
    client.run_command(RESET_COMMAND_KEY)

    while not reset:
        continue


@pytest.mark.timeout(1)
def test_step_command(client_server):
    client, server = client_server
    step_count = server.dynamics.get_number_of_steps()
    client.run_command(STEP_COMMAND_KEY)
    time.sleep(0.1)
    assert server.is_running is False
    assert server.dynamics.get_number_of_steps() == step_count + 1
