import os
import shutil
import tempfile
import time
from collections import Collection
from typing import NamedTuple, List

import pytest
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.md import Langevin

from narupa.ase.trajectory_logger import TrajectoryLogger
from test_imd_reset import fcc_atoms
from ase.io import read, write, Trajectory
import numpy as np

SUPPORTED_EXTENSIONS = ['xyz']
FRAMES = 10


@pytest.fixture()
def atoms():
    atoms = fcc_atoms()
    atoms.set_calculator(LennardJones())
    return atoms


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    path = tmpdir_factory.mktemp('outputs')
    return path


def check_file_images(images: List[Atoms], file, symbols=True, positions=True, cell=True):
    """

    Tests that a trajectory file is consistent with a reference trajectory of ASE atom objects

    :param images: Reference trajectory of ASE atom states.
    :param file: Trajectory file to check
    :param symbols: Whether to test chemical symbols
    :param positions: Whether to test positions
    :param cell: Whether to test unit cell
    """
    with open(file, 'r') as f:
        debug = f.read()

    new_images = read(file, index=':')
    assert len(new_images) == len(images)
    for atoms, new_atoms in zip(images, new_images):
        assert len(atoms) == len(new_atoms)
        if positions:
            assert np.allclose(atoms.positions, new_atoms.positions, atol=0.01)
        if symbols:
            assert atoms.get_chemical_symbols() == new_atoms.get_chemical_symbols()
        if cell:
            assert np.allclose(atoms.get_cell(), new_atoms.get_cell())


@pytest.mark.parametrize('ext', SUPPORTED_EXTENSIONS)
def test_write(atoms, tmp_dir, ext):
    file = os.path.join(tmp_dir, "atoms" + "." + ext)

    logger = TrajectoryLogger(atoms, file)
    logger.write()

    assert os.path.exists(file)
    check_file_images([atoms], file)


@pytest.mark.parametrize('ext', SUPPORTED_EXTENSIONS)
def test_write_multiple_frames(atoms, tmp_dir, ext):
    file = os.path.join(tmp_dir, "atoms" + "." + ext)

    logger = TrajectoryLogger(atoms, file)
    for i in range(FRAMES):
        logger.write()

    assert os.path.exists(file)
    check_file_images([atoms] * FRAMES, file)

@pytest.mark.parametrize('ext', SUPPORTED_EXTENSIONS)
def test_attach_to_md(atoms, tmp_dir, ext):
    file = os.path.join(tmp_dir, "atoms" + "." + ext)

    logger = TrajectoryLogger(atoms, file)

    ase_traj_file = os.path.join(tmp_dir, "test.traj")
    with Trajectory(ase_traj_file, 'w', atoms) as ase_trajectory:
        md = Langevin(atoms, 0.5, 300, 0.001, trajectory=ase_trajectory)
        md.attach(logger)
        md.run(FRAMES)

    assert os.path.exists(ase_traj_file)
    assert os.path.exists(file)

    with Trajectory(ase_traj_file, 'r') as ase_trajectory_read:
        frames = [atoms for atoms in ase_trajectory_read]

    assert len(frames) == FRAMES + 1
    check_file_images(frames, file)