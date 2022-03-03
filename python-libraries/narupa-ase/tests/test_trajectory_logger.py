import os
import time
import re
from typing import List

import pytest
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.md import Langevin

from narupa.ase.trajectory_logger import TrajectoryLogger, UnsupportedFormatError
from test_imd_reset import fcc_atoms
from ase.io import read, write, Trajectory
import numpy as np

SUPPORTED_EXTENSIONS = ['xyz']
NONEXISTANT_EXTENSION = 'babyyoda'
FRAMES = 10


@pytest.fixture()
def atoms():
    atoms = fcc_atoms()
    atoms.calc = LennardJones()
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

    with TrajectoryLogger(atoms, file) as logger:
        logger.write()
        assert os.path.exists(logger.current_path)
        path = logger.current_path
    check_file_images([atoms], path)


@pytest.mark.parametrize('ext', SUPPORTED_EXTENSIONS)
def test_write_multiple_frames(atoms, tmp_dir, ext):
    file = os.path.join(tmp_dir, "atoms" + "." + ext)

    with TrajectoryLogger(atoms, file) as logger:
        path = logger.current_path
        for i in range(FRAMES):
            logger.write()

    assert os.path.exists(path)
    check_file_images([atoms] * FRAMES, path)


def test_append_unsupported_filename(atoms, tmp_dir):
    file = os.path.join(tmp_dir, "atoms" + "." + NONEXISTANT_EXTENSION)

    with pytest.raises(UnsupportedFormatError):
        _ = TrajectoryLogger(atoms, file)


def test_attach_to_md(atoms, tmp_dir):
    file = os.path.join(tmp_dir, "atoms" + ".xyz")

    ase_traj_file = os.path.join(tmp_dir, "test.traj")
    with Trajectory(ase_traj_file, 'w', atoms) as ase_trajectory, TrajectoryLogger(atoms, file) as logger:
        path = logger.current_path
        md = Langevin(atoms, 0.5,
                      temperature_K=300,
                      friction=0.001,
                      trajectory=ase_trajectory)
        md.attach(logger)
        md.run(FRAMES)

    assert os.path.exists(ase_traj_file)
    assert os.path.exists(path)

    with Trajectory(ase_traj_file, 'r') as ase_trajectory_read:
        frames = [atoms for atoms in ase_trajectory_read]

    assert len(frames) == FRAMES + 1
    check_file_images(frames, path)


def test_timestamp(atoms, tmp_dir):
    """
    if using timestamp, the filename the logger uses should have a timestamp in it.
    """
    path_regex = re.compile(r".*\d{4}_\d{2}_\d{2}__\d{2}_\d{2}_\d{2}_\d{2}\.xyz")
    file = os.path.join(tmp_dir, "atoms" + ".xyz")
    with TrajectoryLogger(atoms, file, timestamp=True) as logger:
        assert logger.current_path != file
        assert path_regex.match(logger.current_path) is not None


def test_no_timestamp(atoms, tmp_dir):
    file = os.path.join(tmp_dir, "atoms" + ".xyz")
    with TrajectoryLogger(atoms, file, timestamp=False) as logger:
        assert logger.current_path == file


def test_reset_timestamp(atoms, tmp_dir):
    """
    tests that resetting the logger that's using timestamps sets the frame index back to zero,
    and creates a new file.
    """
    file = os.path.join(tmp_dir, "atoms" + ".xyz")
    with TrajectoryLogger(atoms, file) as logger:
        logger.write()
        initial_name = logger.current_path
        assert os.path.exists(initial_name)
        assert logger.frame_index > 0
        time.sleep(0.01)  # sleep briefly to ensure that timestamp will be different
        logger.reset()
        assert logger.current_path != initial_name
        assert logger.frame_index == 0


def test_reset_no_timestamp(atoms, tmp_dir):
    """
    tests that a logger not using a timestamp will overwrite the old file
    when resetting.
    """
    file_path = os.path.join(tmp_dir, "atoms" + ".xyz")
    with TrajectoryLogger(atoms, file_path, timestamp=False) as logger:
        logger.write()
        initial_name = logger.current_path
        assert os.path.exists(initial_name)
        assert logger.frame_index > 0
        logger.reset()
        assert logger.current_path == initial_name
        assert logger.frame_index == 0
        logger.write()
        path = logger.current_path
    # check only one frame has been written
    check_file_images([atoms], path)
