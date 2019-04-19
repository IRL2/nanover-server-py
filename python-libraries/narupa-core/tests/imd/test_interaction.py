import numpy as np
import pytest

from narupa.imd.interaction import Interaction


@pytest.fixture
def interaction():
    return Interaction()


def test_player_id():
    interaction = Interaction("2")
    assert interaction.player_id == "2"


def test_get_default_position(interaction):
    assert np.allclose(interaction.position, [0, 0, 0])


def test_set_position(interaction):
    interaction.position = [1, 1, 1]
    assert np.allclose(interaction.position, [1, 1, 1])


def test_set_invalid_position(interaction):
    with pytest.raises(ValueError):
        interaction.position = [0, 0]


def test_get_default_atoms(interaction):
    assert len(interaction.atoms) == 0


def test_set_atoms(interaction):
    interaction.atoms = [0, 1, 2, 3, 4]
    assert np.allclose(interaction.atoms, [0, 1, 2, 3, 4])


def test_set_atom_unique(interaction):
    interaction.atoms = [0, 0, 0, 1, 2, 3, 4]
    assert np.allclose(interaction.atoms, [0, 1, 2, 3, 4])
