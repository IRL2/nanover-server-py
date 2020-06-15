import pytest
from narupa.app.selection import *


def test_selection_defaults():
    selection = RenderingSelection("id", "name")
    assert selection.selection_name == "name"
    assert selection.selection_id == "id"
    assert selection.selected_particle_ids == set()
    assert selection.interaction_method == INTERACTION_METHOD_DEFAULT
    assert selection.renderer == RENDERER_DEFAULT
    assert selection.velocity_reset == VELOCITY_RESET_DEFAULT


def test_selection_dict_with_name_and_id():
    dict = {
        KEY_SELECTION_ID: "id",
        KEY_SELECTION_NAME: "name"
    }
    selection = RenderingSelection.from_dictionary(dict)
    assert selection.selection_name == "name"
    assert selection.selection_id == "id"
    assert selection.interaction_method == INTERACTION_METHOD_DEFAULT
    assert selection.renderer == RENDERER_DEFAULT
    assert selection.velocity_reset == VELOCITY_RESET_DEFAULT


def test_selection_dict_with_no_id():
    dict = {}
    with pytest.raises(ValueError):
        selection = RenderingSelection.from_dictionary(dict)


def test_selection_from_dictionary_with_interaction_method():
    dict = {
        KEY_SELECTION_ID: "id",
        KEY_SELECTION_NAME: "name",
        KEY_SELECTION_PROPERTIES: {
            KEY_PROPERTY_INTERACTION_METHOD: INTERACTION_GROUP
        }
    }
    selection = RenderingSelection.from_dictionary(dict)
    assert selection.interaction_method == INTERACTION_GROUP


def test_selection_from_dictionary_with_velocity_reset():
    dict = {
        KEY_SELECTION_ID: "id",
        KEY_SELECTION_NAME: "name",
        KEY_SELECTION_PROPERTIES: {
            KEY_PROPERTY_VELOCITY_RESET: True
        }
    }
    selection = RenderingSelection.from_dictionary(dict)
    assert selection.velocity_reset == True


def test_selection_from_dictionary_with_renderer():
    dict = {
        KEY_SELECTION_ID: "id",
        KEY_SELECTION_NAME: "name",
        KEY_SELECTION_PROPERTIES: {
            KEY_PROPERTY_RENDERER: "some_renderer"
        }
    }
    selection = RenderingSelection.from_dictionary(dict)
    assert selection.renderer == "some_renderer"


def test_set_selection():
    selection = RenderingSelection("id", "name")
    selection.set_particles({0, 1, 2, 3})
    assert selection.selected_particle_ids == {0, 1, 2, 3}
    selection.set_particles({2, 4, 6, 8})
    assert selection.selected_particle_ids == {2, 4, 6, 8}


def test_add_selection():
    selection = RenderingSelection("id", "name")
    selection.add_particles({0, 1, 2, 3})
    assert selection.selected_particle_ids == {0, 1, 2, 3}
    selection.add_particles({2, 4, 6, 8})
    assert selection.selected_particle_ids == {0, 1, 2, 3, 4, 6, 8}


def test_clear_selection():
    selection = RenderingSelection("id", "name")
    selection.add_particles({0, 1, 2, 3})
    assert selection.selected_particle_ids == {0, 1, 2, 3}
    selection.clear_particles()
    assert selection.selected_particle_ids == set()


def test_add_selection_empty():
    selection = RenderingSelection("id", "name")
    selection.add_particles({0, 1, 2, 3})
    assert selection.selected_particle_ids == {0, 1, 2, 3}
    selection.add_particles(())
    assert selection.selected_particle_ids == {0, 1, 2, 3}


def test_set_selection_empty():
    selection = RenderingSelection("id", "name")
    selection.add_particles({0, 1, 2, 3})
    assert selection.selected_particle_ids == {0, 1, 2, 3}
    selection.set_particles(())
    assert not selection.selected_particle_ids


def test_selection_updated():
    selection = RenderingSelection("id", "name")

    def callback(selection):
        callback.call_count += 1

    callback.call_count = 0

    selection.updated.add_callback(callback)

    assert callback.call_count == 0

    with selection.modify():
        selection.set_particles([0, 1, 2])

    assert callback.call_count == 1


def test_selection_removed():
    selection = RenderingSelection("id", "name")

    def callback(selection):
        callback.call_count += 1

    callback.call_count = 0

    selection.removed.add_callback(callback)

    assert callback.call_count == 0

    selection.remove()

    assert callback.call_count == 1