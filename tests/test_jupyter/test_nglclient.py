import numpy as np
import pytest
from ipywidgets import HTML

from nanover.jupyter import nglclient
from nanover.trajectory import FrameData


class _DummyView:
    def __init__(self):
        self.structures = []
        self.removed = []
        self.coordinates = []

    def add_structure(self, structure):
        self.structures.append(structure)
        return f"component-{len(self.structures)}"

    def remove_component(self, component):
        self.removed.append(component)

    def set_coordinates(self, coordinates):
        self.coordinates.append(coordinates)


def _minimum_usable_frame():
    frame = FrameData()
    frame.particle_count = 1
    frame.particle_elements = np.array([6], dtype=np.uint8)
    frame.bond_pairs = np.zeros((0, 2), dtype=np.uint32)
    return frame


def test_refresh_view_updates_coordinates(monkeypatch):
    monkeypatch.setattr(nglclient, "nglview", object())
    monkeypatch.setattr(nglclient, "FrameDataStructure", lambda frame: ("structure", frame))

    frame = _minimum_usable_frame()
    frame.particle_positions = np.array([[0.1, 0.2, 0.3]], dtype=np.float32)
    client = object.__new__(nglclient.NGLClient)
    client._view = _DummyView()
    client._structure = None
    client._current_frame = frame

    assert client.refresh_view() is True
    assert client._structure == "component-1"
    np.testing.assert_allclose(
        client._view.coordinates[0][0], frame.particle_positions * 10, rtol=0, atol=1e-7
    )


def test_refresh_view_resets_structure(monkeypatch):
    monkeypatch.setattr(nglclient, "nglview", object())
    monkeypatch.setattr(nglclient, "FrameDataStructure", lambda frame: ("structure", frame))

    frame = _minimum_usable_frame()
    frame.particle_positions = np.array([[1.0, 2.0, 3.0]], dtype=np.float32)
    client = object.__new__(nglclient.NGLClient)
    client._view = _DummyView()
    client._structure = "component-0"
    client._current_frame = frame

    assert client.refresh_view(reset_structure=True) is True
    assert client._view.removed == ["component-0"]
    assert client._structure == "component-1"


def test_refresh_view_handles_missing_positions(monkeypatch):
    monkeypatch.setattr(nglclient, "nglview", object())
    monkeypatch.setattr(nglclient, "FrameDataStructure", lambda frame: ("structure", frame))

    frame = _minimum_usable_frame()
    client = object.__new__(nglclient.NGLClient)
    client._view = _DummyView()
    client._structure = None
    client._current_frame = frame

    assert client.refresh_view() is False
    assert client._view.coordinates == []


def test_frame_data_to_nglwidget_without_nglview(monkeypatch):
    monkeypatch.setattr(nglclient, "nglview", None)

    widget = nglclient.frame_data_to_nglwidget(FrameData())

    assert isinstance(widget, HTML)
    assert "NGLView is not installed" in widget.value
    assert nglclient.is_nglview_available() is False


def test_nglclient_requires_nglview(monkeypatch):
    monkeypatch.setattr(nglclient, "nglview", None)

    with pytest.raises(ModuleNotFoundError, match="NGLView is required for NGLClient"):
        nglclient.NGLClient()
