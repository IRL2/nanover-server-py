"""
Module providing an nglview client to a narupa frame server.
"""
from narupa.app import NarupaClient


class NGLClient(NarupaClient):
    """
    nglview client to a narupa server.
    """
    def __init__(self, *args, update_callback=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._view = None
        self.update_callback = update_callback

    @property
    def view(self):
        if self._view is None:
            atoms = frame_data_to_ase(
                client.first_frame,
                topology=True,
                positions=False,
            )
            atoms.set_positions(
                np.array(client.latest_frame.particle_positions) * 10
            )
            self._view = nglview.show_ase(atoms)
        return self._view

    def _on_frame_received(self, frame_index: int, frame):
        super()._on_frame_received(frame_index, frame)
        self.view.set_coordinates(
            {0: np.array(self.latest_frame.particle_positions) * 10}
        )
        if self.update_callback is not None:
            self.update_callback(self.universe)