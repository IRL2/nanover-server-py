# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Connects to a Narupa simulation and outputs osc messages for multiplayer number
fields and some metrics from framedata.
Run with:

.. code bash
    narupa-omm-ase nanotube.xml
    python generic.py --osc-port 9000
"""
from narupa.trajectory.frame_data import KINETIC_ENERGY, POTENTIAL_ENERGY
from osc_app import OscApp


def build_frame_generator(osc_client):
    def frame_to_osc_messages(frame):
        for key, value in osc_client.narupa_client.latest_multiplayer_values.items():
            try:
                yield f"/multiplayer/{key}", value.number_value
            except AttributeError:
                pass

        if KINETIC_ENERGY in frame.values:
            yield "/energy/kinetic", frame.kinetic_energy
        if POTENTIAL_ENERGY in frame.values:
            yield "/energy/potential", frame.potential_energy

    return frame_to_osc_messages


if __name__ == '__main__':
    app = OscApp(build_frame_generator)
    app.run()
