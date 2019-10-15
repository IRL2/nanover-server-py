# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Server to run an OpenMM simulation and publish its frames for Narupa.
"""

from typing import Optional

from simtk.openmm.app import Simulation

from narupa.trajectory.frame_server import FrameServer
from .runner import Runner
from .narupareporter import NarupaReporter
from .serializer import deserialize_simulation


class Server(Runner):
    """
    Run and serve an OpenMM simulation for Narupa.

    This server extends the :class:`Runner` class and adds the ability to
    publish the frames for Narupa.

    :param simulation: An instance of OpenMM :class:`Simulation` to run.
    :param address: The address the service will bind to.
    :param port: The port to listen to.
    :param publish_interval: The frequency, in frames, of publishing.

    If the address or the port is set to ``None``, the default value for a
    trajectory service is used. These defaults can be accessed as the
    ``narupa.trajectory.frame_server.DEFAULT_ADDRESS`` and
    ``narupa.trajectory.frame_server.DEFAULT_PORT`` constants.

    Publishing the frames can be activated, or deactivated, by setting the
    value of the :attr:`publishing_frames`, or by using the
    :meth:`make_publish_frames` and :meth:`make_not_publish_frames` methods.
    The publication is activated by default.
    """

    # TODO: Which the IMD, activation of the frame publishing should be coupled
    #  with activation of the IMD service. API may have to change.
    def __init__(
            self, simulation: Simulation, *,
            address: Optional[str] = None, port: Optional[int] = None,
            publish_interval: int = 1
    ):
        super().__init__(simulation)
        self._frame_server = FrameServer(address=address, port=port)
        self._frame_reporter = NarupaReporter(
            report_interval=publish_interval,
            frame_server=self._frame_server,
        )
        self.make_publish_frames()

    @classmethod
    def from_xml_input(
            cls, input_xml, *,
            address: Optional[str] = None, port: Optional[int] = None,
            publish_interval: int = 1
    ):
        """
        Create a runner from a serialized simulation.

        :param input_xml: Path to an XML serialised OpenMM simulation.
        :param address: The address the service will bind to.
        :param port: The port to listen to.
        :param publish_interval: The frequency, in frames, of publishing.
        :return: An instance of the class.

        If the address or the port is set to ``None``, the default value for a
        trajectory service is used. These defaults can be accessed as the
        ``narupa.trajectory.frame_server.DEFAULT_ADDRESS`` and
        ``narupa.trajectory.frame_server.DEFAULT_PORT`` constants.

        .. seealso::

            The XML serialized simulation can be produced by
            :func:`narupa.openmm.serializer.serialize_simulation`.

        """
        with open(str(input_xml)) as infile:
            simulation = deserialize_simulation(infile.read())
        return cls(
            simulation,
            address=address, port=port,
            publish_interval=publish_interval,
        )

    @property
    def trajectory_port(self):
        return self._frame_server.port

    def make_publish_frames(self):
        """
        Activate the publication of the frames.
        """
        if not self.publishing_frames:
            self.simulation.reporters.append(self._frame_reporter)

    def make_not_publish_frames(self):
        """
        Deactivate the publication of the frame.
        """
        if self.publishing_frames:
            self.simulation.reporters.remove(self._frame_reporter)

    @property
    def publishing_frames(self):
        """
        Returns ``True`` if the publication of the frames is activated.
        """
        return self._frame_reporter in self.simulation.reporters

    @publishing_frames.setter
    def publishing_frames(self, value: bool):
        """
        Activate or deactivate the publication of the frames.
        """
        if value:
            self.make_publish_frames()
        else:
            self.make_not_publish_frames()

    def close(self):
        """
        Close the network connection.
        """
        self._frame_server.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
