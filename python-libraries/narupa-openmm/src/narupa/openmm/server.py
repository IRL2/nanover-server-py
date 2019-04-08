"""
Server to run an OpenMM simulation and publish its frames for Narupa.
"""

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
    :param address: The name of this host.
    :param port: The port to listen to.
    :param publish_interval: The frequency, in frames, of publishing.

    Publishing the frames can be activated, or deactivated, by setting the
    value of the :attr:`publishing_frames`, or by using the
    :meth:`make_publish_frames` and :meth:`make_not_publish_frames` methods.
    The publication is activated by default.
    """
    # TODO: The API is not satisfying:
    #  * Should it be possible to deactivate the Narupa reporter? Is it any
    #    useful?
    #  * Host name and port should have a default; that default might even be
    #    a dynamic port where the first available range in a range is used.
    def __init__(
            self, simulation: Simulation, *,
            address: str, port: int,
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
    def from_xml_input(cls, input_xml, *,
                       address, port, publish_interval: int = 1):
        """
        Create a runner from a serialized simulation.

        :param input_xml: Path to an XML serialised OpenMM simulation.
        :param address: The name of this host.
        :param port: The port to listen to.
        :param publish_interval: The frequency, in frames, of publishing.
        :return: An instance of the class.

        .. seealso::

            The XML serialized simulation can be produced by
            :fun:`narupa.openmm.serializer.serialize_simulation`.

        """
        with open(str(input_xml)) as infile:
            simulation = deserialize_simulation(infile.read())
        return cls(
            simulation,
            address=address, port=port,
            publish_interval=publish_interval,
        )

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
            self.make_quiet()
