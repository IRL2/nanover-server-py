"""
Module defining a Service.
"""

import json
from typing import Optional, Tuple
from uuid import uuid4
import nanover.essd

MAXIMUM_MESSAGE_SIZE = 1024
SERVICE_NAME_KEY = "name"
SERVICE_ADDRESS_KEY = "address"
SERVICE_ID_KEY = "id"
SERVICE_SERVICES_KEY = "services"
ESSD_VERSION_KEY = "essd_version"


class ServiceHub:
    """
    A definition of a ServiceHub that can be discovered or broadcast.

    A service hub consists of properties that must at least consist of a name and ip address.
    The payload can optionally include additional information on the services provided.

    :param: name: The name of the service hub.
    :param: address: The address of the service hub.
    :param: id: The unique ID of the service hub. If not specified, it will be generated.
    :param: essd_version: The version of ESSD this service hub uses. If not specified it will be determined automatically.
    :param: services: Dictionary of service names and their ports. Standard NanoVer services include
        imd, trajectory, multiplayer and builder.

    Example
    =======

    >>> # the ID field is usually autogenerated, provided here for completeness.
    >>> hub = ServiceHub(name="Example NanoVer Service Hub", address="localhost", id="12345")
    >>> hub.add_service("imd", 54322)
    >>> hub.add_service("trajectory", 54323)
    >>> hub.services
    {'imd': 54322, 'trajectory': 54323}
    >>> hub.message
    '{"name": "Example NanoVer Service Hub", "address": "localhost", "id": "12345", "essd_version": "1.0.0", "services": {"imd": 54322, "trajectory": 54323}}'

    The IP address of a service can either be a specific IP address of the interface to be broadcast, or it can be
    one of two special values: `localhost` or `[::]`.

    If `localhost` is specified, it will be broadcasted as running at 127.0.0.1, which is the usual translation of
    such a definition.

    If `[::]` is specified, then the appropriate IP address will be broadcast for each interface on the system.

    Example
    =======

    Consider a system with two interfaces, with IP addresses 192.168.1.2 (e.g. ethernet), and 72.34.5.5 (e.g. wireless),
    and broadcast addresses 192.168.1.255 and 72.34.255.255 respectively.
    The address [::] is provided. Then the service will be broadcast at as being at 192.168.1.2 on the ethernet
    network, and as being at 72.34.5.5 on the second network.
    Thus, a client at on the ethernet network will receive the address 192.168.1.2, which it can route to.
    """

    def __init__(self, **properties):
        _validate_service(properties)
        self.properties = properties
        if SERVICE_ID_KEY not in properties:
            self.properties[SERVICE_ID_KEY] = str(uuid4())
        if ESSD_VERSION_KEY not in properties:
            self.properties[ESSD_VERSION_KEY] = nanover.essd.__version__
        if SERVICE_SERVICES_KEY not in properties:
            self.properties[SERVICE_SERVICES_KEY] = {}

    @classmethod
    def from_json(cls, json_properties):
        """
        Constructs an instance of :class:`ServiceHub` from the given json string.

        :param json_properties: The JSON string containing the properties of the ServiceHub
        :return: An instance of :class:`ServiceHub`

        :raises`KeyError`: if the properties do not contain required fields, name and address.
        """
        properties = json.loads(json_properties)
        return cls(**properties)

    @property
    def name(self):
        """
        The name of the service hub.

        :return: The name of the service hub.
        :raises KeyError: if name has not been set.
        """
        return self.properties[SERVICE_NAME_KEY]

    @property
    def address(self):
        """
        The IP address of the service hub.

        :return: The IP address of the service hub.
        :raises KeyError: if name has not been set.
        """
        return self.properties[SERVICE_ADDRESS_KEY]

    @property
    def id(self):
        """
        Gets the unique ID string of this service hub.

        :return: The unique ID string of this service hub.
        """
        return self.properties[SERVICE_ID_KEY]

    @property
    def version(self):
        """
        Gets the version of ESSD this service hub is compatible with.

        :return: Version string of this service hub.
        """
        return self.properties[ESSD_VERSION_KEY]

    @property
    def services(self):
        """
        Gets the services registered at this service hub.

        :return: Dictionary of service definitions.
        """
        return self.properties[SERVICE_SERVICES_KEY]

    @property
    def message(self):
        """
        Returns the message that represents this service hub as a JSON string.

        :return: The message representing this service hub, as a JSON string.
        """
        return _construct_message(self.properties)

    def add_service(self, name, port):
        """
        Adds a service with the given name and port to the service hub definition.

        :param name: Name of the service
        :param port: Port at which the service is running

        """
        self.services[name] = port

    def to_message(self, override_address: Optional[str] = None) -> str:
        """
        Returns the JSON message representing this service hub, with the option to override this address.

        :param override_address: The address to override in the resulting message.
        :return: JSON message representing this service hub.

        When broadcasting services, it is useful to be able to provide human-readable shortcuts for underlying addresses,
        but clients receiving broadcasting need to know the actual address the service hub is running at.

        A typical use case is when using the '[::]' notation for defining a gRPC service, which means it will
        listen on all interfaces. When it comes to broadcasting, the discovery server will broadcast on all interfaces
        with the correct address for that interface.

        >>> hub = ServiceHub(name='Example', address='[::]', id='12345')
        >>> # The resulting message is not particularly useful.
        >>> hub.message
        '{"name": "Example", "address": "[::]", "id": "12345"}'
        >>> # Instead, at transmission, override with an actual address for the target interface.
        >>> hub.to_message(override_address="192.168.1.15")
        '{"name": "Example", "address": "192.168.1.15", "id": "12345"}'
        """
        if override_address is not None:
            hub = ServiceHub(**dict(self.properties))
            hub.properties[SERVICE_ADDRESS_KEY] = override_address
        else:
            hub = self
        return hub.message

    def get_service_address(self, service_name: str) -> Optional[Tuple[str, int]]:
        """
        Gets the address and port of a service, if it exists.
        :param service_name: Service name
        :return: Tuple consisting of address and port of a service, or None if not found.
        """
        if service_name not in self.services:
            return None
        else:
            return self.address, self.services[service_name]

    def __repr__(self):
        return f"{self.__class__.__name__}(**{self.properties})"

    def __str__(self):
        return str(self.properties)

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError(f"Cannot compare {self.__class__} with {other.__class__}")
        return hash(self) == hash(other)


def _construct_message(payload) -> str:
    return json.dumps(payload)


def _validate_field(properties, key):
    if key not in properties:
        raise KeyError(f"Service does not contain the required field: {key}")


def _validate_service(properties):
    _validate_message_length(properties)
    for field in [SERVICE_NAME_KEY, SERVICE_ADDRESS_KEY]:
        _validate_field(properties, field)


def _validate_message_length(properties):
    message = _construct_message(properties).encode()
    if len(message) > MAXIMUM_MESSAGE_SIZE:
        raise ValueError(
            f"Service definition exceeds the maximum message size of {MAXIMUM_MESSAGE_SIZE}"
        )