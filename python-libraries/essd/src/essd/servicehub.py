"""
Module defining a Service.
"""
import json
from socket import socket, AF_INET, SOCK_DGRAM, error

MAXIMUM_MESSAGE_SIZE = 1024


def get_default_ip():
    """
    Portable method for getting the default IP address of the machine.
    If no network connection is available, the loopback IP is returned.
    Note that if the machine is running multiple network interfaces, then this picks the 'primary' IP.

    :return: The primary IP address of the machine.
    """
    s = socket(AF_INET, SOCK_DGRAM)
    try:
        s.connect(('10.255.255.255', 1))
        ip = s.getsockname()[0]
    except error:
        ip = '127.0.0.1'
    finally:
        s.close()
    return ip


def construct_message(payload) -> str:
    return json.dumps(payload)


def validate_field(properties, key):
    if key not in properties:
        raise KeyError(f"Service does not contain the required field: {key}")


def validate_service(properties):
    validate_message_length(properties)
    for field in [SERVICE_NAME_KEY, SERVICE_ADDRESS_KEY]:
        validate_field(properties, field)


def validate_message_length(properties):
    message = construct_message(properties).encode()
    if len(message) > MAXIMUM_MESSAGE_SIZE:
        raise ValueError(f"Service definition exceeds the maximum message size of {MAXIMUM_MESSAGE_SIZE}")


class ServiceHub:
    """
    A definition of a ServiceHub that can be discovered or broadcast.

    A service hub consists of properties that must at least consist of a name and ip address.
    The payload can optionally include additional information on the services provided.
    """

    def __init__(self, **properties):
        validate_service(properties)
        self.properties = properties

    @property
    def name(self):
        return self.properties[SERVICE_NAME_KEY]

    @property
    def address(self):
        return self.properties[SERVICE_ADDRESS_KEY]

    @property
    def message(self):
        return construct_message(self.properties)

    @classmethod
    def from_json(cls, payload):
        properties = json.loads(payload)
        return cls(**properties)

    def __repr__(self):
        return f'{self.name}:{self.address}'

    def __str__(self):
        return str(self.properties)

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError(f"Cannot compare {self.__class__} with {other.__class__}")
        return hash(self) == hash(other)


SERVICE_NAME_KEY = "name"
SERVICE_ADDRESS_KEY = "address"