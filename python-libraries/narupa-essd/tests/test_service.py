import json

import pytest

import narupa.essd
from narupa.essd.servicehub import (ServiceHub, SERVICE_NAME_KEY, SERVICE_ADDRESS_KEY, SERVICE_SERVICES_KEY,
                                    SERVICE_ID_KEY, ESSD_VERSION_KEY)
from narupa.essd.utils import get_broadcast_addresses


def get_broadcastable_ip():
    broadcast_addresses = get_broadcast_addresses()
    if len(broadcast_addresses) == 0:
        raise RuntimeError("No broadcastable IP addresses could be found on the system!")
    return broadcast_addresses[0]['addr']


@pytest.fixture
def properties():
    properties = {
        SERVICE_NAME_KEY: 'test service',
        SERVICE_ADDRESS_KEY: get_broadcastable_ip(),
        SERVICE_SERVICES_KEY: {
            "trajectory": 54321,
            "imd": 54322,
            "multiplayer": 54323,
            "builder": 54324
        },
        ESSD_VERSION_KEY: "1.0.0",
        SERVICE_ID_KEY: "12345"
    }
    return properties


def test_service_message(properties):
    service = ServiceHub(**properties)
    assert service.message == json.dumps(properties)


def test_version(properties):
    del properties['essd_version']
    service = ServiceHub(**properties)
    assert service.version == narupa.essd.__version__


def test_service_generate_uuid(properties):
    del properties['id']
    service = ServiceHub(**properties)
    assert 'id' in service.properties


def test_service_from_json(properties):
    message = json.dumps(properties)
    service = ServiceHub.from_json(message)
    assert service.properties == properties


def test_service_no_name(properties):
    del properties[SERVICE_NAME_KEY]
    with pytest.raises(KeyError):
        _ = ServiceHub(**properties)


def test_service_no_address(properties):
    del properties[SERVICE_ADDRESS_KEY]
    with pytest.raises(KeyError):
        _ = ServiceHub(**properties)


def test_service_too_long(properties):
    properties['some_key'] = '1' * 1024
    with pytest.raises(ValueError):
        _ = ServiceHub(**properties)


def test_repr(properties):
    hub = ServiceHub(**properties)
    # tests that creating a new hub from the evaluation of a representation of a hub creates an equal object.
    assert eval(repr(hub)) == hub


def test_equals_object(properties):
    hub = ServiceHub(**properties)
    x = object()
    with pytest.raises(TypeError):
        _ = hub == x


def test_add_service(properties):
    hub = ServiceHub(**properties)
    assert hub.services == properties['services']
    hub.add_service("test", 54322)
    properties['services'].update({"test": 54322})
    assert hub.services == properties['services']


def test_add_service_replacement(properties):
    hub = ServiceHub(**properties)
    hub.add_service("imd", 88888)
    properties['services'].update({"imd": 88888})
    assert hub.services == properties['services']