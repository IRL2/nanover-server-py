import json

import pytest

from essd.servicehub import ServiceHub, SERVICE_NAME_KEY, SERVICE_ADDRESS_KEY


@pytest.fixture
def properties():
    properties = {
        SERVICE_NAME_KEY: 'test service',
        SERVICE_ADDRESS_KEY: '127.0.0.1',
        "services": {
            "trajectory": 54321,
            "imd": 54322,
            "multiplayer": 54323,
            "builder": 54324
        },
        "essd_version": "1.0"
    }
    return properties


def test_service_message(properties):
    service = ServiceHub(**properties)
    assert service.message == json.dumps(properties)


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
