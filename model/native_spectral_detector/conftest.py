import pytest


@pytest.fixture(scope='session')
def native_packet():
    from . import packet
    return packet._produce()
