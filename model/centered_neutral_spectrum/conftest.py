import pytest
from . import packet

@pytest.fixture(scope='session')
def result():
    return packet.produce()
