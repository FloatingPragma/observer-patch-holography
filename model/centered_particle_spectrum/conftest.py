import pytest
from . import packet

@pytest.fixture(scope='session')
def completed_spectrum():
    return packet._produce()
