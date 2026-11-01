import pytest

from . import packet


@pytest.fixture(scope='session')
def readout():
    return packet._produce()
