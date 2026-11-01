import pytest

from . import packet


@pytest.fixture(scope='session')
def ground_readout():
    return packet._produce()
