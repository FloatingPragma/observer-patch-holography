import numpy as np
import pytest

from oph_openwave.geometry import carrier, flat_ring, REFERENCE_FACE, RING_PORTS


def test_carrier_spectrum_and_geometric_ring():
    c = carrier()
    lap = c.incidence @ c.incidence.T
    expected = sorted([0] + [5-5**.5]*3 + [6]*5 + [5+5**.5]*3)
    np.testing.assert_allclose(np.linalg.eigvalsh(lap), expected, atol=1e-12)
    adjacency = np.diag(np.diag(lap)) - lap
    derived = {v for v in range(12) if v not in REFERENCE_FACE
               and sum(adjacency[v, p] for p in REFERENCE_FACE) == 2}
    assert derived == set(RING_PORTS)
    assert all(adjacency[a, b] == 0 for a in RING_PORTS for b in RING_PORTS)


def test_ring_connectivity_and_native_seam_count():
    for lx, ly in ((1, 1), (2, 3)):
        c = flat_ring(lx, ly)
        assert c.ports == 24*lx*ly and len(c.seams) == 63*lx*ly
        assert len({tuple(sorted(s)) for s in c.seams}) == len(c.seams)
        eigs = np.linalg.eigvalsh(c.incidence @ c.incidence.T)
        assert np.count_nonzero(abs(eigs) < 1e-10) == 1
    for bad in ((0, 1), (-1, 3), (1.5, 2)):
        with pytest.raises((ValueError, TypeError)):
            flat_ring(*bad)
