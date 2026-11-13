from fractions import Fraction as F
from . import source as s

def test_complete_curvature_gauge_and_action_identities():
    for seed in (1,2,4):
        h=s.fixture(seed);G=s.einstein(h)
        assert all(not s.divergence(G,b)for b in range(4))
        other=s.fixture(seed+3)
        assert s.pairing(h,s.einstein(other))==s.pairing(other,G)
    xi={b:{s.parity((b,)):F(b+1)}for b in range(4)}
    assert all(not x for x in s.einstein(s.gauge(xi)).values())

def test_wrong_location_and_omitted_component_are_rejected():
    import pytest
    h=s.fixture(1);h[0,1]=h[1,0]={(0,0,0,0):F(1)}
    with pytest.raises(AssertionError):s.einstein(h)
    G=s.einstein(s.fixture(2));G[0,1]={}
    assert s.divergence(G,1)

def test_checkerboard_and_native_nearest_stencil():
    r=s.run()
    assert all(x['constant_kernel_only']for x in r['complete_periodic_controls'])
    assert r['exact_identities']['native_profile_Laplacian_equals_original_nearest_neighbour']
    assert not r['scope']['co_located_nine_null_observer_instrument_derived']
