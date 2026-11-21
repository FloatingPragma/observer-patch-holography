import numpy as np

from .. import formation


def test_full_carrier_write_leaves_every_other_channel_alone():
    chart = (2,1,3,2)
    n = int(np.prod(chart))
    rnd = np.random.default_rng(25813)
    incoming = rnd.normal(size=(n,8,n,8)) + 1j*rnd.normal(size=(n,8,n,8))
    copy = incoming.copy()
    remaining, written = formation.first_writes(incoming.reshape((*chart,8,*chart,8)))
    target = copy.copy()
    expected = np.zeros(n,complex)
    for site in range(n):
        expected[site] = copy[site,5,site,4]/np.sqrt(2)
        target[site,5,site,4] = expected[site]
    np.testing.assert_allclose(remaining.reshape(target.shape),target,rtol=0,atol=5e-16)
    np.testing.assert_allclose(written.ravel(),expected,rtol=0,atol=5e-16)
    np.testing.assert_array_equal(incoming,copy)


def test_disabled_strip_vertex_keeps_incident_column_and_explicit_zero_branch():
    chart = (1,1,3,1,8)
    state = np.zeros(chart+chart,dtype=int)
    state[0,0,1,0,5,0,0,1,0,4] = 1
    free, branches = formation.capture(state, angle=0.)
    np.testing.assert_array_equal(free,state)
    assert [z for z,_ in branches] == [1]
    assert np.count_nonzero(branches[0][1]) == 0
