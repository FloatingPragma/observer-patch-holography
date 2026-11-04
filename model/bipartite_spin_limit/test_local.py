"""Independent source-matrix controls for the exact local averaging step."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from . import local
from .test_limit import jordan_wigner


def source_rotation(sites, edges, c, U, reverse=False):
    H, P, spin = jordan_wigner(sites, edges, c, U)
    # On Nparticle=sites, the integer onsite hole/doublon penalty is 2D.
    penalty = np.rint(2*np.diag(H)/U).astype(int)
    nu = U/2
    equal = penalty[:, None] == penalty[None, :]
    frequency = penalty[:, None]-penalty[None, :]
    G, Y = H.copy(), np.eye(len(H))
    controls = []
    for _ in range(local.STEPS):
        V = np.where(equal, 0, G)
        S = np.divide(V, nu*frequency, out=np.zeros_like(V), where=~equal)
        if reverse:
            S = -S
        step = expm(S)
        rotated = step@G@step.T
        controls.append((np.linalg.norm(V, 2), np.linalg.norm(np.where(equal, 0, rotated), 2)))
        assert np.linalg.norm(S+S.T, 2) < 1e-13
        G, Y = rotated, step@Y
    return H, P, spin, penalty, G, Y, controls


@pytest.mark.parametrize('sites,edges', [
    (2, ((0, 1, 1),)),
    (3, ((0, 1, -1), (1, 2, 1))),
    (4, ((0, 1, -1), (0, 2, 1), (1, 3, 1), (2, 3, 1))),
])
def test_all_penalty_sectors_and_exact_source_rotation(sites, edges):
    H, P, spin, N, G, Y, controls = source_rotation(sites, edges, 1, 32)
    assert set(N) >= {0, 2}
    assert np.linalg.norm(Y@Y.T-np.eye(len(Y)), 2) < 2e-12
    assert np.linalg.norm(G-Y@H@Y.T, 2) < 1e-10
    assert np.linalg.norm(N[:, None]*P, 2) == 0
    for before, after in controls:
        if before > 1e-10:
            assert after < before/2
    assert np.linalg.norm(P.T@G@P-spin, 2) < .001
    assert np.linalg.norm(P.T@G@P+spin, 2) > .01


def test_wrong_inverse_sign_preserves_neither_decoupling_nor_spin_energy():
    edges = ((0, 1, 1),)
    good = source_rotation(2, edges, 1, 16)
    wrong = source_rotation(2, edges, 1, 16, reverse=True)
    assert good[-1][0][1] < good[-1][0][0]/10
    assert wrong[-1][0][1] > wrong[-1][0][0]


def test_actual_dressed_local_observable_converges_to_interacting_spin_motion():
    edges = ((0, 1, -1), (1, 2, 1))
    errors = []
    for c in (4, 8, 16, 32):
        H, P, spin, N, G, Y, _ = source_rotation(3, edges, c, c*c)
        A = np.diag([1., 1., 1., 1., -1., -1., -1., -1.])
        bare = P@A@P.T
        dressed = Y.T@bare@Y
        initial = Y.T@P
        evolved = expm(-.3j*H)@initial
        actual = evolved.conj().T@dressed@evolved
        expected = expm(.3j*spin)@A@expm(-.3j*spin)
        errors.append(np.linalg.norm(actual-expected, 2))
        assert np.linalg.norm(expected-A, 2) > .1
    assert all(later < earlier*.4 for earlier, later in zip(errors, errors[1:]))


def test_source_spin_bonds_are_noncommuting():
    _, _, first = jordan_wigner(3, ((0, 1, 1),), 1, 1)
    _, _, second = jordan_wigner(3, ((1, 2, -1),), 1, 1)
    assert np.linalg.norm(first@second-second@first, 2) > .1


def test_conservative_scale_certificate_and_fixed_exchange_rate():
    result = local.local_bounds(2**256, 2**512, 1, 8)
    assert result['exchange'] == 1
    assert 0 < result['bare_observable_on_dressed_input_error'] < F(1, 10**6)
    doubled = local.local_bounds(2**257, 2**514, 1, 8)
    ratio = doubled['bare_observable_on_dressed_input_error']/result['bare_observable_on_dressed_input_error']
    assert F(49, 100) < ratio <= F(1, 2)
    assert not result['depends_on_total_volume']
    assert not result['undressed_input_comparison']
    assert not result['interleaved_repair_comparison']
    assert len(local.iteration_majorants()) == 5


@pytest.mark.parametrize('args', [(1, 1), (True, 2**40), (1., 2**40),
                                  (1, 2**40, -1), (1, 2**40, 0, True),
                                  (1, 2**40, 0, 0)])
def test_outside_local_theorem_domain_is_rejected(args):
    with pytest.raises(ValueError):
        local.local_bounds(*args)
