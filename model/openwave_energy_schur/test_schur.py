from fractions import Fraction as F
import json
import numpy as np
import pytest
from .schur import HERE, patch, tail


def test_shifted_resolvent_bound_against_independent_dense_complement():
    # A genuine noncommuting complement, including negative spectral parameters.
    D = np.diag([89., 97., 113.])
    V = np.array([[2., 3., 0.], [3., -1., 2.], [0., 2., 1.]])
    T = np.array([[.3, .2], [-.4, .1], [.2, -.3]])
    b2 = F(1); d = F(80); magnetic = F(9)
    for z in (-11, 0, 9, 57, 79):
        exact = T.T@np.linalg.solve(D-V-z*np.eye(3), T)
        gram = T.T@np.linalg.solve(D-z*np.eye(3), T)
        assert np.linalg.norm(exact-gram, 2) <= float(tail(b2, d, magnetic, z, 0))
        inverse = np.linalg.solve(D-z*np.eye(3), T)
        J = inverse.T@inverse
        bound = 9*np.max(np.sum(abs(J), axis=1))/(1-9/(89-z))
        assert np.linalg.norm(exact-gram, 2) <= bound


def test_domain_and_finite_word_tail():
    for z in (-10, 0, 57, 79):
        values = [tail(F(2), F(80), F(9), z, n) for n in range(5)]
        assert all(a > b > 0 for a, b in zip(values, values[1:]))
    with pytest.raises(AssertionError):
        tail(F(1), F(80), F(9), 80, 0)


def test_third_patch_incidence_and_nonuniformity():
    rows = [patch(n) for n in (1, 2, 3)]
    assert [r['elementary_faces'] for r in rows] == [6, 36, 108]
    assert rows[2]['plaquettes_disjoint_from_fixed_boundary_probe'] >= 60
    for r in rows:
        assert r['same_clock_electric_subdivision_identity']
        assert F(r['full_QHP_norm_squared_lower_exact']) == F(r['beta_exact'])**2*r['plaquettes_disjoint_from_fixed_boundary_probe']/4


def test_authentic_receipt_inertia_and_narrower_intervals():
    receipt = json.loads((HERE/'results.json').read_text())
    assert not receipt['uniform_global_offblock_bound_claimed']
    for row in receipt['records']:
        assert row['lower_endpoint']['negative_eigenvalues_lower'] == 1
        assert row['lower_endpoint']['negative_eigenvalues_upper'] == 1
        assert row['upper_endpoint']['negative_eigenvalues_lower'] >= 2
        assert row['gap_interval_width_improvement_factor'] > 1
        assert not row['words_beyond_degree_zero_executed']
        assert F(row['shifted_bound_exact']) < F(row['old_unshifted_bound_exact'])
