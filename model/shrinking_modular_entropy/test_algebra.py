import numpy as np
import pytest
import sympy as s
from . import algebra,controls


def test_exact_full_Fock_and_observer_matrix():
    c=algebra.certificate()
    assert all(c['exact_controls'].values())
    assert c['observer_response_inverse_norm']=='68/5'


def test_wrong_hole_sign_is_detected():
    e=algebra.finite_example()
    wrong=e['A']+e['theta']*e['F']*e['F'].T
    assert s.simplify(wrong-e['actual_covariance'])!=s.zeros(2)


def test_deleting_exterior_support_breaks_rank_one_formula():
    # If g is internal, its centered number has an additional direct CAR term.
    e=algebra.finite_example();a=e['A'][0,0]
    c=algebra.annihilator(0);n=c.T*c
    true=s.trace((n-a*s.eye(16))*e['global_density']*n)
    assert true==a*(1-a)
    assert true!=-a*a


def test_all_finite_entropy_and_endpoint_controls():
    c=controls.certificate()
    assert c['entropy']['occupation_in_place_of_cell_norm_rejected']
    assert c['endpoints']['arbitrary_unweighted_probe_claim_rejected']
    assert c['projectors']['finite_low_frequency_pairs']==128


def test_quasifree_local_density_positive_at_both_endpoints():
    e=algebra.finite_example()
    for theta in (-1,1):
        rho=np.array(e['rho'].subs(e['theta'],theta),dtype=float)
        assert np.linalg.eigvalsh(rho).min()>=-1e-14
        assert abs(np.trace(rho)-1)<1e-14


def test_invalid_CAR_and_log_inputs():
    for args in ((True,4),(4,4),(0,False)):
        with pytest.raises(ValueError):algebra.annihilator(*args)
    with pytest.raises(ValueError):controls.hermitian_log(np.diag([0,1]))
