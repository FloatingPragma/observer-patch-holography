from fractions import Fraction as Q
from unittest.mock import patch

import numpy as np
import pytest

from . import certificate as c


def dense_poly(t, degree, P):
    out = np.zeros(P.shape, dtype=complex)
    power = np.eye(len(P))
    for re, im in c.coefficients(t, degree):
        out += complex(float(re), float(im))*power
        power = power@P
    return out


def test_binomial_expansion_and_exact_time_sign():
    assert c.coefficients(0, 4) == (c.ONE,)+(c.ZERO,)*4
    assert c.coefficients(Q(1, 2), 1) == ((Q(1), Q(-1, 2)), (Q(), Q(1, 2)))
    assert c.coefficients(Q(-2, 3), 6) == tuple(c.conjugate(z) for z in c.coefficients(Q(2, 3), 6))
    for x in (Q(), Q(1, 7), Q(1, 3)):
        for t in (Q(1, 2), Q(3, 2), Q(-2)):
            value = sum(complex(float(re), float(im))*float(x)**k for k,(re,im) in enumerate(c.x_coefficients(t, 10)))
            assert abs(value-complex(1-float(x))**complex(0,float(t))) <= float(c.tail_bound(t, 10))+1e-14


def test_coefficient_conversion_preserves_nontrivial_polynomial():
    P = np.array([[.8,.2],[.2,.8]])
    X = np.eye(2)-P
    powers = [np.eye(2)]
    for _ in range(5):powers.append(powers[-1]@X)
    expected = sum(complex(float(a),float(b))*Xk for (a,b),Xk in zip(c.x_coefficients(Q(2,3),5),powers))
    assert np.max(np.abs(expected-dense_poly(Q(2,3),5,P))) < 1e-13


def test_counted_words_include_zero_power_constraints_and_nonreturns():
    assert c.word_moment(2, (0,1), (0,0)) == 0
    assert c.word_moment(2, (0,), (0,)) == c.reference_entry(2,0,0,0)
    assert sum(c.word_moment(2,(a,),(0,)) for a in range(c.D)) == 1
    assert 0 < c.source.return_probability(2) < 1


def test_ordered_modular_correlations_match_independent_dense_evaluation():
    counts, denominator = c.source.source_counts()
    P = np.array(counts,dtype=float)/denominator
    vals, vecs = np.linalg.eigh(P)
    rho = (vecs*(vals**2/sum(vals**2)))@vecs.T
    indices, times = (0,1,2), (Q(-1,2),Q(),Q(1,2))
    product = np.eye(c.D,dtype=complex)
    for a,t in zip(indices,times):
        U=(vecs*np.exp(1j*float(t)*np.log(vals)))@vecs.T
        projector=np.zeros_like(P);projector[a,a]=1
        product=product@U@projector@U.conj().T
    exact=np.trace(rho@product)
    re,im=c.correlation(2,indices,times,14)
    assert abs(complex(float(re),float(im))-exact) <= float(c.approximation_bound(times,14))+1e-13


def test_floating_diagnostic_is_checked_but_does_not_change_receipt():
    counts, denominator = c.source.source_counts()
    eigenvalues, vectors = np.linalg.eigh(np.array(counts, dtype=float)/denominator)
    args = (2, (0, 6), (Q(), Q(1, 2)), 20)
    first, second = [], []
    expected = c.numerical_control(*args, first)
    with patch.object(c.np.linalg, 'eigh', return_value=(eigenvalues+1e-15, vectors)):
        assert c.numerical_control(*args, second) == expected
    assert first != second
    assert 'absolute_difference' not in expected
    assert expected['within_analytic_bound_plus_roundoff']
    with patch.object(c.np.linalg, 'eigh', return_value=(eigenvalues+0.1, vectors)):
        with pytest.raises(ValueError, match='spectral diagnostic disagrees'):
            c.numerical_control(*args)


def test_common_time_shift_and_loop_clock_conversion_are_exact():
    indices,times=(0,1),(Q(-1,3),Q(1,2))
    assert c.correlation(2,indices,times,10) == c.correlation(2,indices,tuple(t+Q(7,5) for t in times),10)
    assert c.coefficients(2*Q(1,4),12) == c.coefficients(4*Q(1,8),12)
    assert c.correlation(1,indices,times,10) != c.correlation(2,indices,times,10)


def test_record_and_matrix_stability_charge_normalization_and_amplification():
    times=(Q(),Q(1,2))
    assert c.record_error_bound(times,12,Q()) == 0
    assert c.record_error_bound(times,12,Q(1,10**12)) > Q(1,10**12)
    assert c.entry_error_bound(2,times,12,Q()) == 0
    with pytest.raises(ValueError,match='return lower bound'):
        c.record_error_bound(times,12,Q(1,66))
    with pytest.raises(ValueError,match='normalization'):
        c.entry_error_bound(2,times,12,Q(1,100))
    budget=c.observer_budget(2,times,20)
    assert budget['record_amplification_error_upper'] < budget['polynomial_bias_upper']
    assert budget['total_attempts_max'] > 20000


def test_every_proper_coordinate_cut_has_a_native_crossing():
    rows,_=c.source.source_counts()
    # Singletons and every initial segment exercise cuts without assuming
    # positivity of every matrix entry or complete-graph connectivity.
    for size in range(1,c.D):
        S=set(range(size))
        assert any(rows[i][j] for i in S for j in set(range(c.D))-S)
    assert len(c.source.matrix_unit_paths()) == c.D*c.D


def test_nested_sector_references_and_changed_completion_are_distinct():
    result=c.corners(2)
    assert result['sector_sizes'] == (55,11)
    assert result['nested_sector_size'] == 45
    assert sum(result['sector_weights']) == 1
    assert result['completed_reference_differs_from_original']
    assert result['two_step_compression_leakage_diagonal'] > 0
    assert result['cross_sector_reference_entry'] > 0
    assert result['within_sector_reference_entry'] > 0
    assert not result['general_tensor_region_no_go_claimed']


def test_completed_block_flow_intertwines_and_original_flow_does_not():
    rho=np.array([[float(c.reference_entry(2,0,i,j)) for j in range(c.D)] for i in range(c.D)])
    S=[i for i,state in enumerate(c.source.STATES) if state[0]==1]
    T=[i for i in range(c.D) if i not in S]
    bar=rho.copy();bar[np.ix_(S,T)]=0;bar[np.ix_(T,S)]=0
    v,V=np.linalg.eigh(bar);U=(V*np.exp(1j*.3*np.log(v)))@V.T
    sigma=rho[np.ix_(S,S)];sigma/=np.trace(sigma)
    s,W=np.linalg.eigh(sigma);Us=(W*np.exp(1j*.3*np.log(s)))@W.T
    A=np.zeros((len(S),len(S)));A[0,0]=1
    embedded=np.zeros_like(rho);embedded[np.ix_(S,S)]=A
    expected=np.zeros_like(U);expected[np.ix_(S,S)]=Us@A@Us.conj().T
    assert np.max(np.abs(U@embedded@U.conj().T-expected)) < 1e-12
    Qs=np.diag([int(i in S) for i in range(c.D)])
    assert np.linalg.norm(rho@Qs-Qs@rho)>0
    assert np.linalg.norm(bar@Qs-Qs@bar)==0


def test_tracial_and_one_dimensional_controls_have_trivial_modular_action():
    counts,denominator=c.source.source_counts()
    P=np.array(counts,dtype=float)/denominator
    A=np.diag([1]+[0]*(c.D-1))
    assert np.linalg.norm(P@A-A@P)>0
    trace=np.eye(c.D)/c.D
    assert np.linalg.norm(trace@A-A@trace)==0
    singleton=np.array([[float(c.reference_entry(2,0,0,0))]])
    singleton/=np.trace(singleton)
    assert np.array_equal(singleton,np.eye(1))


@pytest.mark.parametrize('bad',[True,1.0,'1'])
def test_exact_inputs_reject_coercion(bad):
    with pytest.raises(ValueError):c.coefficients(bad,3)
    with pytest.raises(ValueError):c.power_counts(bad)


def test_existing_source_bindings_are_checked_instead_of_resealed():
    old=c.strict_json.load
    def corrupt(path):
        value=old(path)
        value['source_sha256']['model/observer/native.py']='0'*64
        return value
    with patch.object(c.strict_json,'load',side_effect=corrupt):
        with pytest.raises(ValueError,match='inherited source binding changed'):
            c.inherited_pins()


def test_arbitrary_one_body_diagonals_cancel_the_interaction_contrast():
    from . import interaction
    eps=tuple(Q(i*i-3*i+2,7) for i in range(12))
    diagonal=tuple(Q(13,11)+sum(eps[i]*x for i,x in enumerate(state)) for state in c.source.STATES)
    assert interaction.contrast(diagonal)==0
    result=interaction.certificate(12)
    assert result['local_generator_distance_from_all_quadratics_lower']==Q(1,60)
    lo,hi=result['logarithm_contrast_interval']
    assert hi<0 or lo>0
    assert result['logarithm_operator_distance_from_all_quadratics_lower']>0
    assert result['all_native_sector_states_retained']==66


def test_exact_interaction_survives_more_terms_and_rejects_changed_source():
    from . import interaction
    a,b=interaction.certificate(8),interaction.certificate(16)
    assert max(a['logarithm_contrast_interval'][0],b['logarithm_contrast_interval'][0]) <= min(a['logarithm_contrast_interval'][1],b['logarithm_contrast_interval'][1])
    native=interaction.field.generator()
    changed=tuple(tuple(value+int(i==0 and j==0) for j,value in enumerate(row)) for i,row in enumerate(native))
    with patch.object(interaction.field,'generator',return_value=changed):
        with pytest.raises(ValueError,match='differs from native'):
            interaction.certificate()
