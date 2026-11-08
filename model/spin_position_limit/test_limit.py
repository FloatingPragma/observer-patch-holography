"""Admission, exact error bounds and the source finite-section control."""
from fractions import Fraction as F
import numpy as np
import pytest
from .bounds import probability_error
from .finite import controls, _inverse_root


def test_actual_source_finite_compressions():
    report=controls()
    assert len(report['actual_source_character_weights'])==8
    assert F(report['certified_positive_character_gap'])>F(1,100)
    assert all(row['full_rank_agrees_with_original_position'] for row in report['finite_polar_compressions'])
    assert report['compression_does_not_commute_with_inverse_square_root']


@pytest.mark.parametrize('gap,err',[(F(),F()),(F(1,10),F(1,10)),(F(1,10),F(1,9))])
def test_no_inverse_without_certified_gram_gap(gap,err):
    with pytest.raises(ValueError):probability_error(gap,err,F(1),F())


@pytest.mark.parametrize('value',[True,0.1,-1,F(-1,3)])
def test_reject_nonexact_or_negative_errors(value):
    with pytest.raises(ValueError):probability_error(F(1),value,F(1),F())


def test_zero_errors_and_explicit_detector_error():
    a=probability_error(F(1,100),F(),F(2),F(),F(1,17))
    assert a['original_probability_error']==0
    assert a['normalized_input_event_error_upper_bound']==F(1,17)
    assert not a['division_by_success_probability']


def test_event_error_bound_against_noncommuting_matrix_perturbation():
    c=np.array([[2.,.4],[.4,1.5]])
    cp=c+np.array([[.01,.004],[.004,-.007]])
    a=np.array([[.3,.12],[.12,.6]])
    z=np.array([.3+.2j,.1-.1j]);zp=z+np.array([.001j,-.002])
    root=_inverse_root(c);rootp=_inverse_root(cp)
    actual=abs(np.vdot(zp,rootp@a@rootp@zp)-np.vdot(z,root@a@root@z))
    bound=probability_error(F(1),F(1,50),F(1),F(1,100))
    assert np.linalg.norm(cp-c)<.02
    assert np.linalg.norm(zp-z)<.01
    assert actual<float(bound['original_probability_error'])


def test_weak_statistics_do_not_control_position_variance():
    # Probability1/N at distanceN tends weakly to no tail, but contributesN
    # to the second moment. This rules out a variance claim from our limits.
    n=1000
    probability=F(1,n)
    assert probability<F(1,100)
    assert probability*n*n==n


def test_bounded_source_scores_need_amplitude_scaled_accuracy():
    from .bounds import bounded_probe_errors
    first=bounded_probe_errors(3,F(1,10),F(1,1000000),F(1,1000000))
    second=bounded_probe_errors(3,F(1,20),F(1,1000000),F(1,1000000))
    assert second['gram_operator_error']==4*first['gram_operator_error']
    assert second['overlap_vector_error']==2*first['overlap_vector_error']
    assert not first['division_by_observed_acceptance']


@pytest.mark.parametrize('n,s',[(0,F(1)),(True,F(1)),(2,F()),(1,0.5)])
def test_bounded_probe_admission(n,s):
    from .bounds import bounded_probe_errors
    with pytest.raises(ValueError):bounded_probe_errors(n,s,F(),F())


def test_effect_approximation_is_counted_separately():
    exact=probability_error(F(1,2),F(),F(2),F())
    approx=probability_error(F(1,2),F(),F(2),F(),effect_error=F(1,100))
    assert approx['original_probability_error']-exact['original_probability_error']==F(2,25)
