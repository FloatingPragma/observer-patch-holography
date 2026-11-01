from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import inversion as inv


def test_interval_reduction_and_inverse_residual_are_certified_without_float_solver():
    result = inv.matrix_certificate()
    assert F(result['inverse_infinity_norm_upper']) < 2
    assert F(result['inverse_times_true_V_minus_identity_infinity_norm_upper']) < F(1,10**18)
    assert result['floating_linear_solver_used'] is False
    for x in [F(-321),F(27,10),F(1111,3)]:
        c,s = inv.trig_interval(x-F(1,100000),x+F(1,100000))
        assert float(F(c.lo,inv.iv.SCALE)) <= math.cos(float(x)) <= float(F(c.hi,inv.iv.SCALE))
        assert float(F(s.lo,inv.iv.SCALE)) <= -math.sin(float(x)) <= float(F(s.hi,inv.iv.SCALE))


def test_unknown_weights_are_recovered_from_counts_on_the_declared_support():
    # A deliberately different distribution verifies that source weights are
    # not silently substituted by the reconstruction.
    weights = np.arange(1,9,dtype=float)/36
    frequencies = np.array([float(sum(map(F,row['omega_interval']))/2) for row in inv.source_rows()])
    signal = 1/16+15/16*np.exp(-1j*np.outer(inv.TIMES,frequencies))@weights
    trials = 1 << 40
    first = [round(float((1+x.real)/2)*trials) for x in signal]
    second = [round(float((1-x.imag)/2)*trials) for x in signal]
    result = inv.reconstruct(first,second,trials,F(1,10**12),F(1,10**12))
    for expected,row in zip(weights,result['weights'],strict=True):
        lo,hi = map(F,row['weight_interval'])
        assert float(lo) < expected < float(hi)
    assert result['all_eight_weights_strictly_positive']
    assert result['source_weights_used_in_reconstruction'] is False


def test_bad_data_are_not_clipped_into_a_successful_identification():
    result = inv.reconstruct([0]*8,[0]*8,100,0,0)
    assert not result['all_eight_weights_strictly_positive']


@pytest.mark.parametrize('first,second,trials,err', [([0]*7,[0]*8,10,0),([11]*8,[0]*8,10,0),([0]*8,[0]*8,0,0),([0]*8,[0]*8,10,.1),([0]*8,[0]*8,10,-1)])
def test_inexact_or_impossible_sampling_inputs_are_rejected(first,second,trials,err):
    with pytest.raises(ValueError): inv.reconstruct(first,second,trials,err,0)
