from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import identification as ident
from ..particle_identification import vacuum


def test_exact_sampling_matrix_is_invertible_with_controlled_noise_amplification():
    result = ident.spectral_resolution_certificate()
    rows = sorted(vacuum.certificate()['spectral_rows'], key=lambda r: F(r['omega_interval'][0]))
    frequencies = np.array([float(sum(map(F,r['omega_interval']))/2) for r in rows])
    v = np.exp(-1j*np.outer(result['sample_model_times'], frequencies))
    lo, hi = map(F, result['determinant_certified_interval'])
    assert float(lo)-1e-8 < abs(np.linalg.det(v)) < float(hi)+1e-8
    assert lo > 1400
    assert np.linalg.norm(np.linalg.inv(v),np.inf) < float(ident.INVERSE_INFINITY_UPPER)
    assert result['long_time_sampling_protocol_executed_here'] is False
    assert result['arbitrary_unknown_spectrum_identified_by_finitely_many_samples'] is False


@pytest.mark.parametrize('point', [F(-101,3), F(-1,10), F(), F(7,2), F(2345,19)])
def test_pi_reduction_encloses_sine_on_both_signs_and_many_periods(point):
    lo, hi = ident.absolute_sine_interval(point-F(1,10000), point+F(1,10000))
    assert float(lo) <= abs(math.sin(float(point))) <= float(hi)


def test_duplicate_frequency_negative_control_is_not_identifiable():
    assert ident.vandermonde_determinant_bounds([(F(1),F(1)),(F(1),F(1))],39)[0] == 0


@pytest.mark.parametrize('value', [-1, 0.1, float('nan'), True])
def test_inexact_or_negative_probability_errors_are_rejected(value):
    with pytest.raises(ValueError):ident.weight_error_bound(value)


def test_prospective_resolution_budget_keeps_accuracy_and_shots_as_unexecuted_premises():
    result = ident.prospective_sampling_budget()
    assert F(result['joint_failure_probability_strict_upper']) < F(1,100)
    assert result['total_shots'] == 77309411328
    assert result['requires_systematic_error_bound_at_all_16_settings'] is True
    assert result['sampling_budget_executed'] is False
    assert result['abort_records_are_retained_without_postselection'] is True
    assert F(result['reconstructed_weight_error_strict_upper']) < F(1,1000)
