from copy import deepcopy
from fractions import Fraction as F
import json
from unittest import mock

import numpy as np
import pytest

from . import algebra, bounds, packet
from . import __main__ as cli


def test_source_single_bond_and_complex_dual_identity():
    result = algebra.exact_checks()
    assert result['native_single_bond_inverse'] == '1/4'
    assert result['native_single_bond_unit_tangent_cost'] == '1'
    assert result['complex_trial_dual_residual_identity_exact'] is True
    assert result['complete_bounded_sine_probe_exact'] is True
    assert result['bounded_probe_inverse_response_limit_exact'] is True


@pytest.mark.parametrize('channel', ('P', 'Q'))
@pytest.mark.parametrize('radius', (F(1, 8), F(1), F(64)))
def test_actual_continuum_response_has_positive_finite_cost(channel, radius):
    row = bounds.continuum_response(channel, F(1, 1000), radius)
    assert 0 < row['inverse_energy_lower'] <= row['inverse_energy_upper']
    assert 0 < row['minimum_tangent_energy_lower'] <= row['minimum_tangent_energy_upper']
    assert row['sharp_shell_or_speed_inferred'] is False


def test_scaling_of_relaxed_cost_is_preserved():
    for channel in ('P', 'Q'):
        base = bounds.continuum_response(channel, F(1, 10), 1)
        small = bounds.continuum_response(channel, F(1, 10), F(1, 4))
        multiplier = F(1) if channel == 'P' else F(1, 16)
        for key in ('minimum_tangent_energy_lower', 'minimum_tangent_energy_upper'):
            assert small[key] == multiplier*base[key]


@pytest.mark.parametrize('bad', (True, 0.5, '1/2', None))
def test_exact_parameters_reject_type_coercion(bad):
    with pytest.raises(ValueError):
        bounds.continuum_response('P', bad, 1)


@pytest.mark.parametrize('args', (('P', F(0), 1), ('Q', F(3, 4), 1), ('P', F(1, 4), 0), ('other', F(1, 4), 1)))
def test_invalid_continuum_hypotheses_fail(args):
    with pytest.raises(ValueError):
        bounds.continuum_response(*args)


def test_inverse_regularization_and_feasible_energy_on_nonsharp_response():
    weights = (F(1, 3), F(2, 3)); energies = (F(1), F(4))
    chi = sum(w/e for w, e in zip(weights, energies))
    previous = F(1)
    for eta in (F(1), F(1, 4), F(1, 16), F(1, 64)):
        ce = sum(w/(e+eta) for w, e in zip(weights, energies))
        pe = sum(w/(e+eta)**2 for w, e in zip(weights, energies))
        remainder = sum(w/(e*(e+eta)**2) for w, e in zip(weights, energies))
        assert chi == ce+eta*pe+eta*eta*remainder
        assert 0 <= remainder <= pe/min(energies)
        cost = (ce-eta*pe)/(4*ce*ce)
        assert F(1, 4)/chi <= cost <= previous
        previous = cost
        lower = bounds.regularized_response(1, 3, eta)
        assert ce >= lower['regularized_inverse_lower']
        assert pe >= lower['inverse_square_lower']


def test_dual_residual_certificate_without_gap():
    # Include a zero ground eigenvalue; positive eta alone controls the inverse.
    energy = np.diag([0., 1., 3.])
    a = np.array([0., 1., 2j])
    eta = 0.125
    trial = np.array([1+1j, 0.7, 0.2+0.4j])
    shifted = energy+eta*np.eye(3)
    exact = np.vdot(a, np.linalg.solve(shifted, a)).real
    residual = a-shifted@trial
    dual = (2*np.vdot(a, trial).real-np.vdot(trial, shifted@trial).real)
    assert 0 <= exact-dual <= np.vdot(residual, residual).real/eta
    rational = bounds.residual_certificate(F(-1), F(3), F(1, 8))
    assert rational['regularized_inverse_lower'] == 0
    assert rational['regularized_inverse_upper'] == 23
    assert rational['requires_ground_gap'] is False
    with pytest.raises(ValueError):
        bounds.residual_certificate(-1, 0, 1)


def test_complete_filter_budget_retains_time_tail_and_origin_boundary():
    row = bounds.filter_budget(F(1, 4), 1024, F(1, 4096), 28, F(1, 10**6))
    assert row['successful_operator_error_upper'] == F(1, 256)+F(7, 2048)+F(1, 10**6)
    assert row['absolute_ground_origin_used_by_filter'] is False
    assert row['native_compilation_error_included'] is False


def test_complete_score_does_not_normalize_by_success():
    row = bounds.finite_readout_error(F(1, 4), F(1, 1000), F(1, 100), F(1, 10000), F(1, 50))
    assert row['bounded_score_error_upper'] == F(341, 10000)
    assert row['far_label_probability_upper'] == F(1, 50)
    assert row['spectral_label_coupling_certificate_required'] is True
    assert row['all_reject_scores_zero'] is True
    assert row['renormalization_by_success_used'] is False


def test_nominal_spacing_does_not_bound_far_label_error():
    eta = F(1, 4)
    true_energy = F(1)
    ideal_score = eta/(true_energy+eta)
    far_label_score = eta/(F(10)+eta)
    leakage = F(1, 2)
    actual_score = (1-leakage)*ideal_score+leakage*far_label_score
    error = abs(actual_score-ideal_score)
    assert error > 0
    with_leakage = bounds.finite_readout_error(eta, 0, 0, far_label_probability=leakage)
    assert error <= with_leakage['bounded_score_error_upper']
    assert error > bounds.finite_readout_error(eta, 0, 0)['bounded_score_error_upper']
    with pytest.raises(ValueError):
        bounds.finite_readout_error(eta, 0, 0, far_label_probability=F(3, 2))


def test_bounded_source_probe_keeps_limit_order_and_signal_cost():
    row = bounds.bounded_probe_error(F(1, 4), F(1, 8), F(1, 10**8))
    assert row['known_score_factor'] == F(1, 256)
    assert row['regularized_inverse_error_upper'] == F(1, 390625)
    assert row['divides_by_observed_success'] is False
    assert row['continuum_limit_order'][-1] == 'regularizer_to_zero'
    assert row['effective_continuum_modulus_proved'] is False
    with pytest.raises(ValueError):
        bounds.bounded_probe_error(1, 0, 0)


@pytest.mark.parametrize('args', ((0, 1, 1, 1), (1, 0, 1, 1), (1, 1, -1, 1), (1, 1, 1, -1)))
def test_invalid_filter_budget_fails(args):
    with pytest.raises(ValueError):
        bounds.filter_budget(*args)


def test_retained_complete_source_intervals():
    payload = json.loads(packet.OUTPUT.read_text())
    rows = payload['complete_source_resolvent_response']['rows']
    assert len(rows) == 9
    for row in rows:
        lo, hi = map(F, row['regularized_inverse_interval'])
        clo, chi = map(F, row['minimum_unit_tangent_energy_interval'])
        success = tuple(map(F, row['resolvent_filter_success_interval']))
        rejected = tuple(map(F, row['all_complementary_reject_probability_interval']))
        assert 0 < lo <= hi and 0 < clo <= chi
        assert 0 < success[0] <= success[1] < 1
        assert rejected[0]+success[1] <= 1 <= rejected[1]+success[0]
        assert row['ground_mean_zero_from_exact_parity'] is True


@pytest.mark.parametrize('mutation', ('boolean', 'bound', 'drop_reject', 'source'))
def test_scope_or_complete_receipt_mutation_fails(mutation):
    baseline = json.loads(packet.OUTPUT.read_text())
    changed = deepcopy(baseline)
    if mutation == 'boolean':
        changed['scope']['first_energy_uniform_integrability_assumed'] = 0
    elif mutation == 'bound':
        changed['actual_continuum_response_bounds'][0]['inverse_energy_upper'] = '1/13'
    elif mutation == 'drop_reject':
        del changed['complete_source_resolvent_response']['rows'][0]['all_complementary_reject_probability_interval']
    else:
        changed['source_bindings'][next(iter(changed['source_bindings']))] = '0'*64
    with mock.patch.object(packet, 'build', return_value=baseline):
        with pytest.raises(ValueError):
            packet.validate(changed)


def test_current_bound_sources_exclude_docs_and_receipt_self_reference():
    sources = packet.bindings()
    assert sources and all(not x.startswith('docs/') for x in sources)
    assert 'model/spin_relaxed_response/results.json' not in sources
    assert 'model/native_ordered_ground/witness.npz' in sources
    assert 'model/native_ordered_ground/certificate.py' in sources


def test_owner_cli_checks_requested_receipt_path():
    with mock.patch('sys.argv', ['spin_relaxed_response', '--check', '/tmp/selected-response.json']):
        with mock.patch.object(packet, 'verify', return_value=True) as verify:
            cli.main()
    verify.assert_called_once_with('/tmp/selected-response.json')
