from copy import deepcopy
from fractions import Fraction as F
import json

import numpy as np
import pytest
from scipy.linalg import expm

from . import algebra, bounds as b, packet


@pytest.mark.parametrize('value', [True, 0.5, -1, '1/2', None])
def test_rejects_inexact_or_invalid_budget(value):
    with pytest.raises(ValueError):
        b.finite_history(value, 2, [F(1)])


@pytest.mark.parametrize('sites', [True, 0, -1, F(2), 2.0])
def test_rejects_invalid_carrier(sites):
    with pytest.raises(ValueError):
        b.finite_history(F(1), sites, [F(1)])


def test_both_directions_count_and_no_success_normalization():
    result = b.finite_history(F(1, 100), 8, [F(3), F(3)],
                              input_error=F(1, 100), native_error=F(1, 100))
    assert result['absolute_duration'] == 6
    assert result['clock_halftrace_upper'] == F(6, 25)
    assert result['complete_halftrace_upper'] == F(13, 50)
    with pytest.raises(ValueError):
        b.finite_history(F(1), 2, [F(-1)])


def test_zero_time_and_zero_field_are_exact():
    assert b.finite_history(F(), 100, [F(1000)])['complete_halftrace_upper'] == 0
    assert b.finite_history(F(1), 2, [])['complete_halftrace_upper'] == 0
    result = b.local_clock_plan(support_radius=0, support_sites=1, coupling=F(1),
                                duration=F(), interventions=3, error=F(1, 100))
    assert result['total_halftrace_upper'] == 0


def test_energy_bound_uses_HJ_range_and_own_ground():
    result = b.ground_energy_bound(F(3, 8), 2, 1, prepared_halftrace_error=F(1, 100))
    assert result['prepared_HJ_energy_excess_upper'] == F(19, 25)
    assert result['HJ_variance_for_exact_HB_ground_upper'] == F(9, 64)
    assert result['HJ_vacuum_selected'] is False


@pytest.mark.parametrize('duration', [F(1, 100), F(1, 4), F(3)])
def test_local_plan_against_existing_complete_history_bound(duration):
    from model.native_spin_observer.bounds import history_cutoff
    plan = b.local_clock_plan(support_radius=1, support_sites=3, coupling=F(1),
                             duration=duration, interventions=5, error=F(1, 100000))
    # The two halftrace cutoff contributions equal one unhalved cutoff bound.
    exact = history_cutoff(plan['radius'], 1, 3, F(1), duration, 5)
    assert exact <= plan['two_cutoffs_halftrace_upper']
    assert exact+plan['field_upper']*plan['sites']*duration/2 <= F(1, 100000)


def test_growing_time_uses_logarithmic_arithmetic():
    plan = b.local_clock_plan(support_radius=10**50, support_sites=1, coupling=F(1),
                             duration=F(10**100), interventions=7, error=F(1, 10**300))
    assert plan['total_halftrace_upper'] <= F(1, 10**300)
    assert plan['radius'].bit_length() < 400
    assert 0 < plan['field_upper'] < F(1, 10**300)
    assert plan['infinite_ground_marginal_error_supplied'] is False
    assert plan['global_controlled_spectrometer_covered'] is False


def test_source_service_halftrace_convention():
    assert b.spin_port_error(10, 20) == F(290, 2**20)
    assert b.spin_port_error(10, 1) == 1


def test_exact_ground_clock_countercontrol():
    result = algebra.controls()
    assert result['HJ_excitation_mean'] == '1/10'
    assert result['complete_return_record_TV'] == '9/25'
    assert result['full_state_halftrace_at_pi'] == '3/5'
    assert result['HB_ground_is_HJ_stationary'] is False


def test_complete_adaptive_history_including_failure_and_entangled_reference():
    # Two spins plus a spectator. Both clocks receive the same pure entangled
    # input; a complete three-outcome instrument controls the second unitary.
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    z = np.diag([1, -1]).astype(complex)
    h0 = sum(np.kron(a, a) for a in (x, y, z))/4-np.eye(4)/4
    hb = h0-F(3, 8)*(np.kron(x, np.eye(2))-np.kron(np.eye(2), x))/2
    hb = np.asarray(hb, complex)
    vector = np.array([1, 2j, -3, 1, 2, -1j, 1j, 2], complex)
    vector /= np.linalg.norm(vector)
    rho = np.outer(vector, vector.conj())
    p0, p1 = (np.eye(2)+z)/2, (np.eye(2)-z)/2
    kraus = [np.sqrt(4/5)*np.kron(p0, np.eye(2)),
             np.sqrt(4/5)*np.kron(p1, np.eye(2)), np.eye(4)/np.sqrt(5),
             np.zeros((4, 4))]
    assert np.linalg.norm(sum(k.conj().T@k for k in kraus)-np.eye(4)) < 1e-14
    def run(h):
        first = np.kron(expm(-1j*h/5), np.eye(2))
        state = first@rho@first.conj().T
        out = []
        for index, k in enumerate(kraus):
            k = np.kron(k, np.eye(2))
            direction = -1 if index == 2 else 1
            second = np.kron(expm(-1j*direction*h/7), np.eye(2))
            out.append(second@k@state@k.conj().T@second.conj().T)
        return out
    left, right = run(hb), run(h0)
    assert abs(sum(np.trace(a) for a in left)-1) < 1e-13
    assert np.trace(left[2]).real == pytest.approx(1/5)
    assert np.linalg.norm(left[3]) == 0
    error = sum(np.abs(np.linalg.eigvalsh(a-c)).sum()/2 for a, c in zip(left, right))
    bound = b.finite_history(F(3, 8), 2, [F(1, 5), F(1, 7)])
    assert 1e-5 < error <= float(bound['complete_halftrace_upper'])


@pytest.mark.parametrize('field,value', [
    ('same_prepared_state_complete_clock_comparison', False),
    ('effective_infinite_ground_selection_modulus', True),
    ('finite_field_ground_is_zero_field_vacuum', True),
    ('native_pulse_tape_executed_here', True),
    ('sharp_mass_shell_or_relativistic_rest_mass', True),
])
def test_forged_claims_rejected(field, value):
    result = packet.produce()
    result['claims'][field] = value
    with pytest.raises(ValueError):
        packet.validate(result)


def test_forged_records_and_source_bindings_rejected():
    result = packet.produce()
    for mutation in ('drop_failure', 'wrong_clock', 'source', 'bool_as_int'):
        bad = deepcopy(result)
        if mutation == 'drop_failure':
            bad['exact_two_spin_control']['HJ_return_probabilities'] = ['1']
        elif mutation == 'wrong_clock':
            bad['exact_two_spin_control']['HJ_first_gap'] = '9/8'
        elif mutation == 'source':
            key = next(iter(bad['source_bindings']))
            bad['source_bindings'][key] = '0'*64
        else:
            bad['all_pass'] = 1
        with pytest.raises(ValueError):
            packet.validate(bad)


def test_explicit_receipt_cli_path(tmp_path, monkeypatch):
    path = tmp_path/'exact.json'
    monkeypatch.setattr('sys.argv', ['clock', '--write', str(path)])
    packet.main()
    monkeypatch.setattr('sys.argv', ['clock', '--check', str(path)])
    packet.main()
    assert json.loads(path.read_text())['schema'] == 'oph.native_spin_clock_join.v1'


def test_explicit_receipt_cli_rejects_forged_claim(tmp_path, monkeypatch):
    path = tmp_path/'forged.json'
    value = packet.produce()
    value['claims']['finite_field_ground_is_zero_field_vacuum'] = True
    path.write_text(json.dumps(value))
    monkeypatch.setattr('sys.argv', ['clock', '--check', str(path)])
    with pytest.raises(ValueError):
        packet.main()
