"""Independent source, entropy, full-record and modular-control checks."""
from copy import deepcopy
from decimal import Decimal, localcontext
from fractions import Fraction as F
import math

import numpy as np
import pytest

from ..radial_interacting_core import algebra as source_algebra
from . import modular, records, regional


def test_actual_full_source_origin_coupling_and_ground_faithfulness():
    h, sites = source_algebra.source_matrix()
    matrix = np.array(h, dtype=complex)
    origin = sites.index((0, 0, 0))
    center = [2*origin, 2*origin+1]
    complement = [j for j in range(54) if j not in center]
    V = matrix[np.ix_(complement, center)]
    np.testing.assert_array_equal(V.conj().T@V, 1.5*np.eye(2))
    assert np.linalg.norm(matrix[np.ix_(complement, complement)], 2) < 66
    result = regional.ground_faithfulness()
    assert F(result['central_zero_weight_upper']) == F(2904, 2905)
    assert F(1, 4096) < F(result['actual_ground_parity_lower'])
    assert F(result['actual_ground_parity_upper']) < F(1, 1000)
    assert result['actual_ground_local_M4_faithful']
    assert not result['actual_ground_full_Fock_state_faithful']


@pytest.mark.parametrize('q', [F(1, 4), F(1, 3), F(1, 2), F(3, 4)])
def test_actual_source_family_uniform_positive_weights(q):
    for ng, nc in [(F(0), F(1)), (F(1, 2000), F(1999, 2000))]:
        g, c, eps = regional.family_weights(q, ng, nc)
        assert min(g, c) > F(6, 25) and g+c+eps == 1
        assert g*ng+c*nc+eps/2 == q
        gap, eg, einf = F(3, 10**7), F(-123), F(1, 17)
        if q < F(3, 4):
            delta = F(1, 10000)
            g2, c2, _ = regional.family_weights(q+delta, ng, nc)
            energy = g*eg+c*(eg+gap)+eps*einf
            energy2 = g2*eg+c2*(eg+gap)+eps*einf
            assert (energy2-energy)/delta == gap/(nc-ng)


@pytest.mark.parametrize('values', [(0.5, F(0), F(1)),
    (F(1, 8), F(0), F(1)), (F(1, 2), F(1, 1000), F(1)),
    (F(1, 2), F(0), F(999, 1000))])
def test_family_rejects_outside_certified_source_domain(values):
    with pytest.raises(ValueError):
        regional.family_weights(*values)


def test_full_site_entropy_and_arbitrary_coherence_variational_identity():
    D = np.diag([1., 0., 0., 1.])
    rng = np.random.default_rng(517)
    raw = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    omega = raw@raw.conj().T
    omega /= np.trace(omega)
    for r in (0.25, 0.375, 0.5, 0.75):
        tau = np.diag([r, 1-r, 1-r, r])/2
        kappa = math.log((1-r)/r)
        left = records.entropy(omega)-kappa*np.trace(D@omega).real
        right = records.entropy(tau)-kappa*r-records.relative_entropy(omega, tau)
        assert abs(left-right) < 2e-14
        assert abs(records.entropy(tau)-math.log(2)-records.binary_entropy(r)) < 2e-14
        for q in (0.25, 0.4, 0.5, 0.75):
            other = np.diag([q, 1-q, 1-q, q])/2
            divergence = records.relative_entropy(other, tau)
            assert 2*(q-r)**2-2e-14 <= divergence <= (8/3)*(q-r)**2+2e-14
            assert abs(records.entropy(other)-records.entropy(tau)-kappa*(q-r)+divergence) < 2e-14


@pytest.mark.parametrize('q', [F(1, 4096), F(1, 1000), F(1, 4), F(1, 2),
                             F(3, 4), F(4095, 4096)])
def test_logarithm_rational_enclosure_against_independent_decimal(q):
    with localcontext() as context:
        context.prec = 160
        expected = (Decimal((1-q).numerator)/Decimal((1-q).denominator)
                    /(Decimal(q.numerator)/Decimal(q.denominator))).ln()
        for terms in (0, 4, 48):
            lower, upper = modular.log_odds_interval(q, terms)
            lo = Decimal(lower.numerator)/Decimal(lower.denominator)
            hi = Decimal(upper.numerator)/Decimal(upper.denominator)
            assert lo <= expected <= hi
        assert upper-lower < F(1, 10**45)


def test_native_modular_complete_site_execution_and_half_parity_identity():
    result = modular.certificate()
    for row in result['native_examples']:
        assert row['native_pulses'] == 8 and row['native_first_failure_slots'] == 28
        assert row['full_site_rows'] == row['retained_reference_columns'] == 4
        assert F(row['complete_isometry_error_upper']) < F(1, 10**9)
        assert row['all_failure_weights_retained']
        assert not row['success_postselected_or_renormalized']
    assert modular.log_odds_interval(F(1, 2)) == (0, 0)
    assert modular.log_odds_interval(F(3, 4)) == tuple(-x for x in reversed(modular.log_odds_interval(F(1, 4))))


def test_modular_literal_rotation_mutation_is_rejected(monkeypatch):
    original = modular.circuits.emit
    def wrong_angle(template):
        for index, event in enumerate(original(template)):
            yield dict(event, angle=event['angle']+0.01) if index == 0 else event
    monkeypatch.setattr(modular.circuits, 'emit', wrong_angle)
    with pytest.raises(ArithmeticError, match='error budget'):
        modular.word_and_bound()


def test_refined_controls_cover_every_source_port_and_fixed_pi_pulse():
    coarse, fine = modular.refined_commands(bits=40), modular.refined_commands(bits=80)
    assert len(fine['exact_dyadic_port_commands']) == 8
    assert F(fine['complete_isometry_error_upper']) < F(coarse['complete_isometry_error_upper'])/10**10
    for command in fine['exact_dyadic_port_commands']:
        assert len(command['integrated_port_coefficients']) == 12
        assert max(map(F, command['port_error_bounds'])) <= F(1, 2**80)
    assert not fine['legacy_binary64_API_used_for_refined_commands']
    fixed_ports = modular.refined_commands(bits=160, terms=4)
    both_refined = modular.refined_commands(bits=160, terms=48)
    assert F(both_refined['complete_isometry_error_upper']) < F(fixed_ports['complete_isometry_error_upper'])/10**30
    with pytest.raises(ValueError):
        modular.word_and_bound(time=F(1))
    with pytest.raises(ValueError):
        modular.log_odds_interval(F(0))
    with pytest.raises(ArithmeticError):
        modular.word_and_bound(terms=0)


def test_complete_native_history_entropy_map():
    result = records.native_record_certificate()
    assert result['complete_histories_checked'] == 78*1024
    assert result['fixed_coin_and_pointer_native_map_injective']
    assert not result['extra_current_entropy_added']


def test_deleted_current_record_fails_input_recovery(monkeypatch):
    original = records.native.native
    def deleted_current(state, edge, coin):
        output, committed = original(state, edge, coin)
        return output, (committed[0], 0)
    monkeypatch.setattr(records.native, 'native', deleted_current)
    with pytest.raises(ArithmeticError, match='recover'):
        records.native_record_certificate()


def test_correlated_quantum_records_and_same_H_full_Fock_work():
    quantum = records.quantum_record_control()
    assert quantum['arbitrary_input_unread_coherence_change'] > 0.1
    result = records.chain_control()
    assert result['readout_work'] > 0
    assert result['relative_entropy_to_same_H_Gibbs'] > 0
    assert result['full_Fock_rows'] == 64
    assert not result['work_identified_as_heat']


def test_classical_site_record_entropy_direct_eight_probabilities():
    for q in (0.25, 0.5, 0.75):
        initial = np.array([q, 1-q, 1-q, q])/2
        conditional = np.array([[1, 49, 49, 1], [49, 1, 1, 49]])/50
        joint = conditional*initial
        output_entropy = -np.sum(joint*np.log(joint))
        input_entropy = -np.dot(initial, np.log(initial))
        assert abs(output_entropy-input_entropy-records.binary_entropy(1/50)) < 2e-14


def test_scope_and_manifest_mutations_fail():
    from . import packet
    value = packet.report()
    packet.validate(value, value)
    for field in ('actual_ground_full_Fock_faithfulness',
                  'finite_G3_or_G4_Einstein_premise_discharged',
                  'signed_kinetic_bias_identified_with_area'):
        bad = deepcopy(value); bad['scope'][field] = True
        with pytest.raises(ValueError):
            packet.validate(bad, value)
    bad = deepcopy(value); bad['source_bindings'].pop(next(iter(bad['source_bindings'])))
    with pytest.raises(ValueError):
        packet.validate(bad, value)
    bad = deepcopy(value); bad['records']['native_records']['global_complete_entropy_increment'] = 'H(B)'
    with pytest.raises(ValueError):
        packet.validate(bad, value)
    bad = deepcopy(value); bad['physical_comparison_count_added'] = True
    with pytest.raises(ValueError):
        packet.validate(bad, value)


def test_retained_receipt_matches_cold_replay():
    from tools import strict_json
    from . import packet
    packet.validate(strict_json.load(packet.OUTPUT), packet.report())
