"""Independent logarithm, complex-reference and complete-output controls."""
from fractions import Fraction as F
import json
import copy
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm, logm

from . import flow as m
from ..full_fock_field.field import FieldGeometry
from ..field_thermal.thermal import car_restriction, restrict
from ..observer.source_scalar import compiler
from ..reference_observation.observation import ReferenceObservation


@pytest.mark.parametrize('axis', (0, 1, 2))
def test_complete_input_source_conjugation_is_the_full_density_modular_action(axis):
    geometry = FieldGeometry(tuple(2 if i == axis else 1 for i in range(3)))
    service = m.ModularFlow(geometry, F(3, 2))
    s = F(1, 3)
    run, cert = service.certificate(s, 8)
    rho = service.thermal.target()['reference']
    modular_left = expm(1j*float(s)*logm(rho))
    # The source implements the Schrödinger unitary; its adjoint is the
    # modular conjugating matrix, up to an irrelevant scalar phase.
    rng = np.random.default_rng(501+axis)
    x = rng.normal(size=(16, 16))+1j*rng.normal(size=(16, 16))
    x /= np.linalg.norm(x, 2)
    actual = run.success.conj().T@x@run.success
    expected = modular_left@x@modular_left.conj().T
    bound = float(F(cert['ideal_split_operator_bound']))
    assert np.linalg.norm(actual-expected, 2) <= 2*bound+2e-9
    wrong_sign = modular_left.conj().T@x@modular_left
    assert np.linalg.norm(expected-wrong_sign, 2) > .05
    np.testing.assert_allclose(run.success.conj().T@run.success+run.failure_gram,
                               np.eye(16), atol=2e-10, rtol=0)
    assert cert['Hamiltonian_model_duration'] == '-1/2'


def test_complex_reference_columns_are_covered_by_the_complete_input_certificate():
    service = m.ModularFlow()
    s, steps = F(-1, 5), 4
    operator_run, cert = service.certificate(s, steps)
    rng = np.random.default_rng(80)
    columns = rng.normal(size=(16, 3))+1j*rng.normal(size=(16, 3))
    columns /= np.linalg.norm(columns)
    run = m.execute_complete(service.word(s, steps), 4, columns)
    np.testing.assert_allclose(run.success, operator_run.success@columns, atol=2e-13)
    np.testing.assert_allclose(run.failure_gram, columns.conj().T@operator_run.failure_gram@columns,
                               atol=2e-13)
    np.testing.assert_allclose(run.success.conj().T@run.success+run.failure_gram,
                               columns.conj().T@columns, atol=2e-12)
    error = np.sqrt(np.linalg.norm(run.success-service.unitary(s)@columns)**2+run.failure_probability)
    assert error <= cert['complete_dilation_operator_error']+2e-12


def test_two_complete_reads_keep_all_four_noncommuting_histories():
    service = m.ModularFlow()
    segments = ((F(1, 3), 8), (F(-1, 4), 8))
    run, cert = service.history(segments, read_mode=0)
    # Direct matrix branching from log(rho), independent of the producer's
    # Hamiltonian-time conversion and target assembly.
    log_rho = logm(service.thermal.target()['reference'])
    phases = [expm(-1j*float(s)*log_rho) for s, _ in segments]
    p = [np.diag([int((i & 1) == b) for i in range(16)]) for b in (0, 1)]
    actual = run.success.reshape(4, 16, 16)
    for second in (0, 1):
        for first in (0, 1):
            expected = p[second]@phases[1]@p[first]@phases[0]
            # Density channels ignore the scalar Gibbs partition phase.
            x = np.eye(16)/16
            np.testing.assert_allclose(actual[first+2*second]@x@actual[first+2*second].conj().T,
                                       expected@x@expected.conj().T,
                                       atol=2*float(F(cert['ideal_complete_history_bound'])), rtol=0)
    np.testing.assert_allclose(run.success.conj().T@run.success+run.failure_gram,
                               np.eye(16), atol=3e-12)
    assert cert['history_outcomes'] == 4
    assert cert['initialized_record_bits'] == 2
    # The initial number read does not commute with hopping; changed read
    # outcomes are retained instead of postselecting equal records.
    assert np.linalg.norm(actual[1]) > .001
    assert np.linalg.norm(actual[2]) > .001


def test_modular_group_and_inverse_use_the_same_complete_action():
    service = m.ModularFlow()
    a, b = F(1, 7), F(-2, 9)
    ra, ca = service.certificate(a)
    rb, cb = service.certificate(b)
    rc, cc = service.certificate(a+b)
    bound = sum(float(F(c['ideal_split_operator_bound'])) for c in (ca, cb, cc))
    assert np.linalg.norm(rb.success@ra.success-rc.success, 2) < bound+2e-9
    inverse, _ = service.certificate(-a)
    np.testing.assert_allclose(inverse.success@ra.success, np.eye(16), atol=2e-11)


def test_prepared_reference_purifier_read_and_all_abort_slots_remain_present():
    service = m.ModularFlow()
    run, cert = service.prepared_readout(F(1, 4), 4, read_mode=0)
    rho = service.thermal.target()['reference']
    actual = run.success[:, 0].reshape(2, 16, 16)
    for b in (0, 1):
        p = np.diag([int((i & 1) == b) for i in range(16)])
        assert np.linalg.norm(actual[b].T@actual[b].conj()-p@rho@p, 'nuc') < 2*float(F(cert['ideal_split_operator_bound']))+2e-9
    assert abs(cert['complete_probability']-1) < 2e-11
    assert cert['nonselective_read_reference_disturbance'] > .01
    assert cert['initialized_work_bits'] == 9
    assert cert['retained_purifier_bits'] == 4
    assert len(run.word) == cert['resources']['native_pulses']
    assert len(run.failure_weights) == len(run.word)


def test_tracial_and_noncommuting_controls_are_distinct():
    zero = m.ModularFlow(FieldGeometry((2, 1, 1), speed=0., interaction=0.))
    run, cert = zero.certificate(F(1, 2), 1)
    np.testing.assert_allclose(zero.thermal.target()['reference'], np.eye(16)/16, atol=2e-15)
    np.testing.assert_allclose(run.success, np.eye(16), atol=2e-12)
    assert cert['ideal_split_operator_bound'] == '0'
    ordinary = m.ModularFlow()
    number = np.diag([i & 1 for i in range(16)])
    unitary = ordinary.unitary(F(1))
    assert np.linalg.norm(unitary.conj().T@number@unitary-number, 2) > .05


def test_global_modular_implementation_does_not_claim_regional_covariance():
    service = m.ModularFlow()
    rho = service.thermal.target()['reference']
    restriction = car_restriction(4, (0, 1))
    regional = restrict(rho, restriction)
    observe = ReferenceObservation(rho, restriction)
    x = np.diag([i & 1 for i in range(16)])
    left = expm(.7j*logm(rho)); local = expm(.7j*logm(regional))
    first = observe.observe(left@x@left.conj().T)
    second = local@observe.observe(x)@local.conj().T
    assert np.linalg.norm(first-second) > 1e-4


def test_wrong_source_sign_is_rejected_on_complete_basis(monkeypatch):
    service = m.ModularFlow()
    original = service.word
    monkeypatch.setattr(m.ModularFlow, 'word', lambda self, s, n=4: original(-s, n))
    with pytest.raises(ArithmeticError, match='full-algebra'):
        service.certificate(F(1, 4), 8)


def test_export_leakage_cannot_be_hidden_by_success_normalization(monkeypatch):
    original = compiler.native_pulse
    def leaking(kind, angle):
        full = original(kind, angle)
        d = 4 if kind == 'J' else 2
        rotation = np.eye(6)
        rotation[0, 0] = rotation[d, d] = np.cos(.02)
        rotation[0, d] = -np.sin(.02); rotation[d, 0] = np.sin(.02)
        return rotation@full
    monkeypatch.setattr(compiler, 'native_pulse', leaking)
    with pytest.raises(ArithmeticError, match='full-algebra'):
        m.ModularFlow().certificate(F(1, 4), 4)


@pytest.mark.parametrize('beta', (True, .5, F(0), F(-1), F(33)))
def test_invalid_temperature_domain_is_rejected(beta):
    with pytest.raises(ValueError):
        m.ModularFlow(beta=beta)


@pytest.mark.parametrize('parameter,steps', ((.25, 4), (True, 4), (F(5), 4), (F(1), 0), (F(1), True)))
def test_invalid_schedule_domain_is_rejected(parameter, steps):
    with pytest.raises(ValueError):
        m.ModularFlow().certificate(parameter, steps)


@pytest.mark.parametrize('segments,mode', (([], 0), ([(F(1),)], 0), ([(F(1), 4)], True), ([(F(1), 4)], 4)))
def test_invalid_history_domain_is_rejected(segments, mode):
    with pytest.raises(ValueError):
        m.ModularFlow().history(segments, read_mode=mode)


def test_retained_report_replays_with_unchanged_domains_and_resources():
    retained = json.loads(Path(__file__).with_name('results.json').read_text())
    current = m.report()
    m.validate_report_replay(retained, current)


def test_resource_certificate_accepts_equivalent_phase_costs_without_hiding_accounting():
    retained = json.loads(Path(__file__).with_name('results.json').read_text())
    alternate = copy.deepcopy(retained)
    read = alternate['prepared_reference_readout']
    # Observed Linux/Mac difference; both actual costs remain reported and
    # independently satisfy the compiler-derived bound and additive ledger.
    difference = 21673.76874403813-21675.17375156573
    read['resources']['integrated_processor_norm_upper'] += difference
    read['preparation_resource_certificate']['measured_preparation_processor_norm_upper'] += difference
    m.validate_report_replay(retained, alternate)
    read['resources']['integrated_processor_norm_upper'] += 1
    with pytest.raises(ValueError, match='accounting'):
        m.validate_report_replay(retained, alternate)


@pytest.mark.parametrize('value', (-1e-300, float('nan'), float('inf'), -float('inf'), 1.01))
def test_resource_replay_rejects_invalid_probability_before_tolerance(value):
    retained = json.loads(Path(__file__).with_name('results.json').read_text())
    changed = copy.deepcopy(retained)
    changed['prepared_reference_readout']['first_export_abort_probability'] = value
    with pytest.raises(ValueError, match='nonfinite|invalid replay probability'):
        m.validate_report_replay(changed, retained)


@pytest.mark.parametrize('value', (float('inf'), -float('inf')))
def test_resource_replay_rejects_matching_nonfinite_error_diagnostics(value):
    changed = json.loads(Path(__file__).with_name('results.json').read_text())
    changed['full_algebra_certificate']['complete_dilation_operator_error'] = value
    with pytest.raises(ValueError, match='nonfinite'):
        m.validate_report_replay(changed, copy.deepcopy(changed))


@pytest.mark.parametrize('path,value', ((('full_Fock_dimension',), 16.),
    (('canonical_physical_premises_discharged',), False),
    (('prepared_reference_readout', 'retained_pointer_bits'), True),
    (('prepared_reference_readout', 'resources', 'native_pulses'), 11660.)))
def test_resource_replay_preserves_strict_integer_and_boolean_types(path, value):
    retained = json.loads(Path(__file__).with_name('results.json').read_text())
    changed = copy.deepcopy(retained)
    parent = changed
    for key in path[:-1]:
        parent = parent[key]
    parent[path[-1]] = value
    with pytest.raises(ValueError, match='type mismatch'):
        m.validate_report_replay(changed, retained)


@pytest.mark.parametrize('damage', ('pulse_count', 'understated_error', 'over_budget', 'under_cost'))
def test_portable_resource_replay_rejects_regressions(damage):
    retained = json.loads(Path(__file__).with_name('results.json').read_text())
    changed = copy.deepcopy(retained)
    read = changed['prepared_reference_readout']
    if damage == 'pulse_count':
        read['resources']['native_pulses'] += 1
    elif damage == 'understated_error':
        read['complete_preparation_and_readout_error'] /= 2
    elif damage == 'over_budget':
        certificate = read['preparation_resource_certificate']
        excess = float(F(certificate['processor_norm_upper_rational']))+1-certificate['measured_preparation_processor_norm_upper']
        certificate['measured_preparation_processor_norm_upper'] += excess
        read['resources']['integrated_processor_norm_upper'] += excess
    else:
        certificate = read['preparation_resource_certificate']
        read['resources']['integrated_processor_norm_upper'] -= certificate['measured_preparation_processor_norm_upper']
        certificate['measured_preparation_processor_norm_upper'] = 0.
    with pytest.raises(ValueError):
        m.validate_report_replay(retained, changed)


def test_actual_resource_count_is_recomputed_from_executed_immutable_tape(monkeypatch):
    original = compiler.resources
    def undercount(word):
        result = original(word)
        result['integrated_processor_norm_upper'] -= 1
        return result
    monkeypatch.setattr(compiler, 'resources', undercount)
    with pytest.raises(ArithmeticError, match='executed source tape'):
        m.ModularFlow().prepared_readout()


def test_extra_identity_pulses_cannot_hide_in_portable_preparation_budget(monkeypatch):
    original = m.ThermalField.word
    def padded(self):
        return original(self)+[compiler.pulse('I', 2*np.pi, [0])]
    monkeypatch.setattr(m.ThermalField, 'word', padded)
    with pytest.raises(ArithmeticError, match='schedule'):
        m.ModularFlow().prepared_readout()


@pytest.mark.parametrize('bits', (4, 8))
def test_preparation_budget_covers_arbitrary_nonnegative_amplitudes_and_phases(bits):
    rng = np.random.default_rng(bits+709)
    amplitudes = rng.random(1 << bits)
    amplitudes /= np.linalg.norm(amplitudes)
    phases = rng.uniform(-np.pi, np.pi, 1 << bits)
    word = compiler.compile_real_state(amplitudes)+compiler.compile_diagonal(phases)
    result = m.preparation_resources(word, bits)
    assert result['measured_preparation_processor_norm_upper'] == pytest.approx(
        compiler.resources(word)['integrated_processor_norm_upper'], abs=1e-12)
    assert result['measured_preparation_processor_norm_upper'] <= float(F(result['processor_norm_upper_rational']))


def test_actual_tiny_amplitude_phase_variants_keep_outputs_and_audited_costs(monkeypatch):
    retained = json.loads(Path(__file__).with_name('results.json').read_text())
    def alternative_word(self):
        target = self.target()['purification']
        phases = np.angle(target)
        small = abs(target) < 1e-15
        phases[small] = np.random.default_rng(702).uniform(-np.pi, np.pi, np.count_nonzero(small))
        return compiler.compile_real_state(abs(target))+compiler.compile_diagonal(phases)
    monkeypatch.setattr(m.ThermalField, 'word', alternative_word)
    current = m.report()
    difference = (current['prepared_reference_readout']['resources']['integrated_processor_norm_upper']-
                  retained['prepared_reference_readout']['resources']['integrated_processor_norm_upper'])
    assert abs(difference) > .1
    m.validate_report_replay(retained, current)
