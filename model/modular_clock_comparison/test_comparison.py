"""Independent spectra/logarithms, complete-source and invalid-contract tests."""
from fractions import Fraction as F
import copy
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest
from scipy.linalg import expm, logm

from . import comparison as m
from . import intervals as q
from ..field_thermal.thermal import ThermalField, car_restriction, restrict
from ..full_fock_field.field import FieldGeometry, FullFock
from ..reference_observation.observation import ReferenceObservation


def _mp(fraction):
    return mp.mpf(fraction.numerator)/fraction.denominator


@pytest.mark.parametrize('x', (F(-4), F(-3, 8), F(0), F(1, 16), F(1), F(4)))
def test_exact_series_enclosures_against_independent_high_precision(x):
    with mp.workdps(100):
        e = q.exp(x)
        assert _mp(e.lo) <= mp.exp(_mp(x)) <= _mp(e.hi)
        assert e.hi-e.lo < F(1, 10**24)
        if x > 0:
            l = q.log(x)
            assert _mp(l.lo) <= mp.log(_mp(x)) <= _mp(l.hi)


def test_interval_outward_arithmetic_and_domain_rejection():
    x = q.Interval(F(1, 3), F(2, 3))
    y = q.Interval(F(-3, 2), F(-1, 2))
    assert (x*y).lo <= -1 and (x*y).hi >= F(-1, 6)
    assert x.lo <= F(1, 3) and x.hi >= F(2, 3)
    assert q.sqrt_integer(17).lo**2 <= 17 <= q.sqrt_integer(17).hi**2
    for call in (lambda: q.Interval(2, 1), lambda: q.Interval.point(.1),
                 lambda: q.Interval(-1, 1).reciprocal(), lambda: q.log(0),
                 lambda: q.exp(5), lambda: q.sqrt_integer(True)):
        with pytest.raises(ValueError):
            call()


def test_exact_source_spectrum_and_rational_frequency_separation():
    certificate = m.frequency_certificate()
    assert certificate['frequency_pair_checks'] == 243
    assert F(certificate['minimum_frequency_separation_lower']) > F(1, 10000)
    assert certificate['regional_centralizer_complex_dimension'] == 6
    # Independent eigendecomposition and signed CAR reduction; no formula
    # from the producer is used to construct these field references.
    field = ThermalField().target()
    energies = np.linalg.eigvalsh(field['hamiltonian'])
    reduced = restrict(field['reference'], car_restriction(4, (0, 1)))
    diagonal = np.diag(reduced).real
    frequencies = [np.log(diagonal[a]/diagonal[b]) for a, b in ((0, 1), (1, 3), (0, 3))]
    for value, endpoints in zip(frequencies, certificate['positive_regional_frequency_intervals'], strict=True):
        lo, hi = map(F, endpoints)
        assert abs(value-float((lo+hi)/2)) < 3e-14
        assert hi-lo < F(1, 10**24)
    assert min(abs(local-(left-right)) for local in frequencies
               for left in energies for right in energies) > 1e-4
    # The small p0-p1 splitting is physical in this finite model and cannot
    # be discarded by an eigensolver's relative grouping tolerance.
    assert diagonal[0]-diagonal[1] > 3e-5


def test_spectral_attachment_rejects_mutated_source_action(monkeypatch):
    original = FullFock.hamiltonian
    def shifted(self, geometry):
        from scipy import sparse
        return original(self, geometry)+sparse.eye(self.dimension)/16
    monkeypatch.setattr(FullFock, 'hamiltonian', shifted)
    with pytest.raises(ArithmeticError, match='characteristic polynomial'):
        m.exact_spectral_attachment()


def test_local_modular_actions_use_actual_reduced_density_and_keep_exterior():
    service = m.ClockComparison()
    rng = np.random.default_rng(1148)
    local = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    outside = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    s = F(2, 5)
    rho = service.local_references[0]
    modular_left = expm(1j*float(s)*logm(rho))
    w = service.unitary('site0', s)
    embedded = np.kron(np.eye(4), local)
    expected = np.kron(np.eye(4), modular_left@local@modular_left.conj().T)
    np.testing.assert_allclose(w.conj().T@embedded@w, expected, atol=2e-13)
    np.testing.assert_allclose(w@np.kron(outside, np.eye(4)), np.kron(outside, np.eye(4))@w, atol=2e-13)
    # A changed regional action is an active perturbation of this correlated
    # global reference, even though the regional density itself is fixed.
    assert np.linalg.norm(w@service.reference@w.conj().T-service.reference, 'nuc') > .005
    np.testing.assert_allclose(restrict(w@service.reference@w.conj().T, car_restriction(4, (0, 1))), rho, atol=2e-13)


def test_time_cocycle_and_regional_composition_hold_on_arbitrary_reference_columns():
    service = m.ClockComparison()
    a, b = F(2, 5), F(-1, 7)
    w = service.unitary('global', a)
    c_a = service.comparison('site0', 'global', a)
    c_b = service.comparison('site0', 'global', b)
    np.testing.assert_allclose(c_a@w@c_b@w.conj().T,
                               service.comparison('site0', 'global', a+b), atol=2e-14)
    # A cocycle is generally not itself a one-parameter group.
    assert np.linalg.norm(c_a@c_b-service.comparison('site0', 'global', a+b), 2) > .001
    rng = np.random.default_rng(83)
    columns = rng.normal(size=(16, 3))+1j*rng.normal(size=(16, 3))
    columns /= np.linalg.norm(columns)
    chain = service.comparison('site1', 'site0', a)@c_a
    np.testing.assert_allclose(chain@columns,
                               service.comparison('site1', 'global', a)@columns, atol=2e-14)


def test_source_comparison_keeps_all_input_columns_exports_and_complex_references():
    service = m.ClockComparison()
    run, certificate = service.certificate()
    assert certificate['input_basis_columns'] == 16
    assert certificate['resources']['native_pulses'] == 2826
    assert certificate['first_export_failure_slots_retained'] == len(run.word)
    assert certificate['complete_dilation_operator_error'] < float(F(certificate['ideal_split_operator_bound']))
    np.testing.assert_allclose(run.success.conj().T@run.success+run.failure_gram,
                               np.eye(16), atol=4e-12)
    rng = np.random.default_rng(81)
    columns = rng.normal(size=(16, 3))+1j*rng.normal(size=(16, 3))
    columns /= np.linalg.norm(columns)
    second = m.execute_complete(service.comparison_word('site0', 'global', F(1, 2)), 4, columns)
    np.testing.assert_allclose(second.success, run.success@columns, atol=2e-13)
    np.testing.assert_allclose(second.failure_gram, columns.conj().T@run.failure_gram@columns, atol=2e-13)
    wrong = service.unitary('site0', F(-1, 2))@service.unitary('global', F(1, 2)).conj().T
    assert np.linalg.norm(run.success-wrong, 2) > .5


def test_twenty_pulse_interregional_comparison_and_no_hopping_covariance_control():
    service = m.ClockComparison()
    _, result = service.certificate('site1', 'site0')
    assert result['resources']['native_pulses'] == 20
    assert result['ideal_split_operator_bound'] == '0'
    # With zero boundary hopping, the product reference has a covariant
    # nontrivial Petz observation: this is outside the nonresonance target.
    product = m.ClockComparison(FieldGeometry((2, 1, 1), speed=0.))
    obs = ReferenceObservation(product.reference, car_restriction(4, (0, 1)))
    x = np.zeros((16, 16), complex); x[0, 3] = 1
    s = F(1, 2)
    source = product.unitary('global', s)
    local = expm(-1j*float(s)*logm(product.local_references[0]))
    np.testing.assert_allclose(obs.observe(source.conj().T@x@source),
                               local.conj().T@obs.observe(x)@local, atol=2e-13)
    assert np.linalg.norm(local.conj().T@obs.observe(x)@local-obs.observe(x)) > .01


def test_adversarial_source_leakage_is_not_conditioned_away(monkeypatch):
    service = m.ClockComparison()
    original = m.compiler.native_pulse
    rotation = np.eye(6, dtype=complex)
    theta = .001
    rotation[0, 0] = rotation[4, 4] = np.cos(theta)
    rotation[4, 0], rotation[0, 4] = np.sin(theta), -np.sin(theta)
    monkeypatch.setattr(m.compiler, 'native_pulse', lambda kind, angle: rotation@original(kind, angle))
    with pytest.raises(ArithmeticError, match='error budget'):
        service.certificate('site1', 'site0')


@pytest.mark.parametrize('label,s,steps', [('bad', F(1), 8), ('site0', .5, 8),
                                        ('site0', F(5), 8), ('site1', F(1), 0),
                                        ('site1', F(1), True)])
def test_invalid_clock_contracts_fail_before_source_execution(label, s, steps):
    with pytest.raises(ValueError):
        m.ClockComparison().word(label, s, steps)


def test_retained_exact_certificate_is_replayed():
    retained = json.loads(Path(__file__).with_name('results.json').read_text())
    assert retained['frequency_certificate'] == m.frequency_certificate()
    assert retained['canonical_physical_premises_discharged'] == 0
    assert not retained['physical_time_or_energy_calibration_established']


def test_cli_replays_the_whole_receipt_and_rejects_tampering(tmp_path, monkeypatch):
    from . import __main__ as cli
    current = m.report()
    monkeypatch.setattr(cli, 'report', lambda: current)
    receipt = tmp_path/'receipt.json'
    receipt.write_text(json.dumps(current))
    cli.main(['--check', str(receipt)])
    corrupted = copy.deepcopy(current)
    corrupted['source_clock_comparisons'][0]['resources']['native_pulses'] -= 1
    receipt.write_text(json.dumps(corrupted))
    with pytest.raises(ValueError, match='Receipt mismatch'):
        cli.main(['--check', str(receipt)])
    corrupted = copy.deepcopy(current)
    corrupted['frequency_certificate']['simple_strict_separation_lower'] = '0'
    receipt.write_text(json.dumps(corrupted))
    with pytest.raises(ValueError, match='Receipt mismatch'):
        cli.main(['--check', str(receipt)])
