"""Tests target nontrivial bridge boundaries and complete records."""
from fractions import Fraction as F
import numpy as np
import pytest
import sympy as s
from scipy.linalg import expm

from .construction import (yukawa_code, yukawa_pulse, source_pulse, fine_pointer,
                    reset_receipt, source, link_audit, coefficient_precision)
from .finite_compiler import decompose, product_word, prepare, matrix


def test_actual_native_H_and_all_gauge_charges_intertwine():
    code, h, support = yukawa_code()
    assert support == [0, 53, 64, 88, 169]
    assert code.shape == (320, 2) and h == s.Matrix([[0, -2], [-2, 1]])


@pytest.mark.parametrize('time', [F(0), F(1, 3), F(-2, 5), F(3)])
def test_source_executes_full_phase_and_both_columns_without_trotter(time):
    _, h, _ = yukawa_code()
    full, receipt = yukawa_pulse(time)
    expected = expm(-1j * float(time) * np.asarray(h.tolist(), complex))
    assert np.linalg.norm(full[:2, :2] - expected, 2) < 8e-14
    assert receipt['source_pulses'] == 1
    assert np.linalg.norm(full[2:, :2], 2) < 1e-13
    assert np.linalg.norm(sum(k.conj().T @ k for k in source.transfer_kraus(2)) - np.eye(6), 2) == 0


@pytest.mark.parametrize('law', ['C', 'E'])
def test_complete_fine_pointer_preserves_reference_and_all_bit_sector_coherence(law):
    v = np.array([1, 2j, 3, -4j, 2, 1j], complex) / np.sqrt(35)
    rho = np.outer(v, v.conj())
    branches = fine_pointer(rho, law)
    assert abs(sum(b['probability'] for b in branches) - 1) < 1e-13
    for bit in (0, 1):
        out = sum(b['output'] for b in branches if b['record']['bit'] == bit)
        projector = np.diag([int(j // 3 == bit) for j in range(6)])
        assert np.linalg.norm(out - projector @ rho @ projector, 2) < 1e-13
    assert any(b['record']['coin'] == 0 for b in branches)
    assert any(b['record']['coin'] == 1 for b in branches)


def test_reset_retains_failed_paths_and_actual_native_records():
    result = reset_receipt(2)
    assert result['certificate']['input_columns'] == 78
    assert all(x['attempts'] == len(x['fine_record_tape']) for x in result['replay_examples'])
    assert any(x['success'] for x in result['replay_examples'])
    assert any(not x['success'] for x in result['replay_examples'])
    assert F(result['certificate']['exact_max_failure']) <= F(result['certificate']['failure_upper'])


def test_all_compensated_CAR_columns_and_exact_generation_state():
    result = link_audit()
    assert result['operator_intertwining_columns'] == 15 * 2 ** 15
    assert result['resources']['native_pulses'] == 15
    assert result['exact_electric_readout']['mean_electric'] == '-18144/7225'
    assert result['prepared_state_error_norm_diagnostic'] < 2e-13


def test_complex_H_general_compiler_phase_padding_and_error_bound():
    h = s.zeros(4)
    h[0, 0] = 2
    h[0, 1], h[1, 0] = s.I / 3, -s.I / 3
    h[1, 2] = h[2, 1] = s.Rational(2, 5)
    terms = decompose(h)
    assert sum((c * matrix(p) for p, c in terms), s.zeros(4)) == h
    word, report = product_word(h, F(2, 7), 5)
    actual = source.execute(word, 2)
    expected = expm(-2j / 7 * np.asarray(h.tolist(), complex))
    assert np.linalg.norm(actual - expected, 2) < float(s.sympify(report['operator_error_upper']))
    assert actual.shape == (4, 4)  # All padding and full phase remain.


def test_preparation_emits_source_word_for_relative_and_global_complex_phase():
    v = np.array([1j, -1, 2j, 0], complex) / np.sqrt(6)
    word = prepare(v)
    assert np.linalg.norm(source.execute(word, 2, np.eye(4)[:, 0]) - v) < 1e-12


def test_real_reflection_cannot_make_required_complex_density():
    # Exact rational Cayley hop on the native gauge code: real input becomes
    # density with nonzero imaginary off-diagonal. Every real orthogonal
    # word on that same real code preserves real density instead.
    v = s.Matrix([s.Rational(15, 17), -8 * s.I / 17])
    assert (v * v.H)[0, 1] == 120 * s.I / 289


def test_exact_port_precision_bound_and_all_export_failure_outputs():
    result = coefficient_precision(24, port_nudges={0: 1})
    assert result['source_roundoff_residual_separate_diagnostic'] <= F(result['single_pulse_operator_error_upper'])
    assert len(result['complete_export_outputs']) == 5
    assert abs(result['complete_raw_output_mass'] - 1) < 1e-12
    assert all((1 << 24) % F(x).denominator == 0 for x in result['quantized_coefficients'])
    assert any(x['probability'] > 0 for x in result['complete_export_outputs'][1:])


def test_rejects_unsupported_or_inexact_domains():
    with pytest.raises(ValueError):
        yukawa_pulse(0.5)
    with pytest.raises(ValueError):
        source_pulse({'unsupported': 1})
    with pytest.raises(ValueError):
        decompose(s.Matrix([[0., 1], [1, 0]]))
    with pytest.raises(ValueError):
        prepare([1, 1])


def test_package_authenticates_sources_and_rejects_altered_exact_receipt():
    from .__main__ import verify_pins, compare
    assert verify_pins() > 10
    compare({'exact': '1/3', 'diagnostic': 1e-16}, {'exact': '1/3', 'diagnostic': 2e-16})
    with pytest.raises(ValueError):
        compare({'exact': '1/4'}, {'exact': '1/3'})
    with pytest.raises(ValueError):
        compare({'diagnostic': float('nan')}, {'diagnostic': 0.})
