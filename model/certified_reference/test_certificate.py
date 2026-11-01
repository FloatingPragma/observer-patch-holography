"""Independent enclosures, exact-source controls, and rejected countercontrols."""
from dataclasses import replace
from fractions import Fraction as F
import math

import mpmath as mp
import numpy as np
import pytest

from ..observer.source_scalar import compiler
from ..plaquette_reference import preparation as p
from . import interval as i, source, certificate as c, continuation, design_check
from .packet import modular_bounds


@pytest.fixture(scope='module')
def proved():
    return c.certify_preparation()


def mp_fraction(value):
    return mp.mpf(value.numerator)/value.denominator


@pytest.mark.parametrize('value', [F(-4), F(-1, 7), F(), F(1, 10**50),
                                  F(1, 3), F(3), F(4)])
def test_trig_outward_enclosures_against_high_precision(value):
    with mp.workdps(100):
        for enclosure, fn in zip(i.sine_cosine(value), (mp.sin, mp.cos), strict=True):
            expected = fn(mp_fraction(value))
            assert mp.mpf(enclosure.lo)/i.SCALE <= expected <= mp.mpf(enclosure.hi)/i.SCALE
            assert enclosure.hi-enclosure.lo < 10**8


def test_integer_sqrt_and_pi_enclose_independently():
    with mp.workdps(100):
        for value in (F(), F(1, 2), F(5), F(10**70+3, 10**81+7)):
            bounds = i.sqrt_fraction(value)
            exact = mp.sqrt(mp_fraction(value))
            assert mp.mpf(bounds.lo)/i.SCALE <= exact <= mp.mpf(bounds.hi)/i.SCALE
        bounds = i.pi()
        assert mp.mpf(bounds.lo)/i.SCALE < mp.pi < mp.mpf(bounds.hi)/i.SCALE
    assert 0 < i.pi_angle_error(F(math.pi), F(1)) < F(2, 10**16)


@pytest.mark.parametrize('bad', [F(5), F(-5), 0.1, True])
def test_out_of_domain_trigonometric_inputs_rejected(bad):
    with pytest.raises(ValueError):
        i.sine_cosine(bad)


def test_dyadic_interval_product_contains_signed_rational_corners():
    a = i.Interval.point(F(-7, 11))+i.Interval.point(F(1, 13))
    b = i.Interval(-2*i.SCALE, 3*i.SCALE)
    product = a*b
    for x in (F(a.lo, i.SCALE), F(a.hi, i.SCALE)):
        for y in (-2, 3):
            assert F(product.lo, i.SCALE) <= x*y <= F(product.hi, i.SCALE)


def test_independent_inverse_walsh_is_complete():
    values = [F((-1)**j*j, 31) for j in range(16)]
    explicit = [sum((v*(-1)**((p&m).bit_count()) for m, v in enumerate(values)), F())
                for p in range(16)]
    assert c.inverse_walsh(values) == explicit


def test_signed_final_leaf_is_not_replaced_by_populations():
    correct = c.two_vector_error(F(-math.pi/4), 1, 1, True)
    wrong = c.two_vector_error(F(-math.pi/4), 1, -1, True)
    assert correct < F(1, 10**15)
    assert wrong > 1


def test_actual_source_control_error_encloses_independent_high_precision():
    spec = source.specification()
    with mp.workdps(90):
        for kind, angle in [('H', math.pi/2), ('X', math.pi/4), ('Y', -.3),
                            ('Z', .001), ('I', 1.), ('J', -math.pi)]:
            controls = source.literal_controls(kind, angle)
            actual = mp.matrix(6)
            for row in range(6):
                for col in range(6):
                    for value, matrix in zip(controls, spec, strict=True):
                        real, imag = matrix[row][col]
                        actual[row, col] += mp_fraction(value)*(
                            mp_fraction(real[0])+mp.sqrt(5)*mp_fraction(real[1])+
                            1j*(mp_fraction(imag[0])+mp.sqrt(5)*mp_fraction(imag[1])))
            target = mp.matrix(6); a = mp_fraction(F(angle))
            if kind == 'X': target[0, 1] = target[1, 0] = 1j*a
            if kind == 'Y': target[0, 1], target[1, 0] = a, -a
            if kind in ('Z', 'I'):
                target[0, 0], target[1, 1] = 1j*a, (1j*a if kind == 'I' else -1j*a)
            if kind == 'J':
                for j in range(3): target[j, j] = 1j*a
            if kind == 'H':
                target[0, 0] = target[0, 1] = target[1, 0] = 1j*a/mp.sqrt(2)
                target[1, 1] = -1j*a/mp.sqrt(2)
            exact = max(sum(abs(actual[r, s]-target[r, s]) for s in range(6)) for r in range(6))
            bound = source.generator_error(kind, F(angle), controls)
            assert exact <= mp_fraction(bound)+mp.mpf('1e-85')


def test_wrong_port_control_has_positive_full_source_error():
    controls = list(source.literal_controls('H', math.pi/2))
    controls[0] += F(1, 1000)
    error = source.generator_error('H', F(math.pi/2), tuple(controls))
    assert error > F(1, 10000)


def test_complete_export_cannot_omit_zero_weight_complement(monkeypatch):
    source.complete_export_certificate()
    original = compiler.transfer_kraus
    monkeypatch.setattr(compiler, 'transfer_kraus', lambda *a, **kw: original(*a, **kw)[:-1])
    with pytest.raises(ArithmeticError, match='branch was dropped'):
        source.complete_export_certificate()


def test_gate_routing_corruption_cannot_hide_in_scalar_angles():
    word = p.probability_word()
    word[0] = {**word[0], 'targets': [0]}
    with pytest.raises(ArithmeticError, match='gate routing'):
        c.amplitude_certificate(word, p.design().column_norm_squared, probability=True)


def test_actual_tree_coefficient_corruption_is_detected():
    word = p.probability_word()
    # At the first level the only Y factor is X(-pi/4), Z(theta), X(pi/4).
    word[1] = {**word[1], 'angle': word[1]['angle']+.01}
    result = c.amplitude_certificate(word, p.design().column_norm_squared, probability=True)
    assert result['tree_error'] > F(1, 1000)


def test_missing_blanking_phase_and_inverse_frame_rejected():
    word = p.conditional_word(85)
    del word[4]
    with pytest.raises(ArithmeticError, match='gate routing'):
        c.consume(word, 0, c.blank_template(85))
    word = p.inverse_frame_word(); word.pop()
    with pytest.raises(ArithmeticError, match='ended inside'):
        c.consume(word, 0, c.frame_template())


def test_independent_full_action_horner_and_wrong_polynomial():
    assert design_check.certify()['all_column_probabilities_exactly_checked'] == 256
    data = p.design(); value = [list(row) for row in data.numerator]; value[0][0] += 1
    with pytest.raises(ArithmeticError, match='Horner polynomial'):
        design_check.certify(replace(data, numerator=tuple(map(tuple, value))))


def test_all_record_trajectory_bounds_and_modular_consequence(proved):
    assert len(proved['conditional_records']) == 256
    assert proved['all_conditional_tree_nodes'] == 65535
    assert proved['unique_literal_pulse_controls'] == 7612
    assert proved['complete_preparation_output_error'] < F(17, 10**13)
    assert proved['complete_preparation_output_error'] < p.EXECUTION_ADMISSION
    assert proved['complete_preparation_output_error'] == (
        proved['probability_complete_output_error']+
        max(row['complete_conditional_output_error'] for row in proved['conditional_records']))
    continued = continuation.certify()
    bounds = modular_bounds(proved['complete_preparation_output_error'],
                            continued['full_action_complete_operator_error'])
    assert bounds['regional_P8_total_operator_and_trace_zero_tangent_bound'] < F(2, 10**9)
    assert bounds['all_abort_probability_upper'] < F(3, 10**24)
    assert bounds['actual_regional_success_minimum_weight'] > F(1, 20)


def test_continued_field_word_all_input_operator_bound():
    result = continuation.certify()
    assert result['native_pulses'] == 4320
    assert result['literal_source_complete_operator_error'] < F(538, 10**15)
    assert result['ideal_Strang_error'] == F(1463, 805306368)
    assert result['all_Fock_inputs_and_untouched_references']


def test_packet_path_check_rejects_malformed_receipt(tmp_path, monkeypatch):
    # Exercise the actual CLI argument/check path without repeating the
    # independently tested full certificate merely to reject a bad payload.
    import sys
    from . import packet
    bad = tmp_path/'bad.json'; bad.write_text('{"all_pass":true}')
    monkeypatch.setattr(packet, 'run', lambda: {'all_pass': True, 'required_exact_bound': '1/10'})
    monkeypatch.setattr(sys, 'argv', ['certificate', '--check', str(bad)])
    with pytest.raises((ValueError, KeyError, AssertionError)):
        packet.main()
