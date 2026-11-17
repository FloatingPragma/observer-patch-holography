"""Scientific predicates exercised by valid identities and changed processes."""
from fractions import Fraction as F
import json
import math

import numpy as np
import pytest
import sympy as s

from . import bounds as iv
from . import check as j
from . import crosscheck


@pytest.fixture(scope='module')
def certificate():
    return j.result()


def test_exact_source_modes_and_bounded_scalar_rejection(certificate):
    residual = tuple(map(F, certificate['second_energy_residual_interval']))
    assert F(6515, 10000) < residual[0] < residual[1] < F(6516, 10000)
    assert certificate['positive_degenerate_mode_control']
    # The same certified disjointness predicate also accepts equal data.
    a = tuple(map(F, certificate['energy_intervals'][0]))
    assert not iv.disjoint(a, a)
    assert iv.disjoint(tuple(map(F, certificate['energy_intervals'][1])),
                       tuple(map(F, certificate['second_energy_prediction_interval'])))


def test_same_record_with_certified_one_word_difference(certificate):
    lo, hi = map(F, certificate['one_word_coordinate_probability_gap_interval'])
    assert 0 < lo < hi < F(1, 100000)
    assert all(row['witness_arc'] == 10243 for row in certificate['records'])
    assert all(row['theta_cosine_coefficient_display'] != 0 for row in certificate['records'])
    f = j.modes()[1]
    transfer = sum(map(F, certificate['transfer_intervals'][1]))/2
    u, v = j.embed(f), j.reference()
    su = u[j.source()[2]]
    w = (su-float(transfer)*u)/math.sqrt(1-float(transfer)**2)
    phi = float(sum(map(F, certificate['second_phase_prediction_interval']))/2)
    a = j.law(j.walk((u+v)/math.sqrt(2)))
    b = j.law((math.cos(phi)*u+math.sin(phi)*w+v)/math.sqrt(2))
    actual = a[1+10243]-b[1+10243]
    assert abs(actual-float((lo+hi)/2)) < 1e-20


def test_nonidentical_cell_mutation_fails_structural_admission(tmp_path):
    data = json.loads(j.GRAPH.read_text())
    data['seams'][0] = [0, 4]
    path = tmp_path/'graph.json'
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='^identical ordered internal cells required$'):
        j.source(path)


def test_intercell_mutation_passes_source_admission_and_fails_anchor_equation(tmp_path):
    data = json.loads(j.GRAPH.read_text())
    # Keep every internal cell intact; attach this external seam to a nonzero
    # anchor port instead. The simple graph and maximum degree six survive.
    index = data['seams'].index([1, 82])
    data['seams'][index] = [4, 82]
    path = tmp_path/'graph.json'
    path.write_text(json.dumps(data))
    j.source(path)
    with pytest.raises(ValueError, match='^actual compact anchor eigen-equation$'):
        j.modes(path)


def test_full_law_predicate_keeps_failures_and_coherent_reference(certificate):
    for row in certificate['records']:
        delta = row['negative_control_max_coordinate_differences_display']
        assert set(delta) >= {'double_frequency', 'reversed_phase', 'reference_quarter_turn',
                              'reference_dephased', 'reference_removed', 'failures_removed'}
        assert min(delta.values()) > 1e-14
        assert row['integer_word_2pi_alias_non_discriminant']
    assert 'anchored_log_phase_prediction' in certificate['records'][1]['negative_control_max_coordinate_differences_display']
    assert F(certificate['preparation_success_exact'])+F(certificate['preparation_fail_exact']) == 1


def test_one_scale_rule_detects_per_mode_recalibration(certificate):
    k = tuple(map(F, certificate['anchor_scale_interval']))
    e = tuple(map(F, certificate['energy_intervals'][1]))
    theta = tuple(map(F, certificate['phase_intervals'][1]))
    hidden = iv.div(e, theta)
    assert j.same_conversion(k, k)
    assert not j.same_conversion(k, hidden)


def test_reflection_is_part_of_evolution_word():
    f = j.modes()[0]
    u, v = j.embed(f), j.reference()
    psi = (u+v)/math.sqrt(2)
    correct = j.law(j.walk(j.walk(psi)))
    reverse = j.source()[2]
    changed = j.law(psi[reverse][reverse])
    assert not j.same_law(correct, changed)


def test_padding_is_reducing_and_excluded_from_preparation():
    for position, event in ((0, j.K), (j.V, 0)):
        state = np.zeros(j.D)
        state[position*j.J+event] = 1
        assert not j.admitted(state)
        with pytest.raises(ValueError):
            j.law(state)
        image = j.walk(state)
        assert np.array_equal(abs(image), state)
        assert image[position*j.J+event] == (-1 if position < j.V else 1)


def test_all_input_compression_positive_identity_and_boundaries():
    _, L, reverse, _ = j.source()
    assert np.array_equal(reverse[reverse], np.arange(j.D))
    # Integer Chebyshev recurrence on the actual source: K² A*W²A=2C²-K²I.
    C = j.K*s.eye(j.V)-L
    for profile, ell in zip(j.modes()[:2], (s.Integer(6), 5+s.sqrt(5))):
        assert j.zero((2*C*C-j.K**2*s.eye(j.V))*profile
                      -(2*(j.K-ell)**2-j.K**2)*profile)
    assert iv.acos_point(F(1)) == (0, 0)
    assert iv.minus_log_point(F(1)) == (0, 0)
    for bad in (F(-1), F(0), F(3, 2)):
        with pytest.raises(ValueError):
            iv.spectral((bad, bad))
    with pytest.raises(ValueError):
        iv.minus_log_point(F(209, 210), 1)
    constant = j.embed(s.ones(j.V, 1))
    stationary = (constant+j.reference())/math.sqrt(2)
    assert j.same_law(j.law(j.walk(stationary)), j.law(stationary))


def test_internal_opposition_uses_a_different_spectral_and_interval_method(certificate):
    alternate = crosscheck.oppose()
    for key in ('energy_intervals', 'phase_intervals'):
        for primary, other in zip(certificate[key], alternate[key]):
            a, b = tuple(map(F, primary)), tuple(map(F, other))
            assert not iv.disjoint(a, b)
    assert not alternate['independent_reviewer']
