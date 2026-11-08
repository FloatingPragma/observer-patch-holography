"""Independent covariance, dark-sector and source-leakage regression controls."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp

from . import bounds, localization as loc, packet


@pytest.fixture(scope='module')
def source_data():
    return loc.source_dictionary()


@pytest.fixture(scope='module')
def result():
    return packet.produce()


def test_full_character_table_is_orthogonal():
    c = np.array([[loc.character(k, x) for x in range(8)] for k in range(8)], dtype=object)
    assert np.array_equal(c@c.T, 8*np.eye(8, dtype=object))


@pytest.mark.parametrize('support', [(), (0,), (0, 3, 5), tuple(range(8))])
def test_missing_characters_give_complete_tight_frame_not_false_pvm(support):
    gram = np.array(loc.frame_gram(support), dtype=object)
    assert np.array_equal(gram@gram, gram)
    assert np.trace(gram) == len(support)
    if len(support) < 8:
        assert not np.array_equal(gram, np.eye(8, dtype=object))


@pytest.mark.parametrize('support', [(0, 0), (8,), (True,), [0, 1]])
def test_reject_invalid_character_support(support):
    with pytest.raises(ValueError):
        loc.frame_gram(support)


def test_actual_source_detector_full_covariance_and_dark_completeness(source_data):
    columns = source_data['columns']
    psi = np.column_stack([np.array(c, dtype=float)/np.sqrt(float(n))
                           for c, n in zip(columns, source_data['norms'], strict=True)])
    fourier = np.array([[loc.character(k, x)/np.sqrt(8) for x in range(8)] for k in range(8)])
    eta = psi@fourier
    projection = eta@eta.T
    dark = np.eye(256)-projection
    assert np.linalg.norm(eta.T@eta-np.eye(8)) < 1e-12
    assert np.linalg.norm(dark@dark-dark) < 1e-12
    assert abs(np.trace(dark)-248) < 1e-12
    for shift, permutation in enumerate(loc.algebra.reflections()):
        assert np.linalg.norm(eta[list(permutation)]-eta[:, [x ^ shift for x in range(8)]]) < 1e-12
        assert np.linalg.norm(dark[np.ix_(permutation, permutation)]-dark) < 1e-12


def test_literal_probe_orbit_is_not_orthogonal(source_data):
    weights = source_data['weights']
    orbit_gram = np.array([[sum(F(loc.character(k, x)*loc.character(k, y))*weights[k]
                               for k in range(8)) for y in range(8)] for x in range(8)], dtype=object)
    assert not np.array_equal(orbit_gram, np.eye(8, dtype=object))
    assert max(weights) > 10*min(weights)


def test_energy_leakage_independently_from_full_source_matrix(source_data):
    psi = np.column_stack([np.array(c, dtype=float)/np.sqrt(float(n))
                           for c, n in zip(source_data['columns'], source_data['norms'], strict=True)])
    h = np.array(source_data['h'], dtype=float)/source_data['scale']
    v = (np.eye(256)-psi@psi.T)@h@psi[:, 7]
    variance = float(np.dot(h@psi[:, 7], h@psi[:, 7])-np.dot(psi[:, 7], h@psi[:, 7])**2)
    assert np.dot(v, v) > 0.01
    assert abs(np.dot(v, v)-variance) < 1e-12


def test_source_probability_intervals_retain_all_outcomes(result):
    finite = result['finite_source_readout']
    intervals = [tuple(map(F, row['original_probability_interval'])) for row in finite['position_probabilities']]
    assert len(intervals) == 8
    assert sum(x[0] for x in intervals) <= 1 <= sum(x[1] for x in intervals)
    assert finite['actual_cyclic_probe_dark_probability'] == '0'
    assert finite['cyclic_probe_space_invariant_under_source_H'] is False


def test_discrete_uncertainty_has_zero_and_saturating_states():
    x = np.diag([0., 1.]).astype(complex)
    r = (np.eye(2)-np.array([[0, 1], [1, 0]]))/2
    eigen_x = np.array([1., 0.], dtype=complex)
    eigen_y = np.array([1., 1j])/np.sqrt(2)
    def variance(a, state):
        mean = np.vdot(state, a@state)
        return float(np.real(np.vdot(state, a@a@state)-mean*mean))
    assert variance(x, eigen_x) == 0
    assert abs(variance(x, eigen_y)*variance(r, eigen_y)-1/16) < 1e-14
    assert abs(np.vdot(eigen_y, (x@r-r@x)@eigen_y)) > 0.49


def test_fourier_coordinate_covariance_and_ccr_sign_on_schwartz_core():
    k, a, hbar = sp.symbols('k a hbar', real=True)
    f = (1+k)*sp.exp(-k*k)
    x = lambda g: sp.I*sp.diff(g, k)
    p = lambda g: hbar*k*g
    translated = sp.exp(-sp.I*a*k)*f
    assert sp.simplify(sp.exp(sp.I*a*k)*x(translated)-x(f)-a*f) == 0
    assert sp.simplify(x(p(f))-p(x(f))-sp.I*hbar*f) == 0


def test_gram_bounds_control_finite_compression():
    data = bounds.gram_bounds()
    assert data['weighted_field_gram_lower'] == F(1, 1000000)
    assert data['gram_condition_number_upper'] == 3000000
    actual = bounds.frame_perturbation(F(1, 1000), F(1, 2000000))
    assert actual['estimated_gram_lower'] == F(1, 2000000)


@pytest.mark.parametrize('value', [F(0), F(-1), F(1), 0.001, True])
def test_reject_unproved_order_lower_bound(value):
    with pytest.raises(ValueError):
        bounds.gram_bounds(value)


@pytest.mark.parametrize('error', [F(-1), F(1, 1000000), F(1), 1e-7])
def test_reject_gram_error_at_or_above_gap(error):
    with pytest.raises(ValueError):
        bounds.frame_perturbation(F(1, 1000), error)


def test_complete_receipt_roundtrip_and_scope_mutation(result):
    assert packet.validate(result, result)
    changed = deepcopy(result)
    changed['scope']['carrier_rank3_record_equals_continuum_position'] = True
    with pytest.raises(ArithmeticError):
        packet.validate(changed, result)


def test_supplied_receipt_path_is_read_and_malformed_receipt_fails(tmp_path, result):
    proposed = tmp_path/'supplied.json'
    proposed.write_text(packet.canonical(result))
    root = str(packet.ROOT)
    passed = subprocess.run([sys.executable, '-m', 'model.spin_covariant_localization', '--check', str(proposed)],
                            cwd=root, capture_output=True, text=True, check=False)
    assert passed.returncode == 0, passed.stderr
    proposed.write_text('{"schema":"fabricated"}')
    failed = subprocess.run([sys.executable, '-m', 'model.spin_covariant_localization', '--check', str(proposed)],
                            cwd=root, capture_output=True, text=True, check=False)
    assert failed.returncode != 0
    assert 'does not replay' in failed.stderr
