"""Independent covariance constructions and shared-reference failure controls."""
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
import pytest

from . import comparison_budget as cb
from .certificate import certificate, continuum_fit
from .__main__ import same


def test_psd_extrema_are_attainable_with_identical_marginals():
    # Rank-one positive correlation attains the upper endpoint. A Gram
    # triangle attains zero with the same diagonal (4,9,16).
    low = np.array([[4, 1.5, -5.5], [1.5, 9, -10.5], [-5.5, -10.5, 16]])
    high = np.outer([2, 3, 4], [2, 3, 4])
    assert cb.variance_range([1, 1, 1], [2, 3, 4]) == (0, 81)
    assert np.linalg.eigvalsh(low).min() > -1e-12
    assert low.sum() == 0 and high.sum() == 81
    assert np.array_equal(low.diagonal(), high.diagonal())
    assert cb.variance_range([1, -1], [5, 2]) == (9, 49)


def test_independent_sets_keep_only_the_actual_shared_covariance():
    # Y0 independent of (Yref,T), whose standardized correlation is +/-1.
    coeff = [F(3, 2), F(-3, 2), F(1)]
    errors = [F(2), F(3), F(3)]
    low, high = cb.variance_range(coeff, errors, ['a', 'ref', 'ref'])
    assert (low, high) == (F(45, 4), F(261, 4))
    for rho, expected in [(1, low), (-1, high)]:
        cov = [[4, 0, 0], [0, 9, 9*rho], [0, 9*rho, 9]]
        value = sum(coeff[i]*cov[i][j]*coeff[j] for i in range(3) for j in range(3))
        assert value == expected


def test_shared_peak_exact_cancellation_rejects_naive_independence():
    # The exact map T=Yref leaves c0=1.5Y0-.5Yref. Treating T and c0
    # independently misses the subtraction of Yref in the actual map.
    correct = F(3, 2)**2*4 + F(-1, 2)**2*9
    naive = 9 + correct
    assert correct == F(45, 4) and naive == F(81, 4)
    assert correct != naive


@pytest.mark.parametrize('coeff,errors', [([], []), ([1], [1, 2]), ([1], [-1]), ([float('nan')], [1]), ([True], [1])])
def test_invalid_variance_inputs_rejected(coeff, errors):
    with pytest.raises(ValueError):
        cb.variance_range(coeff, errors)


def test_exact_fit_influence_matches_perturbation_of_owning_fit():
    data = json.loads((cb.HERE/'production.json').read_text())['chain_point']
    points = [dict(beta_G=v['bG'], size=4*v['L']/v['bG'], n_side=v['L'],
                   value=v['yscan']['y_max'], err=v['yscan']['y_max_err'])
              for v in data.values() if abs(v['yscan']['x']-.296) < 1e-9]
    weights = cb.intercept_weights(points)
    central = continuum_fit(points)['continuum']
    for index, weight in enumerate(weights):
        changed = deepcopy(points)
        changed[index]['value'] += .00001
        observed = (continuum_fit(changed)['continuum']-central)/.00001
        assert abs(observed-float(weight)) < 1e-8
    with pytest.raises(ValueError, match='rank-deficient'):
        cb.solve([[1, 1], [1, 1]], [1, 0])


def test_rational_sd_enclosures_and_numerical_derivative_scope():
    row = cb.encode_range((F(2), F(7)))
    lower, upper = map(F, row['sd_bounds_exact_GeV'])
    assert lower*lower <= 2 and upper*upper >= 7
    value, difference = cb.slope_at(lambda t: t**3, 2.)
    assert abs(value-12) < 1e-10 and difference < 1e-9
    with pytest.raises(ValueError, match='positive stable'):
        cb.slope_at(lambda t: -t, 2.)


def test_input_custody_rejects_changed_or_new_joint_data(tmp_path, monkeypatch):
    shutil.copytree(cb.HERE/'inputs', tmp_path/'inputs')
    monkeypatch.setattr(cb, 'HERE', tmp_path)
    assert not cb.input_custody()[0]['paired_peak_replicas_retained']
    path = tmp_path/'inputs/EW/data/oph_set_b8_L24.json'
    value = json.loads(path.read_text()); value['paired_peak_replicas'] = [[1, 2], [2, 3]]
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match='comparison input differs'):
        cb.input_custody()
    manifest_path = tmp_path/'inputs/input_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    for row in manifest['files']:
        if row['path'] == 'EW/data/oph_set_b8_L24.json':
            row['bytes'] = path.stat().st_size
            row['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='inspect retained joint data'):
        cb.input_custody()


def test_receipt_preserves_null_uncertainty_and_separate_domains():
    value = certificate()
    budget = value['comparison_budget']
    assert same(budget, json.loads((cb.HERE/'comparison_budget.json').read_text()))
    row = value['rows']['declared']
    assert row['T_c_err_lattice'] is None and row['T_c_err_total'] is None
    assert value['comparison_not_blind']['pull'] is None
    assert row['T_c_err_lattice_range'][0] < row['T_c_err_lattice_independence_control'] < row['T_c_err_lattice_range'][1]
    assert row['truncation_error'] == .01*row['T_c']
    assert budget['finite_design_sensitivity_GeV'] == value['finite_size_cutoff']['matched_size_cutoff']['T_max']['sensitivity_envelope']
    assert not budget['rows']['declared']['confidence_coverage_certified']
    assert not budget['rows']['declared']['nonlinear_sampling_distribution_identified']
    # The x=.36 temperature comes from a separate run set, so under set
    # independence its linearized variance is identified from the marginals.
    b = budget['rows']['one_doublet_B']['independent_run_sets']['variance_bounds_exact']
    assert b[0] == b[1]
    for name, digest in budget['input_sha256'].items():
        assert hashlib.sha256((cb.HERE/name).read_bytes()).hexdigest() == digest
