"""Stored-result validation without recomputing the scientific experiment."""
from __future__ import annotations

import copy
import json

import numpy as np
import pytest

from . import verify


def failures(expected, actual, path='', full=True):
    out = []
    verify.compare(expected, actual, path, out, full)
    return out


@pytest.mark.parametrize('path', [
    'theorem.level0.adjacent_block_largest_level',
    'collar.reading.sampled.rate_symmetric',
    'particle_stability.single_readings.no_pair_runs',
])
@pytest.mark.parametrize('bad', [float('nan'), float('inf'), -float('inf')])
@pytest.mark.parametrize('stored', [True, False])
def test_nonfinite_values_fail_every_numeric_comparison_branch(path, bad, stored):
    value = 786 if path in verify.SAMPLED else 0.5
    expected, actual = (bad, value) if stored else (value, bad)
    assert failures(expected, actual, path)


@pytest.mark.parametrize('path', [
    'carrier.cell.ports',
    'collar.reading.sampled.rate_symmetric',
    'particle_stability.single_readings.no_pair_runs',
])
@pytest.mark.parametrize('boolean', [False, True])
def test_boolean_cannot_replace_a_numerical_record(path, boolean):
    actual = float(boolean) if path.startswith('collar') else int(boolean)
    assert failures(boolean, actual, path)
    assert failures(actual, boolean, path)


@pytest.mark.parametrize('stored', [12.0, 12.000000001])
def test_integer_records_require_integers(stored):
    assert failures(stored, 12, 'carrier.cell.ports')


@pytest.mark.parametrize('path,tolerance', verify.SAMPLED.items())
def test_sampled_count_tolerance_is_retained(path, tolerance):
    assert not failures(100 + tolerance, 100, path)
    assert failures(100 + tolerance + 1, 100, path)
    assert failures(100.0, 100, path)


def test_complete_stored_result_replays_and_nonfinite_sample_is_rejected():
    actual = json.loads(verify.RESULTS.read_text(encoding='utf-8'))
    assert not failures(actual, actual)
    stored = copy.deepcopy(actual)
    stored['collar']['reading']['sampled']['rate_symmetric'] = float('nan')
    # Python's JSON parser admits NaN; comparison must reject it independently.
    decoded = json.loads(json.dumps(stored))
    assert failures(decoded, actual) == ['collar.reading.sampled.rate_symmetric: non-finite value']


def test_full_comparison_rejects_missing_and_unknown_scientific_outputs():
    actual = {'scientific': {'value': 0.5}, 'mode': 'default'}
    assert failures({'mode': 'full'}, actual)
    stored = dict(actual, unexpected={'value': 0.5})
    assert failures(stored, actual)
    assert not failures(stored, actual, full=False)


def test_float_tolerances_remain_unchanged():
    assert not failures(0.5 + 5e-9, 0.5, 'scientific.value')
    assert failures(0.5 + 2e-8, 0.5, 'scientific.value')
    assert not failures(0.5 + 5e-5, 0.5, 'collar.reading.sampled.rate_symmetric')
    assert failures(0.5 + 2e-4, 0.5, 'collar.reading.sampled.rate_symmetric')


def test_doubled_complement_excludes_negative_coupling_from_block_bound():
    from .carrier import cell
    from .sectors import MultisetSector

    sector = MultisetSector(cell(), 2)
    doubled = np.array([2 in occupation for occupation in sector.states])
    uniform = doubled.astype(float) / np.sqrt(doubled.sum())
    for kappa in (-1.0, 0.0, 1.0):
        action = sector.hamiltonian(kappa) @ uniform
        assert np.linalg.norm(action - 24 * kappa * uniform) < 1e-12
        energy = float(uniform @ action)
        if kappa < 0:
            assert energy < 12 * kappa
        else:
            assert energy >= 12 * kappa - 1e-12
