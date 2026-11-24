import json
import math

import pytest
from scipy.integrate import quad

from . import measure


def test_boundary_chord_reconstruction_and_launch_weights():
    import numpy as np
    from . import chords

    replay = chords.run()
    saved = json.loads((measure.HERE / "chords.json").read_text())
    for panel, archive in zip(replay["panels"], saved["panels"], strict=True):
        for case, original in zip(panel["cases"], archive["cases"], strict=True):
            assert case["chord_round_trip_ticks"] == original["chord_round_trip_ticks"]
            assert abs(case["boundary_volume"] - original["boundary_volume"]) < 1e-13
        assert panel["maximum_clock_ratio_error_per_radius"] < 1e-7
    errors = [row["maximum_fractional_difference"] for row in replay["panels"]]
    assert all(a > b for a, b in zip(errors, errors[1:]))
    assert errors[-1] < 4.3e-8
    nodes, weights = np.polynomial.legendre.leggauss(16)
    mu, weights = (1 + nodes) / 2, weights / 2
    radius, tick = 0.5, 2.0**-42
    expected, area = measure.spatial_measures(0, radius)
    tapes = chords.clock_tapes(0., radius, tick, mu)
    assert chords.boundary_volume(area, tick, tapes, mu, weights)[0] == pytest.approx(expected)
    assert area / 4 * (weights @ tapes) * tick / 2 == pytest.approx(0.75 * expected)
    central = chords.clock_tapes(0., radius, tick, np.ones_like(mu))
    assert chords.boundary_volume(area, tick, central, mu, weights)[0] == pytest.approx(1.5 * expected)
    with pytest.raises(ValueError, match="strictly convex"):
        chords.clock_tapes(1., 2., tick, mu)

# Poisson draws of the record. NumPy samples them (and the arrival Gamma) by
# rejection, and the accept decisions follow the platform's libm: one platform
# replays the stored record bit for bit, another can take a different stream
# path with the same law.
DRAWS = ('held_interior_count', 'surface_count')


def _without_draws(value):
    if isinstance(value, dict):
        return {k: _without_draws(v) for k, v in value.items() if k not in DRAWS and k != 'inferred'}
    if isinstance(value, list):
        return [_without_draws(v) for v in value]
    return value


def _draws(value, path='$'):
    if isinstance(value, dict):
        for key, item in value.items():
            if key in DRAWS:
                yield f'{path}.{key}', item
            else:
                yield from _draws(item, f'{path}.{key}')
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from _draws(item, f'{path}[{index}]')


def test_three_channel_replay_and_refinement():
    receipt = measure.run()
    stored = json.loads((measure.HERE / 'measurements.json').read_text())
    measure.verify_replay(receipt, stored)
    # Everything that is no draw replays exactly, the stored readouts follow
    # from the stored draws, and each fresh draw lies within six standard
    # deviations of the stored one.
    assert _without_draws(receipt) == _without_draws(stored)
    for panel in stored['panels']:
        for case in panel['cases']:
            assert case['inferred'] == measure.interpret(panel['flat_radius'], panel['count_budget'],
                                                         panel['proper_clock_tick'], case['record'],
                                                         panel['flat_comparison'])
    for (path, fresh), (_, kept) in zip(_draws(receipt), _draws(stored), strict=True):
        assert abs(fresh - kept) <= 6 * math.sqrt(2 * kept), path
    maximum_residuals = []
    for panel in receipt['panels']:
        count = panel['flat_comparison']['regions']['volume_threshold']['held_interior_count']
        assert abs(count - panel['count_budget']) < 6 * math.sqrt(2 * panel['count_budget'])
        comparisons = []
        for case in panel['cases']:
            read = case['inferred']
            matched = read['volume_threshold']
            assert abs(matched['held_volume_count_standardized']) < 6
            comparisons.append(abs(matched['area_einstein_contraction'] - read['radar_einstein_contraction']))
        maximum_residuals.append(max(comparisons))
    assert all(a > b for a, b in zip(maximum_residuals, maximum_residuals[1:]))
    assert maximum_residuals[-1] < 0.015


def test_equal_radius_leading_coefficient_misses_equal_volume():
    fine = measure.run()['panels'][-1]
    for case in fine['cases']:
        truth = 3 * case['prepared_sectional_curvature']
        matched = case['inferred']['volume_threshold']['area_einstein_contraction']
        wrong_region = case['inferred']['fixed_radius']['area_einstein_contraction']
        assert abs(matched - truth) < 0.006
        if truth:
            assert abs(wrong_region / truth - 5 / 3) < 0.025
            assert abs(wrong_region - truth) > 0.5


def test_same_total_count_mutation_and_changed_decoder_fail_even_when_gates_pass():
    raw = (measure.HERE/'measurements.json').read_bytes()
    archived, altered = json.loads(raw), json.loads(raw)
    case = altered['panels'][-1]['cases'][0]
    regions = case['record']['regions']
    regions['volume_threshold']['surface_count'] += 1
    regions['fixed_radius']['surface_count'] -= 1
    with pytest.raises(ValueError, match='count stream'):
        measure.verify_replay(altered, archived)
    modified = json.loads(raw)
    modified['panels'][0]['cases'][0]['inferred']['volume_threshold']['area_contrast'] += 1e-12
    with pytest.raises(ValueError, match='decoded'):
        measure.verify_replay(modified, archived)


def test_volume_jacobian_and_independence_of_decoder_channels():
    for curvature in measure.CURVATURES:
        radius = 0.4
        volume, _ = measure.spatial_measures(curvature, radius)
        integrated_area = quad(lambda r: measure.spatial_measures(curvature, r)[1], 0, radius)[0]
        assert abs(volume - integrated_area) < 1e-12
    panel = measure.run()['panels'][0]
    case = panel['cases'][-1]
    record = case['record']
    reference = panel['flat_comparison']
    tick, ell, budget = panel['proper_clock_tick'], panel['flat_radius'], panel['count_budget']
    geometry = dict(record, triangle=reference['triangle'])
    comparison = measure.interpret(ell, budget, tick, geometry, reference)
    assert comparison['radar_einstein_contraction'] != case['inferred']['radar_einstein_contraction']
    assert comparison['volume_threshold'] == case['inferred']['volume_threshold']
