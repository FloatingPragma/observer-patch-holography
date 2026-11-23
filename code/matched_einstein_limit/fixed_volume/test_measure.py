import json
import math

import pytest
from scipy.integrate import quad

from . import measure


def test_three_channel_replay_and_refinement():
    receipt = measure.run()
    measure.verify_replay(receipt, json.loads((measure.HERE / 'measurements.json').read_text()))
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
