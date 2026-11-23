import json
import math

from scipy.integrate import quad

from . import measure

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
