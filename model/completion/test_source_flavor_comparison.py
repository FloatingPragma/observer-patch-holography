"""Independent rational controls and exact-custody rejection for source ratios."""
import ast
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import shutil

import pytest

from tools import strict_json
from . import source_flavor_comparison as comparison


def test_exact_interval_division_and_zero_denominator():
    assert comparison.ratio_interval(['2', '3'], ['5', '7']) == ['2/7', '3/5']
    assert comparison.ratio_interval(['0', '0'], ['5', '7']) == ['0', '0']
    assert comparison.ratio_interval(['5', '7'], ['0', '0']) is None
    with pytest.raises(ValueError, match='denominator'):
        comparison.ratio_interval(['2', '3'], ['0', '7'])


@pytest.mark.parametrize('interval,target,expected', [
    (['1', '1'], '9', '4/5'),
    (['8', '10'], '2', '3/5'),
    (['2', '3'], '5/2', '0'),
    (['0', '0'], '1', None),
])
def test_relative_correction_bound_by_independent_scalar_examples(interval, target, expected):
    assert comparison.necessary_correction(interval, target) == expected


def test_correction_bound_is_invariant_under_common_positive_scaling():
    numerator, denominator = ['2', '3'], ['5', '7']
    scale = F(17, 13)
    scaled = [[str(F(x)*scale) for x in box] for box in (numerator, denominator)]
    assert comparison.ratio_interval(*scaled) == comparison.ratio_interval(numerator, denominator)
    # At delta=4/5, e_i=4/5 and e_j=-4/5 give the extremal ratio9.
    assert (1+F(4, 5))/(1-F(4, 5)) == 9


def test_all_source_policies_use_independently_divided_retained_intervals():
    result = comparison.build()
    source = strict_json.load(comparison.ROOT/comparison.SOURCE)
    benchmark = strict_json.load(comparison.ROOT/comparison.BENCHMARK)
    y = list(map(F, benchmark['matched_SM']['yukawa']['charged_lepton']))
    for row, original in zip(result['comparisons'], source['current_gated_joint_spectrum']['policies'], strict=True):
        boxes = next(x for x in original['channels'] if x['channel'] == 'charged_lepton')['dimensionless_singular_values']
        assert row['policy'] == original['policy']
        for name, i, j in comparison.RATIOS:
            item = row['ratios'][name]
            assert F(item['benchmark_exact']) == y[i]/y[j]
            if F(boxes[j][0]) == 0:
                assert item['source_interval_exact'] is None
                continue
            expected = [F(boxes[i][0])/F(boxes[j][1]), F(boxes[i][1])/F(boxes[j][0])]
            assert list(map(F, item['source_interval_exact'])) == expected
            assert not expected[0] <= y[i]/y[j] <= expected[1]
        assert row['direct_common_scale_identification_excluded'] is True
    always, current, event, joint = result['comparisons']
    assert F(always['necessary_uniform_delta_lower_bound_exact']) > F(79, 100)
    assert F(current['necessary_uniform_delta_lower_bound_exact']) > F(59, 100)
    for row in (event, joint):
        assert row['necessary_uniform_delta_lower_bound_exact'] is None
        assert row['finite_multiplicative_correction_possible'] is False
        assert row['exact_zero_source_mode_against_positive_benchmark'] is True


@pytest.mark.parametrize('damage', ['ratio', 'delta', 'zero', 'policy', 'benchmark',
                                   'source_pin', 'missing_pin', 'extra', 'discharge',
                                   'prediction_type', 'correction_domain', 'url'])
def test_replay_rejects_data_contract_and_custody_mutations(damage):
    packet = deepcopy(comparison.build())
    if damage == 'ratio':
        packet['comparisons'][0]['ratios']['lightest_over_heaviest']['source_interval_exact'][0] = '1/2'
    elif damage == 'delta':
        packet['comparisons'][0]['necessary_uniform_delta_lower_bound_exact'] = '0'
    elif damage == 'zero':
        packet['comparisons'][2]['finite_multiplicative_correction_possible'] = True
    elif damage == 'policy':
        packet['comparisons'].pop()
    elif damage == 'benchmark':
        packet['benchmark']['charged_lepton_yukawas_exact'][0] = '1'
    elif damage == 'source_pin':
        packet['source_bindings'][comparison.SOURCE] = '0'*64
    elif damage == 'missing_pin':
        packet['source_bindings'].pop(comparison.BENCHMARK)
    elif damage == 'extra':
        packet['physical_species_identification'] = True
    elif damage == 'discharge':
        packet['canonical_physical_premises_discharged'] = ['assumed physical matching']
    elif damage == 'prediction_type':
        packet['physical_predictions_added'] = False
    elif damage == 'correction_domain':
        packet['correction_model']['domain'] = 'three fitted independent family scales'
    else:
        packet['primary_source_metadata'][0]['url'] = 'https://example.invalid'
    with pytest.raises(ValueError, match='replay'):
        comparison.validate(packet)


@pytest.mark.parametrize('name', tuple(comparison.SCOPE))
def test_every_scope_boundary_is_bound_to_fresh_arithmetic(name):
    packet = comparison.build()
    packet['scope'][name] = not packet['scope'][name]
    with pytest.raises(ValueError, match='replay'):
        comparison.validate(packet)


def test_changed_input_bytes_invalidate_retained_comparison(tmp_path):
    for name in comparison.DIRECT:
        target = tmp_path/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(comparison.ROOT/name, target)
    packet = comparison.build(tmp_path)
    path = tmp_path/comparison.BENCHMARK
    path.write_text(path.read_text()+'\n')
    with pytest.raises(ValueError, match='replay'):
        comparison.validate(packet, tmp_path)


@pytest.mark.parametrize('value', [True, 0.1, 'nan', 'inf', '1/0'])
def test_nonrational_or_untyped_inputs_rejected(value):
    with pytest.raises(ValueError):
        comparison.rational(value)


def test_calculation_imports_no_scientific_producer():
    tree = ast.parse(Path(comparison.__file__).read_text())
    imports = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend(x.name for x in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.append(node.module)
    assert set(imports) <= {'argparse', 'copy', 'fractions', 'hashlib', 'json', 'pathlib', 'tools'}


def test_round_trip_and_duplicate_json_key_rejection(tmp_path):
    path = tmp_path/'comparison.json'
    comparison.write(path)
    assert comparison.check(path)['all_policies_conditionally_disjoint'] is True
    text = path.read_text()
    path.write_text(text.replace('"schema":', '"schema":"duplicate","schema":', 1))
    with pytest.raises(strict_json.DuplicateKeyError):
        comparison.check(path)


def test_retained_comparison_matches_exact_replay():
    assert comparison.check()['policies'] == list(comparison.POLICIES)
