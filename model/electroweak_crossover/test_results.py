"""The stored results equal their recomputation from the production record, and the named checks pass."""
import json
from pathlib import Path

from .__main__ import same
from .certificate import REFERENCE, X_C, certificate

RESULT = certificate()
STORED = json.loads(Path(__file__).with_name('results.json').read_text())


def test_every_named_check_passes():
    failed = [name for name, ok in RESULT['checks'].items() if not ok]
    assert not failed and RESULT['all_checks_pass']
    assert {name.split('.')[0] for name in RESULT['checks']} >= {
        '0_gates', '1_constants', '2_relations', '3_map', '4_lattice', '5_rows', '6_endpoint', '7_validation'}


def test_stored_results_are_the_recomputation():
    assert same(STORED, RESULT)


def test_declared_chain_is_a_crossover():
    row = RESULT['rows']['declared']
    assert row['kind'] == 'crossover'
    assert row['x_at_T_c'] > X_C['value'] + 5 * X_C['err']
    assert 150.0 < row['T_c'] < 170.0
    assert row['T_c_err_lattice'] < 1.0 < row['truncation_error'] < 2.0


def test_peak_grows_at_the_first_order_point_and_saturates_in_the_crossover():
    v = RESULT['validation']
    fo = v['first_order']['sets']
    sides = sorted(int(k) for k in fo)
    assert fo[str(sides[-1])]['chi3_max'] > 2.0 * fo[str(sides[0])]['chi3_max']
    lattice = RESULT['lattice_sets']
    heights = sorted((s['n_side'], s['chi3_max']) for s in lattice.values() if s['beta_G'] == 8 and abs(s['x'] - 0.296) < 1e-9)
    assert heights[-1][1] < 1.25 * heights[0][1]


def test_pipeline_reproduces_the_published_standard_model_peak():
    c = RESULT['comparison_not_blind']['standard_model_peak_beta_G_8']
    assert abs(c['this_pipeline'] - c['published_interpolated']) < 0.5
    assert REFERENCE['T_c'] == 159.5
