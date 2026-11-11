"""Integer lift of the Z6 defect register: exact gates on a reference run, kinematics, readings, imported counts."""
import json
from pathlib import Path

from . import model as m
from .census import CENSUS, parse_word
from .certificate import REFERENCE, result

RESULT = result()


def test_every_named_check_passes():
    failed = [name for name, ok in RESULT['checks'].items() if not ok]
    assert not failed and RESULT['all_checks_pass']
    assert {name.split('.')[0] for name in RESULT['checks']} == {
        '1_closed_carriers', '2_exact_gates', '3_kinematics', '4_readings', '5_reference_run', '6_imported_census'}


def test_closed_carrier_and_stokes():
    fed = m.Federation(1)
    s = fed.structure()
    assert (s['ports'], s['seams'], s['faces'], s['euler_characteristic']) == (12, 30, 20, 2)
    assert s['every_directed_seam_once'] and s['screen_euler_characteristic'] == 2
    link = [(7 * k * k + 3 * k) % 23 - 11 for k in range(fed.seams)]            # any integer link field
    values = m.curl(fed, link)
    assert all(sum(values[20 * c:20 * c + 20]) == 0 for c in range(fed.carriers))
    assert sum(values[fed.intra_faces:]) == 0


def test_exact_gates_on_the_reference_run():
    ref = RESULT['reference_run']
    assert all(ref['gates_at_end'].values()) and ref['sweeps_with_all_gates'] == REFERENCE['sweeps']
    assert set(ref['cuts']) == {str(d) for d in REFERENCE['cuts']}
    for cut in ref['cuts'].values():
        assert all(cut['gates'].values())
        assert cut['twists_are_those_of_depth_at_most_the_cut'] and cut['values_are_the_sum_over_those_twists']
    assert ref['births'] == ref['births_minimal'] > 0
    assert set(ref['applied_by_abs_j']) <= {'1', '2', '3'}


def test_reference_run_equals_the_public_engine():
    ref = RESULT['reference_run']
    stored = CENSUS['reference_cross_check']
    assert (stored['level'], stored['sweeps']) == (REFERENCE['level'], REFERENCE['sweeps'])
    assert ref['digests'] == stored['digests'] and ref['applied_twists'] == stored['applied_twists']
    assert all(stored['public_engine_equal'].values())


def test_the_law_by_hand_and_a_second_execution():
    assert m.attempt(5, 0, 1) == (3, 2, 2) and m.attempt(5, 0, 0) == (2, 3, 3)
    assert m.attempt(2, 3, 1) == (3, 2, -1) and m.attempt(2, 3, 0) == (2, 3, 0) and m.attempt(3, 3, 1) == (3, 3, 0)
    assert m.accepts(0, 0, 3, 5) and m.accepts(0, 0, 1, 3) and not m.accepts(0, 0, 1, 2)
    assert m.accepts(1, 0, -1, -1) and not m.accepts(1, 0, 1, 1) and not m.accepts(3, 0, 1, 1)
    fed = m.Federation(0)
    sweeps = 40
    state = m.run(fed, sweeps)
    x = m.prepare(fed)
    cls, val = [0] * fed.faces, [0] * fed.faces
    for seq, coin in m.schedule(fed, sweeps):                         # the same schedule, executed with the two functions
        for s, c in zip(seq, coin):
            a, b = fed.end_a[s], fed.end_b[s]
            d = x[a] - x[b]
            na, nb, j = m.attempt(x[a], x[b], c)
            x[a], x[b] = na, nb
            p, q = fed.face_plus[s], fed.face_minus[s]
            if j and m.accepts(cls[p], cls[q], j, d):
                cls[p], cls[q] = (cls[p] + j) % 6, (cls[q] - j) % 6
                val[p] += j
                val[q] -= j
    assert x == state['loads'] and cls == state['classes'] and val == state['values']
    assert sum(x) == sum(state['loads_initial']) and len(state['twists']) > 0
    assert all(m.gates(fed, state).values())


def test_kinematics():
    kin = RESULT['kinematics']
    assert kin['released_potential_by_contrast'] == {'0': [], '1': [0], '2': [1], '3': [2], '4': [4], '5': [6]}
    assert kin['births_accepted'] == [[-5, -3], [-5, -2], [-4, -2], [-3, -1], [3, 1], [4, 2], [5, 2], [5, 3]]
    after = kin['value_after_emptying']
    assert after['1'] == after['2'] == after['-1'] == after['-2'] == [0]
    assert after['4'] == after['5'] == [6] and after['-4'] == after['-5'] == [-6]
    assert after['3'] == [0, 6] and after['-3'] == [-6, 0]
    rule = kin['band_rule']
    assert rule['1,0'] == [{'j': -1, 'classes_after': [0, 1]}] and rule['0,0'] == []
    assert all(rule[f'{k},0'] == [] for k in (2, 3, 4))
    assert rule['3,3'] and rule['2,4']                      # adjacent heavy faces react in the band


def test_readings():
    read = RESULT['readings']
    assert read['hypercharge_integers'] == [-6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6]
    assert [r['value'] for r in read['values_by_class']['4']] == [-2, 4]
    assert [r['six_Y_of'] for r in read['values_by_class']['4']] == ['d_R', 'u_R']
    assert [r['value'] for r in read['values_by_class']['1']] == [1]
    assert read['words_with_total_zero']['{4,4,4}'] == ['(+4,-2,-2)']
    assert read['minimal_words_with_total_zero']['{4,4,4}'] == []
    assert 'boson' in read['statistics'] and 'T3' in read['scope']


def test_imported_counts():
    runs = CENSUS['runs']
    assert set(runs) == {'L4', 'L5', 'L6', 'L7'}
    assert CENSUS['source']['status'].startswith('imported run counts')
    for run in runs.values():
        assert run['all_gates_pass']
        for row in run['words_end'].values():
            for word in row['words']:
                values, hidden = parse_word(word)
                assert sum(values) + sum(hidden) == 0
    derived = RESULT['imported_census']['derived']
    for name in ('L5', 'L6', 'L7'):
        d = derived[name]['internal']
        assert 0.996 < d['share_at_a_hypercharge_integer'] < 1 and 0.92 < d['share_minimal'] < 0.94
        assert 0.07 < d['class_4_share_plus_4'] < 0.09
        assert 0.5 < derived[name]['words']['{4,4,4}']['share'] < 0.67
        assert d['null_random_representative']['share_hypercharge_integers'] < 0.66
    assert all(abs(p['standard_errors']) < 3 for p in RESULT['imported_census']['pooled_mirror_pairs'].values())


def test_receipt_matches_recomputation():
    stored = json.loads(Path(__file__).with_name('results.json').read_text())
    for key in ('checks', 'closed_carriers', 'reference_run', 'kinematics', 'readings', 'imported_census', 'input_sha256'):
        assert stored[key] == json.loads(json.dumps(RESULT[key], sort_keys=True)), key
