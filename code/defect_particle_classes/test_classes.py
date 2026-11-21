import json
from pathlib import Path

from . import model as m
from . import statistics as st
from .certificate import result

RESULT = result()


def test_every_named_check_passes():
    failed = [name for name, ok in RESULT['checks'].items() if not ok]
    assert not failed and RESULT['all_checks_pass']
    assert {name.split('.')[0] for name in RESULT['checks']} == {
        '1_class_rule', '2_z6_kernel', '3_class_table', '4_aharonov_bohm', '5_reactions', '6_closed_carriers',
        '7_spin_statistics', '8_imported_census', '9_central_flux_costs'}


def test_central_flux_costs():
    costs = RESULT['central_flux_costs']
    assert [costs[str(k)]['beta'] for k in range(6)] == [None, '12/5', '3', '4', '3', '12/5']
    assert costs['3']['eigenphases_over_pi'] in (['0', '0', '0', '1', '-1'], ['0', '0', '0', '-1', '1'])
    assert costs['3']['trace_over_pi_squared'] == '2'


def test_class_rule_and_multiplicities():
    rows = m.patterns()
    assert len(rows) == 32
    for r in rows:
        assert r['k'] == (3 * (r['n2'] % 2) - 2 * (r['n3'] % 3)) % 6
        assert (r['k'] % 3, r['k'] % 2) == (r['n3'] % 3, r['n2'] % 2)
    assert RESULT['class_multiplicities'] == [4, 6, 6, 4, 6, 6]


def test_class_table():
    t = RESULT['class_table']
    assert t['1']['fermion_multiplets'] == ['Q'] and t['5']['fermion_multiplets'] == ['anti-Q']
    assert t['4']['fermion_multiplets'] == ['d_R', 'u_R'] and t['2']['left_handed_names'] == ['d^c', 'u^c']
    assert t['3']['left_handed_names'] == ['L'] and t['3']['bosons'] == ['H', 'anti-H']
    assert t['0']['left_handed_names'] == ['e^c', 'nu^c'] and t['0']['bosons'] == ['B', 'G', 'W']
    assert RESULT['left_handed_generation_classes'] == {'Q': 1, 'u^c': 2, 'd^c': 2, 'L': 3, 'e^c': 0, 'nu^c': 0, 'H': 3}


def test_flux_visibility():
    seen = RESULT['aharonov_bohm']['classes_seeing_flux']
    assert seen['2'] == seen['4'] == [1, 2, 4, 5]
    assert seen['3'] == [1, 3, 5]
    assert seen['1'] == seen['5'] == [1, 2, 3, 4, 5] and seen['0'] == []


def test_reactions_and_yukawa_vertices():
    react = RESULT['reactions']
    assert sum(r['in_sixteen_types'] for r in react.values()) == 16
    assert all(r['conserves_class'] and r['realizations'] for r in react.values())
    for name in ('1+3->4', '3+5->2', '1+2->3', '4+5->3', '3+3->0'):
        assert 'Yukawa' in react[name]['realization_vertex_types']
    assert react['1+1->2']['realizations'][0]['fields'] == ['Q', 'Q', 'anti-d_R']
    assert all(v['neutral'] and v['hypercharge_sum'] == '0' for v in RESULT['yukawa_operators'].values())


def test_closed_carriers():
    stokes = RESULT['closed_carriers']['stokes']
    assert (stokes['faces'], stokes['seams'], stokes['seam_column_sums']) == (20, 30, [0])
    for content in m.CARRIER_CONTENTS:
        assert sum(content) % 6 == 0


def test_central_involution_and_patterns():
    spin = RESULT['spin_statistics']
    assert spin['cover']['order'] == 120 and spin['cover']['centre_size'] == 2
    assert spin['z_on_code'] == [1, 1, 1, -1, -1]
    assert all(p['z'] == p['minus_one_to_n2'] == p['class3_holonomy'] for p in spin['z_on_patterns'])
    assert st.central_involutions_of_G() == [[1, 1, 1, -1, -1], [1, 1, 1, 1, 1]]


def test_rotation_grading():
    g = RESULT['spin_statistics']['gradings']
    assert g['number_of_gradings'] == 16
    assert g['defect_class_gradings'] == [[0, 0, 0, 0, 0], [1, 0, 1, 0, 1]]
    rot = g['rotation_grading']
    assert set(rot['quanta_by_block'].values()) == {1}
    assert rot['register'] == [0, 0] and set(rot['defects'].values()) == {0} and rot['links'] == [0]
    assert g['rules_that_are_gradings']['degree'] is False
    comp = RESULT['spin_statistics']['composites']
    assert comp['baryon_three_quanta'] == [1] and comp['meson_quark_antiquark'] == [0]


def test_imported_counts_are_labelled():
    census = RESULT['imported_census']
    assert census['source']['status'].startswith('imported run counts')
    assert all(census['carrier_contents_class_neutral'].values())


def test_receipt_matches_recomputation():
    stored = json.loads(Path(__file__).with_name('results.json').read_text())
    for key in ('checks', 'class_table', 'reactions', 'spin_statistics', 'closed_carriers', 'imported_census'):
        assert stored[key] == json.loads(json.dumps(RESULT[key], sort_keys=True)), key
