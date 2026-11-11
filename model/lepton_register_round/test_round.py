import json
from pathlib import Path

from . import classes, walk as wk
from .certificate import result

VALUE = result()


def test_every_named_check_passes():
    failed = [name for name, ok in VALUE['checks'].items() if not ok]
    assert not failed and VALUE['all_checks_pass']
    assert {name.split('.')[0] for name in VALUE['checks']} == {
        '1_register', '2_double_loads', '3_weak_round', '4_clock_unit', '5_kick', '6_carrier', '7_classes'}


def test_register_compression_and_read_count():
    reg = VALUE['register']
    comp = reg['compression']
    assert comp['0']['A'] == [[32, 0, 0], [0, 32, 0], [0, 0, 32]]
    assert comp['+']['B_sqrt2'] == [[0, 0, 16], [16, 0, 0], [0, 16, 0]]
    assert comp['-']['B_sqrt2'] == [[0, 16, 0], [0, 0, 16], [16, 0, 0]]
    assert reg['plain_record_compression']['-']['B_sqrt2'] == [[0, -16, 0], [0, 0, -16], [-16, 0, 0]]
    assert reg['b_over_a_squared'] == '1/2' and reg['eigenvalues']['koide_value'] == '2/3'
    assert reg['object_weights']['C'] == {'idle': '2/3', 'agree': '1/6', 'repair +': '1/12', 'repair -': '1/12'}
    assert reg['read_counts'] == {'P(R+)': '1/12', 'P(H)': '1/6', 'P(R+ and H)': '1/72', 'P(R+ xor H)': '2/9'}
    assert reg['register_outputs']['mu_over_e'].startswith('206.77031')
    assert reg['register_outputs']['tau_over_e'].startswith('3477.472')


def test_geometry_of_the_register():
    geo = VALUE['register']['geometry']
    assert geo['hexagon_in_screen_sense'] == [2, 3, 9, 1, 0, 10]
    assert geo['c_seam_screen_orientation'] == [2, 3] and geo['e_seam_screen_orientation'] == [10, 2]
    assert geo['e_type_seams'] == [[0, 1], [2, 10], [3, 9]]


def test_double_loads_and_weak_round():
    tick = VALUE['double_loads']['tick_one']
    for frame in ('plain', 'oriented'):
        assert tick[frame]['symmetric double'] == {'backward': '1/12', 'double': '2/3', 'forward': '1/12',
                                                   'in-dimer': '1/6'}
        assert tick[frame]['P(R+ xor H) of two objects'] == '2/9'
    rnd = VALUE['weak_round']
    for frame in ('plain', 'oriented'):
        r = rnd['round'][frame]
        assert (r['P'], r['<SS|X>'], r['<X|X>']) == ('2/9', '2/9', '2/9')
        assert r['non_return_A']['norm2'] == '11/144' and r['non_return_B']['norm2'] == '5/36'
    assert rnd['rounds']['1']['phase per unit clock time'] == '0.205334030970278'
    assert rnd['rounds']['10']['phase per unit clock time'] == '0.222062099670503'
    assert rnd['rounds']['1000']['phase per unit clock time'] == '0.222222206218563'
    assert rnd['phase_series'] == '2/9 - 35*x**2/2187' and rnd['log_modulus_times_N'] == '-7/81'
    assert [r['count'] for r in rnd['forward_only_oriented']] == ['2/9', '35/81', '38/81', '35/81', '2/9', '0']
    assert rnd['forward_only_mean'] == '8/27'
    assert rnd['public_one_load_round']['return amplitude'] == '11/12 + I/12'
    assert rnd['public_one_load_round']['no weight'] == '11/72'


def test_clock_unit():
    clk = VALUE['clock_unit']
    gate = {r['tau'] for r in clk['edge_rows'] if r['gate']}
    assert gate == {'1', '5', '-3', '9', '-7'}
    assert clk['principal_branch_members_of_1_plus_4Z'] == [1]
    assert clk['one_minus_overlap_and_record_one_weight'] == [['1/12', '1/12']]


def test_kick_and_crossing():
    k = VALUE['kick']
    assert k['current_on_flags'] == {'C idle': [0], 'C repairs +': [1], 'C repairs -': [-1]}
    assert k['crossing']['plain_law_crosses'] == k['crossing']['starts'] == 120
    assert k['one_port_round']['mean coefficient'] == '2/9'


def test_carrier_counts():
    car = VALUE['carrier']
    assert car['double_load_at_port_2'] == {'blank': '5/6', 'record 2': '1/12', 'record 3': '1/12'}
    assert car['two_cells_at_port_2']['count'] == '2/9'
    assert car['one_register_two_doubles'] == [{'P(R+)': '1/12', 'P(H)': '1/6', 'P(R+ and H)': '0', 'count': '1/4'}]
    assert car['two_nine_times_240'] == '160/3'


def test_class_selection():
    table = VALUE['classes']['consensus_states']
    assert [(table[str(k)]['separated loads'], table[str(k)]['pair on one port']) for k in range(6)] == [
        (2, 1), (0, 0), (0, 0), (0, 1), (0, 0), (0, 0)]
    assert VALUE['classes']['uniform_double_count'] == '2/9'


def test_eisenstein_arithmetic():
    w = classes.Eis(0, 1)
    assert w * w == w - 1 and classes.w_power(6) == classes.Eis(1) and classes.w_power(3) == classes.Eis(-1)
    assert (w * w.conjugate()) == classes.Eis(1) and w.norm() == 1
    assert classes.Eis(2, 3) / classes.Eis(2, 3) == classes.Eis(1)


def test_completion_is_the_public_operator():
    u = wk.completion()
    assert len(u) == 24 and all(len(row) == 24 for row in u)
    assert {str(x) for row in u for x in row} == {'-1/2', '0', '1', '1/2'}


def test_receipt_matches_recomputation():
    stored = json.loads(Path(__file__).with_name('results.json').read_text())
    for key in ('checks', 'register', 'double_loads', 'weak_round', 'clock_unit', 'kick', 'carrier', 'classes',
                'premises'):
        assert stored[key] == json.loads(json.dumps(VALUE[key], sort_keys=True)), key
