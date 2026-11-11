"""Fast checks of section 8 (many rounds on the sixty-address carrier)."""
from fractions import Fraction as F

try:
    from . import rounds
except ImportError:
    import rounds

DATA, CHECKS = rounds.result()


def test_every_named_check_passes():
    failed = [name for name, ok in CHECKS.items() if not ok]
    assert not failed, failed
    assert len(CHECKS) == 19


def test_complete_record_of_the_return_tests():
    assert DATA['test_outcome_weights'] == {'yes_yes': '4/81', 'x_no': '11/324', 'y_no': '125/1296',
                                            'both_no': '55/1296'}
    total = sum(F(v) for v in DATA['test_outcome_weights'].values())
    assert total == F(2, 9)


def test_round_without_reset():
    assert DATA['read_weights_of_a_cell_after_a_no'] == {'X': '11/12', 'Y': '5/6'}
    assert DATA['counts_in_the_four_test_states'] == ['2/9', '7/9', '7/9', '2/9']


def test_partners_that_share_the_kick():
    assert DATA['shared_kick_generator_means'] == {'1': '2/9', '2': '2/9', '3': '2/9', '6': '2/9'}
    assert DATA['shared_kick_polynomials']['1'] == ['7/9', '2/9']
    assert DATA['at_once_counts']['2'] == '73/216' and DATA['at_once_counts']['3'] == '563/1296'
    assert DATA['shared_kick_phases']['n=1, N=1'].startswith('0.20533403097')
    assert DATA['shared_kick_phases']['n=2, N=1'].startswith('0.21573295580')


def test_forward_only_control():
    assert DATA['forward_only_counts'][0] == '2/9'
    assert [round(float(F(c)), 3) for c in DATA['forward_only_counts'][1:]] == [0.486, 0.536, 0.467]


def test_register_on_sixty_addresses():
    assert DATA['double_load_register'] == {'a^2': '5/108', 'b^2': '5/216', 'trace': '5/18'}
    assert DATA['native_pair_register'] == {'a^2': '1/27', 'b^2': '1/54', 'trace': '2/9'}


def test_executed_rounds():
    run = DATA['executed_rounds']
    assert abs(float(run['phase_per_unit_clock_time']) - float(run['closed_form'])) < 1e-11
    assert abs(float(run['closed_form']) - (2 / 9 - 35 / (2187 * 100 ** 2))) < 1e-11
