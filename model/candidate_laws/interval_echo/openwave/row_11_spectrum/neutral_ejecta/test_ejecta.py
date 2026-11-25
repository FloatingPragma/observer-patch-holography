from fractions import Fraction
from math import comb, factorial

from model.candidate_laws.interval_echo.openwave.row_11_spectrum.neutral_ejecta import (
    measure,
)


def test_full_board_first_word_and_geodesic_count():
    for level in range(3):
        result = measure.full_board(level)
        # Each arm takes one inward and one outward internal move between two
        # bridge crossings. These internal moves can occur in either order.
        paths = 2**2 * comb(8, 4)
        assert (result["first_nonzero_amplitude_order"], result["shortest_paths"]) == (
            8,
            paths,
        )
        probability = Fraction(paths, 2**8 * factorial(8)) ** 2
        assert probability == Fraction(
            result["leading_full_word_probability_coefficient"]
        )
        holes, doubles = result["observed_holes_doubles"]
        assert sum(v < 12 for v in holes) == 1 and all(v >= 12 for v in doubles)
        assert sorted(v // 12 for v in holes if v >= 12) == sorted(
            v // 12 for v in doubles
        )
        assert len({v // 12 for v in doubles}) == 2


def test_native_support_and_changed_law_with_same_readout():
    certificate = measure.local_certificate()
    assert certificate["terminal_port_minus_zero"] == [0, 1, 2]
    bench = measure.bench()
    for setting in bench["settings"]:
        for panel in setting["panels"]:
            assert (
                sum(panel["none_left_right_both_counts"]) == bench["shots_per_setting"]
            )
            if setting["conversion_strength"] or (
                setting["seeded_neutral_pairs"] == 2 and not setting["bridge_cut"]
            ):
                assert panel["both_probability"] > 0
            else:
                assert panel["both_probability"] == 0
                assert panel["none_left_right_both_counts"][-1] == 0


def test_rest_ratios_need_large_neutral_receivers_or_an_exterior_record():
    rows = measure.formation_receiver_bounds()
    for ratio, pairs, minimum_in_one in (("206.7682827", 41, 21), ("3477.0", 693, 347)):
        assert rows[ratio] == dict(
            minimum_defect_ports=2 * pairs + 1,
            minimum_pairs_in_two_exhaustive_receivers=pairs,
            minimum_pairs_in_one_of_the_two_receivers=minimum_in_one,
            minimum_ports_in_that_neutral_receiver=2 * minimum_in_one,
            minimum_doubles_outside_two_single_cell_receivers=pairs - 12,
        )
        assert (
            min(max(left, pairs - left) for left in range(pairs + 1)) == minimum_in_one
        )
