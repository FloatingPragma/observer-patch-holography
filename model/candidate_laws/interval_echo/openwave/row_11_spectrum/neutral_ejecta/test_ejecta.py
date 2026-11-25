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
