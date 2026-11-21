"""Check the published geometry against both local records and the source."""

from fractions import Fraction
from itertools import combinations
from math import sqrt

from .six_station import verify_archive
from .station_ranging import report


def test_source_and_count_geometry_on_the_frozen_acquisition():
    arrays = verify_archive()
    reads = report(arrays)
    for pulse, row in enumerate(reads):
        assert sum(row["range_totals"].values()) == 30 * 2**34
        assert row["range_totals"]["both_clicked"] == (171311326, 109477123)[pulse]
        zero, three, five, seven, ten = row["geometry"]
        assert (zero["six_station_verdict"], three["six_station_verdict"],
                five["six_station_verdict"], seven["six_station_verdict"]) == ("excluded",) * 4
        assert ten["six_station_verdict"] == ("excluded", "fit")[pulse]
        assert five["excluded_subsets"] == 0
        x = five["any_dimension_obstruction"]["zero_sum_coefficients"]
        assert sum(x) == 0
        quadratic = Fraction(0)
        for i, j in combinations(range(6), 2):
            lo, hi = map(Fraction, row["ranges"][f"{i}-{j}"])
            c = 2*x[i]*x[j]
            worst = lo*Fraction(19,20) if c > 0 else hi*Fraction(21,20)
            quadratic += c*worst**2
        assert quadratic == Fraction(five["any_dimension_obstruction"]["quadratic_form_lower"]) > 0

        source = row["exact_source_comparison"]
        assert source["source_ranges_outside_or_missing_from_count_bands"] == []
        assert source["indefinite"]
        a, b = source["adjacent_ticks"], source["opposite_ticks"]
        best = (b - sqrt(2)*a)/(b + sqrt(2)*a)
        source_low, source_high = map(Fraction, source["optimal_euclidean_distortion_enclosure"])
        assert float(source_low) <= best <= float(source_high)
        bracket = row["distortion_bracket"]
        assert Fraction(bracket["lower"]) <= source_low <= source_high <= Fraction(bracket["upper"])
        expected_min, expected_max = ((.131, .133), (.079, .080))[pulse]
        assert expected_min < best < expected_max
        for name in ("parked", "detector_disabled"):
            assert all(m["both_clicked"] == 0 and m["receiver_missed_sender_clicked"] > 0
                       for m in row["tag_totals"][name])
        assert all(m["returned_receiver_letter"] + m["changed_receiver_letter"] == m["both_clicked"]
                   for m in row["tag_totals"]["active"])
