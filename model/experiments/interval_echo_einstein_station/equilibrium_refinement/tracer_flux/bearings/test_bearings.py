from fractions import Fraction
import numpy as np
from . import displacements, measure


def test_short_flight_classification_on_all_three_stocks():
    receipt = measure.simulate()
    assert receipt["word_count"] == 660
    assert Fraction(receipt["rotating_frame_bound"]) == Fraction(441, 4)
    for release in receipt["exposures"]:
        threshold = release["decision_probability"]
        assert (
            Fraction(release["universal_nonedge_upper"])
            < threshold
            < Fraction(release["universal_edge_lower"])
        )
        for stock in release["counts"].values():
            assert stock["recovered_edges"] == 30 and stock["gram_rank"] == 3
            assert (
                stock["largest_nonedge_upper"]
                < threshold
                < stock["smallest_edge_lower"]
            )
            assert all(
                sum(row) == receipt["shots_per_launch"] for row in stock["counts"]
            )


def test_record_stencil_reconstructs_the_inner_bearing_mean():
    r = displacements.run()
    small, large = r["shells"]["-1"], r["shells"]["1"]
    assert len(r["shared_origin_records"]) == 5
    assert (
        len(small["antipodal_word_pairs"]) == len(large["antipodal_word_pairs"]) == 10
    )
    np.testing.assert_allclose(small["tilted_total_weight"], (5 - np.sqrt(5)) / 6)
    assert small["positive_isotropic_weights"]
    assert not large["positive_isotropic_weights"]
    assert large["max_paired_squared_axial_cosine"] < 1 / 3
    n, d = r["apparent_scalar_times_radius_squared"]
    assert Fraction(n, d) == Fraction(-264, 65)
    lo, hi = r["simultaneous_99_count_interval"]
    assert lo < n / d < hi < 0
