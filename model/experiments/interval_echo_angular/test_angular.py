import numpy as np

from model.executable_model import cell
from .__main__ import angular_coverage, match_detectors, reconstruct


def test_detector_matching_recovers_the_chart_from_counts_and_abstains_on_ambiguity():
    # Different row and column addresses, with unequal per-source budgets.
    counts = np.array([[1, 15, 2], [6, 0, 1], [3, 4, 7]], dtype=np.uint64)
    order = match_detectors(counts)
    np.testing.assert_array_equal(order, [1, 0, 2])
    np.testing.assert_array_equal(match_detectors(np.iinfo(np.uint64).max * np.eye(3, dtype=np.uint64)),
                                  np.arange(3))
    counts[0, 0] = 3
    assert match_detectors(counts) is None  # exactly half the column
    counts[0, 0] = 20
    assert match_detectors(counts) is None  # two sources claim the same row


def test_reciprocal_readout_abstains_on_uniform_outcomes():
    assert reconstruct(np.ones((12, 12), dtype=int)) is None


def test_omitting_an_antipodal_pair_loses_directional_information():
    counts = 20 * np.eye(12, dtype=int) + cell().adjacency().toarray().astype(int)
    _, frame = reconstruct(counts)
    report = angular_coverage(frame, range(12))
    assert report["subset_size_to_rank_counts"]["10"] == {"8": 6, "9": 60}
    assert report["antipodal_average_rank"] == 6


def test_overlap_logarithms_close_only_on_the_declared_slow_sector():
    from .modular import run as measure

    result = measure()
    assert result["slow_state_variation_rank"] == 9
    assert result["pure_triplet_observable_rank"] == 6
    assert result["full_one_load_observable_rank"] == 12
    assert result["max_axial_form_error"] < 1e-9
    assert result["slow_sector_held_operator_error"] < 1e-9
    assert result["single_service_order_slow_sector_held_error"] > .01
    controls = [row for row in result["readings"] if row["source"] == "sharp_port"]
    assert len(controls) == 2
    assert min(row["max_abs_exact_held_miss"] for row in controls) > 5.
