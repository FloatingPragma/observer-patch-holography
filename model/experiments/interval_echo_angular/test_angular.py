import numpy as np

from model.executable_model import cell
from .__main__ import angular_coverage, reconstruct


def test_reciprocal_readout_abstains_on_uniform_outcomes():
    assert reconstruct(np.ones((12, 12), dtype=int)) is None


def test_omitting_an_antipodal_pair_loses_directional_information():
    counts = 20 * np.eye(12, dtype=int) + cell().adjacency().toarray().astype(int)
    _, frame = reconstruct(counts)
    report = angular_coverage(frame, range(12))
    assert report["subset_size_to_rank_counts"]["10"] == {"8": 6, "9": 60}
    assert report["antipodal_average_rank"] == 6
