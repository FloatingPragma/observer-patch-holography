from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo.openwave.row_27_metric.free_source_clock import (
    mixture,
)


def test_short_time_cannot_explain_a_large_flip():
    assert (
        mixture.witness(np.array([-0.9, -0.9]), np.array([-0.8, -0.8]), 0.001) is None
    )


def test_joint_records_and_constructive_witnesses():
    with np.load(Path(__file__).with_name("counts.npz")) as archive:
        result = mixture.compare(archive)
    assert set(result["count_bands_with_witness"].values()) == {193}
    example = result["example_free1_index108"]
    assert min(example["weights"]) >= 0
    np.testing.assert_allclose(sum(example["weights"]), 1)
    np.testing.assert_allclose(
        example["flip_probabilities"], [0.98043, 0.97796], atol=1e-5
    )
    assert result["source_species_tests"]["free1"]["bonferroni_p"] < 4e-10
    assert result["source_species_tests"]["free2"]["bonferroni_p"] > 0.01
    assert result["source_species_tests"]["park1"]["bonferroni_p"] == 1
