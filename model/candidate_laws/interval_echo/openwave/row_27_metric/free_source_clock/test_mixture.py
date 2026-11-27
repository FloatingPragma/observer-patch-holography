from pathlib import Path

import numpy as np
import pytest

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


def test_missing_timer_or_outcome_cannot_shorten_the_registered_experiment():
    with np.load(Path(__file__).with_name("counts.npz")) as archive:
        key = "k8_even_blank"
        records = {key: archive[key], key + "_timer": archive[key + "_timer"][:-1]}
        with pytest.raises(ValueError, match="Incomplete timing record"):
            mixture.compare(records)
        records[key] = np.column_stack((archive[key], np.zeros(len(archive[key]), dtype=int)))
        with pytest.raises(ValueError, match="Incomplete coin and endpoint record"):
            mixture.compare(records)
