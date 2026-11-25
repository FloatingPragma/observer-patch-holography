"""Source checks from the endpoint loads and the historical census."""

import numpy as np

from . import source


def test_unobserved_neutral_tag_prepares_the_committed_density_for_both_signs():
    prepared, duration = source.prepare()
    target = np.diag([111, 68]) / 179
    for signs in prepared:
        np.testing.assert_allclose(signs @ signs.conj().T, target, atol=2e-13)
        np.testing.assert_allclose(np.abs(signs) ** 2, target, atol=2e-13)
    np.testing.assert_allclose(prepared[0], prepared[1], atol=2e-13)
    np.testing.assert_allclose(np.cos(duration / 2) ** 2, 111 / 179, atol=2e-14)


def test_skipping_the_copy_preserves_diagonal_and_exposes_the_relative_phase():
    prepared, _ = source.prepare(tagged=False)
    rho = prepared[1] @ prepared[1].conj().T
    np.testing.assert_allclose(np.diag(rho), np.array([111, 68]) / 179, atol=2e-13)
    np.testing.assert_allclose(rho[0, 1], -1j * np.sqrt(111 * 68) / 179, atol=2e-13)
    np.testing.assert_allclose(prepared[:, :, 1], 0, atol=2e-13)


def test_native_pointer_cohorts_discriminate_the_tag_only_on_the_phase_context():
    from . import measure

    r = measure.collect()
    assert r["cohorts"] == 64
    assert (
        r["largest_sample_minus_model_absolute"]
        < r["simultaneous_99pct_yes_frequency_halfwidth"]
    )
    for row in r["contexts"]:
        for mode, outputs in row["instruments"].items():
            for sign in outputs.values():
                tagged, control = sign["tagged"], sign["copy_omitted"]
                bins = np.vstack(
                    (tagged["joint_counts"], tagged["all_other_tag_counts"])
                )
                assert bins.shape == (17, 2) and bins.sum() == r["launches_per_cohort"]
                gap = control["yes_probability"] - tagged["yes_probability"]
                expect = 0
                if row["context"] == "phase":
                    expect = np.sqrt(111 * 68) / 179 * (1 if mode == "sharp" else 0.5)
                np.testing.assert_allclose(gap, expect, atol=2e-10)
