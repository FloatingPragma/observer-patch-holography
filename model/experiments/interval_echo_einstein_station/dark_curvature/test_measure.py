from hashlib import sha256
import json

import numpy as np

from .measure import HERE, SHOTS, TIMES, oriented


def test_signed_volume_is_affine_in_anchor_and_each_terminal_vertex():
    rng = np.random.default_rng(455)
    v, other = rng.normal(size=(4, 3)), rng.normal(size=3)
    for j in range(4):
        held, mixed = v.copy(), v.copy()
        held[j] = other
        mixed[j] = 0.3 * other + 0.7 * v[j]
        np.testing.assert_allclose(oriented(mixed), 0.3 * oriented(held) + 0.7 * oriented(v), atol=1e-14)


def test_records_recover_initial_bound_and_held_signed_distortion():
    p = json.loads((HERE / "results.json").read_text())
    assert p["protocol_sha256"] == sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest()
    assert p["source_table_sha256"] == sha256((HERE.parent / "shared_variance.json").read_bytes()).hexdigest()
    counts = np.array(p["histogram_counts"])
    np.testing.assert_array_equal(counts.sum(axis=-1), SHOTS * np.ones((4, len(TIMES)), dtype=int))
    means = counts @ np.array(p["signed_volume_outcomes"]) / (SHOTS * p["initial_volume"])
    np.testing.assert_allclose(means, p["record_signed_volume_over_release_by_arm_time"], atol=1e-14)
    clk = np.array(p["clock_times"])
    w = np.linalg.solve(np.array([np.ones(3), clk[:3]**2, clk[:3]**4]), [0, 2, 0])
    recovered = -means[:, :3] @ w
    np.testing.assert_allclose(recovered, p["record_focusing"], atol=1e-12)
    corrected = recovered - p["conditional_stencil_bias"]
    corrected -= corrected[2]
    np.testing.assert_allclose(corrected, p["corrected_record_focus_contrasts_vs_fresh"], atol=1e-12)
    focus_radius = p["simultaneous_99pct_radius_for_focus_contrast_given_clock"]
    assert (abs(corrected[3]) + focus_radius) / (abs(corrected[0]) - focus_radius) < 0.00249
    np.testing.assert_allclose(p["exact_launch_focusing"][1:], p["exact_launch_focusing"][2], atol=1e-11)
    vol_radius = p["simultaneous_99pct_radius_for_record_volume_difference"]
    assert means[1, -1] - means[2, -1] > vol_radius
    assert means[3, -1] - means[2, -1] < -vol_radius
