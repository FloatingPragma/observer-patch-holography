import numpy as np
import pytest

from . import intervention, material_radar


def test_control_with_equal_clock_and_propagation_laws_has_no_certified_response():
    counts = np.zeros((2, 2, 7, 2, 2, 2, 7), dtype=np.int64)
    counts[..., 0, 0] = 500
    counts[..., 1, 0] = 500
    spec = {"endpoints": [0, 5, 11, 16, 22, 29], "mark": 2,
            "phase_kick": {"shots": 1000, "alpha": 1e-6, "post_waits": [2, 6]}}
    assert intervention.read(spec, counts)["strongest_joint_lower_bound"] == 0
    counts[..., 0, 0] += 1
    with pytest.raises(ValueError):
        intervention.read(spec, counts)


def test_receipt_rejects_missing_switched_and_nonfinite_outcomes(tmp_path):
    expected = {"radar_joint_counts": np.array([4, 1], dtype=np.int64),
                "coherent_radar_overlap": np.array([.5j])}
    for changed in ({"radar_joint_counts": expected["radar_joint_counts"]},
                    dict(expected, radar_joint_counts=np.array([1, 4], dtype=np.int64)),
                    dict(expected, coherent_radar_overlap=np.array([complex(float("nan"), 0)]))):
        path = tmp_path / "archive.npz"
        np.savez(path, **changed)
        with pytest.raises(ValueError):
            material_radar.verify_counts(path, expected)
