import hashlib
import json

import numpy as np
from scipy.integrate import quad

from .measure import ARCHIVE, HERE, doubling, gram_check, infer_angle, profile


def test_volume_angle_and_positive_spherical_distance_control():
    angles = [0.001, 0.1, 0.5, 1.3, np.pi / 2]
    for angle in angles:
        integral = quad(lambda t: 2 * np.sin(t) ** 2, 0, angle)[0]
        np.testing.assert_allclose(profile(angle), integral, atol=1e-15, rtol=1e-14)
        np.testing.assert_allclose(infer_angle(doubling(angle)), angle, atol=8e-11)
    rng = np.random.default_rng(351)
    points = rng.normal(size=(7, 4))
    points /= np.linalg.norm(points, axis=1)[:, None]
    distances = np.arccos(np.clip(points @ points.T, -1, 1))
    np.fill_diagonal(distances, 0.0)
    check = gram_check(distances, 1.0, list(range(7)))
    assert not check["excludes_five_percent_ranges"]


def test_archived_five_range_interval_certificates_and_transfer():
    data = json.loads((HERE / "results.json").read_bytes())
    original = json.loads(ARCHIVE.read_bytes())
    assert (
        data["native_archive_sha256"]
        == hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
    )
    assert (
        data["protocol_sha256"]
        == hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest()
    )
    for survey, row in zip(original["panels"], data["panels"]):
        witness = row["distance_test"]["strongest_submatrix"]
        ids = witness["indices"]
        ang = (
            np.asarray(survey["ball_landmark_ranges"])[np.ix_(ids, ids)]
            * row["fitted_angle"]
            / row["radius_swaps"]
        )
        vec = np.array(witness["interval_witness_vector"])
        product = np.outer(vec, vec)
        for tolerance, target in (
            (0.05, witness["quadratic_form_ceiling_at_five_percent"]),
            (witness["interval_relative_error_lower_bound"], 0),
        ):
            upper = float(
                np.sum(
                    np.maximum(
                        product * np.cos(ang * (1 - tolerance)),
                        product * np.cos(np.minimum(np.pi, ang * (1 + tolerance))),
                    )
                )
            )
            np.testing.assert_allclose(upper, target, atol=3e-12)
        assert witness["quadratic_form_ceiling_at_five_percent"] < 0
    assert all(row["relative_error"] > 0.05 for row in data["fixed_curvature_transfer"])
