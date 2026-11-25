import hashlib
import json
from pathlib import Path

import numpy as np

from . import certificate, current_area, measure, transport


def receipt(filename):
    return json.loads(Path(__file__).with_name(filename).read_text())


def test_joint_rank_and_exact_phase():
    r = certificate.run()
    assert r["stacked_right_inverse_defect"] < 0.000005
    assert max(r["left_inverse_defects"].values()) < 0.00233
    path = Path(__file__).with_name("certificate.json")
    assert (
        hashlib.sha256(path.read_bytes()).hexdigest()
        == "a8d00a55b0dc034ee1bd1cd1ecb6f5ff3a8698a0c371574a2376a99fb085eb0a"
    )
    saved = receipt("certificate.json")["invisible_phase"].copy()
    live = r["invisible_phase"].copy()
    # Different eigensolver proposals have independently certified logarithm
    # errors. Both give the same signed bound and overlap the matched result.
    archived = json.loads(saved.pop("independent_modular_slope"))
    current = json.loads(live.pop("independent_modular_slope"))
    direct = json.loads(live["matched_collar_slope"])
    assert live == saved
    for lower, upper in (archived, current):
        assert 0.00030329 < lower <= direct[0] <= direct[1] <= upper < 0.00034506


def test_separate_probe_measurements():
    result = measure.experiment()
    np.testing.assert_allclose(
        result["common_population_area_volume"],
        [14.2988345301, 2.7634718368],
        rtol=0,
        atol=1e-9,
    )
    archive = receipt("results.json")
    for q in ("-1", "1"):
        r = result["charges"][q]["exposures"][str(2**50)]
        for k, expected in (("ordinary", -0.0177220052), ("pairing", -0.0100027630)):
            arm = r["matched"][k]
            ea, ev = arm["simultaneous_sampling_halfwidth"]
            assert expected - ea < arm["survey_area_slope"] < expected + ea
            assert arm["survey_area_slope"] + ea < 0
            assert abs(arm["survey_volume_slope"]) < ev
            archived = archive["charges"][q]["exposures"][str(2**50)]["matched"][k]
            np.testing.assert_allclose(
                ea, archived["simultaneous_sampling_halfwidth"][0], rtol=1e-4
            )


def test_spatial_transport():
    r = transport.run()
    for actual, archived in zip(
        r["uniform_vertex_ball_reads"],
        receipt("transport.json")["uniform_vertex_ball_reads"],
    ):
        assert actual["ball_counts"] == archived["ball_counts"]
        np.testing.assert_allclose(
            actual["curvature"], archived["curvature"], atol=1e-10
        )


def test_conserved_current_against_independent_area():
    record = current_area.run()
    for q in ("-1", "1"):
        for source, arm in record["charges"][q].items():
            np.testing.assert_allclose(
                arm["area_prediction"], arm["source_prediction"], atol=2e-12
            )
            assert (
                abs(arm["population_current_per_time"][-1] - arm["area_prediction"])
                < 2e-6
            )
            assert (
                abs(
                    arm["marginal_contribution"]
                    + arm["collar_contribution"]
                    - arm["source_prediction"]
                )
                < 1e-12
            )
            if source == "pairing":
                assert abs(arm["marginal_contribution"]) < 1e-12
