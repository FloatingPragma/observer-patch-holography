"""Independent clock marginals, calibrated kinetic mass and held sign tests."""

import json

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.regional_reference.recovery import (
    equipped_mass,
    inward_check,
    joint_records as joint,
    mobile_standard,
)


def archived(name):
    return json.loads((joint.HERE / (name + ".json")).read_text())


def test_kinetic_calibration_predicts_both_held_fringes():
    output = mobile_standard.run()
    joint.source.region.check(archived("mobile_standard"), output)
    for law, expected in (("native", 1), ("half_exchange", 2)):
        res = output["laws"][law]
        lo, hi = res["kinetic_mass_interval"]
        assert lo < expected < hi
        assert all(row["passed"] for row in res["held"])


def test_blank_kinetic_coefficient_on_the_clock_equipped_reference():
    result = equipped_mass.run()
    joint.source.region.check(archived("equipped_mass"), result)
    assert sum(sum(r["excluded_components"]) for r in result["rows"]) == 11
    assert all(r["pure_phase_excluded"] for r in result["rows"])
    assert result["rows"][-1]["length_interval"][0] > 1


def test_joint_clock_marginals_match_the_separate_timing_instrument():
    run = archived("joint_records")
    clock = json.loads((joint.source.timing.HERE / "results.json").read_text())
    scan = np.array(clock["scan"]["signed_return_by_hand"])
    ticks = [clock["scan"]["times"].index(t) for t in (2, 4)]
    probs = {
        name: ((1 + scan[:, ticks, i]) / 2).mean(axis=1)
        for i, name in enumerate(clock["names"])
    }
    for k, x, y in (
        ("baseline", "empty", "paired"),
        ("unit", "mass+", "paired"),
        ("matched", "empty", "feedback"),
    ):
        h = np.array(run["laws"]["8"]["counts"][k])[0] / run["trials_per_bank"]
        recorded = np.array([h[0].sum(), h[:, 0].sum()])
        np.testing.assert_allclose(
            recorded, (probs[x] + probs[y]) / 2, atol=0.00015, rtol=0
        )


def test_all_agreements_and_held_inward_reads():
    run = archived("joint_records")
    native = run["laws"]["8"]
    delta = np.array(run["changes_simultaneous_99_half_width"])
    weights = np.array([1, 0, 10, 0, 0]) / 11
    scores = {k: weights @ v for k, v in native["measured_changes"].items()}
    assert scores["unit"] - weights @ delta > 0.0001997
    assert max(scores["matched"], scores["unmatched"]) + weights @ delta < -0.0001374
    assert native["certified_both_inward"]["u_bounds"][1] > 0.69
    for law in run["laws"].values():
        assert all(v["passed"] for v in law["held_forecasts"].values())
        assert law["half_load_control"]
        for rows in law["counts"].values():
            assert np.all(np.sum(rows, axis=(1, 2, 3)) == joint.SHOTS)
    held = archived("inward_check")
    assert (
        0.67
        < held["source_agreement_interval"][0]
        < held["source_agreement_interval"][1]
        < 0.69
    )
    assert max(held["bare_connected_radius_hops"]) + held["bare_half_width"] < 0
    load, source = (held["equipped"][k] for k in ("unit", "held"))
    assert (
        max(source["connected_radius_shifts_hops"]) + source["radius_half_width_hops"]
        < 0
    )
    assert (
        load["witness"] - load["witness_half_width"]
        > 0
        > source["witness"] + source["witness_half_width"]
    )


def test_new_archives_bind_their_protocols():
    from hashlib import sha256

    for module, stem, protocol in (
        (joint, "joint_records", "JOINT_RECORDS.md"),
        (inward_check, "inward_check", "INWARD_CHECK.md"),
    ):
        assert (
            archived(stem)["protocol_sha256"]
            == sha256((module.HERE / protocol).read_bytes()).hexdigest()
        )
