"""Check the matched preparation, exact observability and actual port records."""

import json

import numpy as np

from . import certify, counts, measure as m


def test_common_start_and_physical_phase_feedback():
    states = m.stocks()
    target = states["target"]
    feedback = states["phase_feedback"]
    np.testing.assert_array_equal(target[0][1], feedback[0][1])
    p = certify.preparations()
    assert p["phase_route"] == [1, 0, 1]
    assert p["endpoint_gated_seam_energy_erased_feedback"] == ["51/16", "17/8"]
    assert p["high_marginal_trace_distances"] == dict(
        coin=0, A=0, B=0, **{"A+B": 0, "coin+A": 0.5, "coin+B": 0.5, "coin+A+B": 0.5}
    )


def test_exact_full_source_law_onset_joint_information_and_all_clock_identity():
    exact = certify.run()
    assert exact == json.loads((m.HERE / "certificates.json").read_text())
    assert [s["dimension"] for s in exact["sectors"]] == [4, 96, 1104]
    for key in ("phase_feedback", "exchange_1_1", "exchange_-1_-1"):
        row = exact["standards"][key]
        assert row["spectral_moments"][0][1] == dict(
            orders=100, first_different_order=None, difference_moment_of_4H=None
        )
        assert (
            row["taylor_even"] == [["0", "0", "0", "-1445/2304", "516125/516096"]] * 2
        )
    expected = ["0", "0", "0", "-1445/2304", "1060235/172032"]
    assert exact["standards"]["target"]["taylor_even"] == [expected] * 2
    assert exact["standards"]["balanced_coherent"]["taylor_even"] == [expected] * 2
    for p in ("2", "13"):
        assert [
            v["power"] for v in exact["real_source_tomography"][p]["independent"]
        ] == [0, 1, 5, 6, 7]
    assert (
        exact["standards"]["phase_feedback"]["spectral_moments"][1][1][
            "difference_moment_of_4H"
        ]
        == "-6288640"
    )


def test_all_stocks_replay_against_native_flights_and_control_records():
    data = m.run()
    archived = json.loads((m.HERE / "results.json").read_text())
    m.fp.check(archived, data)
    stocks = data["standards"]
    target = stocks["target"]
    feedback = stocks["phase_feedback"]
    truth = json.loads((m.fp.HERE / "results.json").read_text())["written"]["8"][
        "probes"
    ]
    for i, t in enumerate(truth):
        np.testing.assert_allclose(
            t["port_contrast"],
            target["probes"][i]["port_response"][:5],
            atol=4e-11,
            rtol=0,
        )
        assert (
            target["probes"][i]["radial"][2] > 0.069
            and feedback["probes"][i]["radial"][2] < -0.134
        )
        assert (
            target["probes"][i]["radial"][4] < -0.34
            and feedback["probes"][i]["radial"][4] > 0.051
        )
        assert feedback["probes"][i]["with_extra_energy_TV"][1] > 0.49
        for a, b in (
            ("exchange_1_1", "phase_feedback"),
            ("exchange_-1_-1", "phase_feedback"),
            ("balanced_coherent", "target"),
            ("balanced_sharp", "target"),
        ):
            np.testing.assert_allclose(
                stocks[a]["probes"][i]["port_response"],
                stocks[b]["probes"][i]["port_response"],
                atol=5e-13,
            )
    for row in (feedback, stocks["exchange_1_1"]):
        np.testing.assert_allclose(
            row["source_only_energy_by_port"],
            target["source_only_energy_by_port"],
            atol=5e-14,
        )
        np.testing.assert_allclose(
            row["local_sector_weights"], target["local_sector_weights"], atol=5e-14
        )
    read = counts.acquire(data)
    assert read == json.loads((m.HERE / "counts.json").read_text())
    for name, panel in read["records"].items():
        for sample, record in zip(panel, stocks[name]["probes"], strict=True):
            assert (
                max(abs(np.array(sample["radial"]) - record["radial"][:5]))
                < read["simultaneous_99_percent_half_width"]
            )
