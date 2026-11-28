"""Cross-check the finite source, spatial solver and independent material flight receipts."""

import json
import numpy as np

from . import exact_guide, measure


def test_native_encoding_energy_and_closed_transfer():
    saved = json.loads((measure.HERE / "results.json").read_text())
    z = saved["input_return"]["8"]["d"][-1]
    gain, rows = measure.local_source(z)
    assert gain == saved["ordinary_unlike_minus_equal"]
    for orientation in (1, -1):
        group = [
            row
            for row in rows
            if row["orientation"] == orientation and row["switches"] == "independent"
        ]
        pulse = [row["pulse_probability"] for row in group]
        np.testing.assert_allclose(
            np.array(pulse) - pulse[0],
            gain / (2 * 512) * orientation * z * z * np.array([0, 0.75, 1]),
            atol=2e-17,
        )
    transfer = measure.conservative_preparation(z)
    assert transfer["reader_energy_change"] > 0
    assert abs(transfer["summed_ticket_change"]) < 5e-16


def test_all_exact_shells_against_independent_guide_and_probe_receipts():
    exact = exact_guide.run()
    saved = json.loads((measure.HERE / "results.json").read_text())
    archived = json.loads((measure.HERE / "exact_guide.json").read_text())
    assert exact == archived
    z = saved["input_return"]["8"]["d"][-1]
    short, _ = measure.guide_panel(10, 1 / 8, z, 8.5)
    gain = 17 * z * z / 2048
    for s in short["shells"]:
        rational = exact["shells"][s["k"] - 1]
        assert s["edges"] == rational["edges"]
        np.testing.assert_allclose(
            s["enclosed_tickets"],
            gain * rational["enclosed_squared_admission"],
            rtol=3e-14,
        )
    for row in saved["paired_acquisition"]["panels"]:
        counts = row["arms"]
        assert (
            sum(sum(arm["starts"]) for arm in counts)
            == saved["paired_acquisition"]["visits"]
        )
        assert row["signed_per_expected_excess"] + row["simultaneous_99_width"] < 0


def test_original_head_covariance_and_receiver_sign_from_time_four_evolution():
    from . import formation_certificate, formed_heads

    actual = formed_heads.run()
    for row in actual["source"]:
        if row["orientation"] != "matching":
            continue
        certificate = formation_certificate.certify(row["source_sign"])
        for key, value in [
            ("covariance_billionths", row["covariance_of_original_load_bits"]),
            ("native_energy_increment_billionths", row["erased_minus_formed_energy"]),
        ]:
            lower, upper = certificate[key]
            assert 0 < lower < 1e9 * value < upper
        opposite = next(
            r
            for r in actual["source"]
            if r["source_sign"] == row["source_sign"] and r["orientation"] == "opposed"
        )
        np.testing.assert_allclose(
            opposite["modular_change_fixed_formed"], row["modular_change_fixed_formed"]
        )
        np.testing.assert_allclose(
            opposite["erased_minus_formed_energy"], -row["erased_minus_formed_energy"]
        )
    assert all(
        p["signed_per_expected_excess"] + p["simultaneous_99_width"] < 0
        for source in actual["native_receiver_acquisitions"]
        for p in source["panels"]
    )
