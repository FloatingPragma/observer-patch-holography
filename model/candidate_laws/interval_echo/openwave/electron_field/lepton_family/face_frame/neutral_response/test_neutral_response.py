"""Test signed flux, held source classification and transmission of the same event."""

import json
from collections import Counter
import numpy as np
from . import measure


def test_calibrated_currents_and_native_material_report():
    fresh = measure.run()
    saved = json.loads((measure.HERE / "records.json").read_text())
    for panel in (saved, fresh):
        assert panel["report_owners"][0] != panel["report_owners"][1]
        rows = panel["held_sources"]
        assert len(rows) == len(panel["source_key_after_scoring"]) == 11
        for station in ("alice", "bob"):
            refs = panel["source_calibration"][station]
            for charge in (-1, 0, 1):
                assert (
                    refs["charge_standards"][str(charge)][charge + 1] == measure.SHOTS
                )
            for sign in (-1, 1):
                ref = refs["reciprocal_charge_to_neutral"][str(sign)]
                assert ref["charged_standard_input"] == -sign
                lo, hi = ref["efficiency99"]
                assert lo <= np.sin(measure.EXPOSURES[station][sign]) ** 2 <= hi
        for r in rows:
            np.testing.assert_allclose(r["neutral_before_both"], 1, atol=3e-14)
            assert r["charge_preflight"] == [[0, measure.SHOTS, 0, 0]] * 2
            for bank in r["banks"]:
                sign = bank["sign"]
                idx = 2 * (sign + 1)
                for site in ("alice", "bob"):
                    counts = np.array(bank[site])
                    assert counts.sum() == counts[idx] + counts[3] == measure.SHOTS
                    lo, hi = measure.interval(counts[idx])
                    assert lo <= bank[site + "_hit"] <= hi
                np.testing.assert_array_equal(
                    bank["same_event_material_joint"], np.diag(bank["bob"])
                )
        assert (
            sum(
                b["bob"][2 * (b["sign"] + 1)] + b["bob"][3]
                for r in rows
                for b in r["banks"]
            )
            == 23068672
        )
        assert Counter(b["classification"] for r in rows for b in r["banks"]) == {
            "compatible": 16,
            "held_mismatch": 4,
            "excluded_at_alice": 2,
        }
        assert [r["blind_source"] for r in rows if not r["alice_mixture_feasible"]] == [
            6,
            7,
        ]
        assert [b["classification"] for b in rows[10]["banks"]] == ["held_mismatch"] * 2
        assert max(b["unchanged_alice_device_tv"] for b in rows[0]["banks"]) > 0.70
        for j, sign in enumerate((-1, 1)):
            for site in ("alice", "bob"):
                dwell = measure.EXPOSURES[site][sign]
                active = np.sin(dwell) ** 2
                pair = np.sin(np.sqrt(2) * dwell) ** 2 / 2
                expected = (
                    [active, 0, 0, active, 0, 0.31 * active]
                    if sign == -1
                    else [0, active, 0, 0, active, 0.46 * active]
                ) + [pair, 2 * pair, 0, 0, 0.4 * pair]
                coefficients = [r["banks"][j][site + "_hit"] for r in rows]
                np.testing.assert_allclose(coefficients, expected, atol=2e-13, rtol=0)
