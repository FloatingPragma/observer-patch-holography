"""Reconstruct the alternating recoil bin directly from all incident detector records."""

import json
from pathlib import Path

from model.candidate_laws.interval_echo.openwave.row_21_boost import measure as boost


def test_two_visits_from_independent_single_and_gap_records():
    source = json.loads((boost.HERE / "results.json").read_text())
    read = json.loads(Path(__file__).with_name("records.json").read_text())
    n = read["shots"]
    for law in read["studies"]:
        length = len(law["route"])
        assert law["empty"] == [n] + [0] * (length + 2)
        standard = law["single"]
        assert sum(standard) == n
        gain = sum(standard[1:length]) * length**2 / ((length - 1) * n)
        prior = next(
            s
            for s in source["studies"]
            if s["board_cells"] == 80 and s["contrast"] == law["contrast"]
        )
        for panel in law["panels"]:
            shape = next(p for p in prior["panels"] if p["winding"] == panel["winding"])
            hist = shape["separation"]["16"]["counts"]
            adjacency = hist[0] / sum(hist[:-1])
            assert sum(panel["one_contact"]) == sum(panel["two_contact"]) == n
            assert panel["one_contact"][-2] == 0
            signal = panel["two_contact"][length // 2]
            prediction = gain**2 * adjacency / length**2
            lo, hi = panel["reading"]["two_contact_pi_interval"]
            assert lo <= prediction <= hi
            if law["contrast"] == 8:
                assert signal > 0 and lo > 0
            else:
                assert signal == 0
