from math import comb, log
import json

from . import measure


def test_kl_inversion_has_finite_binomial_coverage_including_boundary_counts():
    for n in (5, 17, 31):
        windows = [measure.interval(k, n, log(40)) for k in range(n + 1)]
        for p in (0, 0.001, 0.1, 0.5, 0.9, 0.999, 1):
            miss = sum(
                comb(n, k) * p**k * (1 - p) ** (n - k)
                for k, (lo, hi) in enumerate(windows)
                if not lo <= p <= hi
            )
            assert miss <= 0.05 + 1e-13


def test_count_only_replay_uses_every_profile_and_discriminates_controls():
    raw = json.loads(measure.SOURCE.read_text())
    sm = measure.read(raw)
    rows = sm["regions"]
    assert len(rows) == 540
    whole = sm["full_band_at_its_launch_cell"]
    assert whole["binomial_parameters"] == 18
    assert [r["owner"] for r in whole["regions"]] == [0, 5, 10]
    assert all(r["region_charge_interval"][1] < 0.000955 for r in whole["regions"])
    for row in rows:
        if row["control"] == "native_arrival":
            assert row["interval"][1] < 0.001031
        elif row["radius"] == 0 and row["center"] == row["owner"]:
            assert row["interval"][0] > 0.27
    for s in raw["charged_standards"]:
        s.pop("probabilities", None)
    for launch in raw["launches"]:
        for panels in launch["panels"].values():
            for panel in panels:
                for s in panel["complete_records"]:
                    del s["probabilities"]
                for k in list(panel):
                    if k != "complete_records":
                        del panel[k]
    assert measure.read(raw) == sm
