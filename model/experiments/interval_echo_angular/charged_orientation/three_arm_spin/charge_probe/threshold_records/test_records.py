"""Check the detector decision rule and replay the deposited sample."""

from collections import defaultdict
from hashlib import sha256
import json

import numpy as np

from . import measure


def local_marginals(bank):
    values = bank["computed_probabilities"]
    totals = [defaultdict(float) for _ in range(3)]
    multiple = 0.0
    for key, mass in values.items():
        reads = key.split("/")
        if len(reads) != 3:
            continue
        deposits = sum(record.split(",")[0] != "0" for record in reads)
        if deposits > 1:
            multiple += mass
        for output, record in zip(totals, reads, strict=True):
            output[record] += mass
    return totals, multiple


def test_deposited_joint_records_and_independent_marginals():
    data = json.loads((measure.HERE / "results.json").read_text())
    assert (
        data["source_threshold_sha256"]
        == sha256(measure.THRESHOLD.read_bytes()).hexdigest()
    )
    assert measure.calibrated_echo()[0, 0] == (-2, 0)
    assert data["source_double_deposit_control"] == [-2, 0]
    assert [row["length"] for row in data["lengths"]] == [19, 20, 21]
    for row in data["lengths"]:
        for sign, readings in row["boundary"].items():
            d = -int(sign)
            counts = readings["synthetic_counts"]
            assert sum(counts.values()) == measure.BOUNDARY_STARTS
            assert counts.get(f"{sign},{d}/idle", 0) > 900
            assert counts.get(f"{sign},{d}/changed", 0) == 0
            assert readings["boundary_change_given_signed_tip_upper99"] < 0.01
        for event in row["middle"]:
            native, split_control = event["native"], event["independent_local_control"]
            for bank in (native, split_control):
                assert sum(bank["synthetic_counts"].values()) == measure.CUT_STARTS
            single_records, single_multiple = local_marginals(native)
            split_records, split_multiple = local_marginals(split_control)
            for a, b in zip(single_records, split_records, strict=True):
                assert a.keys() == b.keys()
                np.testing.assert_allclose(
                    list(a.values()), [b[key] for key in a], atol=1e-20, rtol=1e-12
                )
            assert single_multiple == 0.0 and split_multiple > 0.0


def test_replay_agrees_with_the_archive():
    reference = json.loads((measure.HERE / "results.json").read_text())
    current = measure.run()
    # Each joint bank has one idle outcome and up to two on each of three arms.
    comparisons = 7 * sum(len(row["middle"]) for row in reference["lengths"])
    tail = 0.01 / (2 * comparisons)
    for a, b in zip(current["lengths"], reference["lengths"], strict=True):
        assert a["length"] == b["length"] and a["cuts"] == b["cuts"]
        for sign in ("-1", "1"):
            current_read, archived_read = a["boundary"][sign], b["boundary"][sign]
            np.testing.assert_allclose(
                current_read["signed_tip_arrival_per_original_start"],
                archived_read["signed_tip_arrival_per_original_start"],
                rtol=2e-6,
                atol=0,
            )
        for now, old in zip(a["middle"], b["middle"], strict=True):
            for key in ("preparation", "time", "kappa"):
                assert now[key] == old[key]
            counts = old["native"]["synthetic_counts"]
            outcomes = {key: n for key, n in counts.items() if key.count("/") == 2}
            admitted = sum(outcomes.values())
            assert admitted > 0
            for event in (now, old):
                law = event["native"]["computed_probabilities"]
                assert law.keys() == counts.keys()
                for outcome, n in outcomes.items():
                    lo, hi = measure.bin_interval(n, admitted, tail=tail)
                    probability = law[outcome] / event["admission"]
                    assert lo <= probability <= hi, (a["length"], now["time"], outcome)
