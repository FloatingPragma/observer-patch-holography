import json

import numpy as np
import pytest
import sympy as sp

from . import observer_bound as ob


def test_matching_projectors_with_an_additional_leakage_output():
    # Two reference coordinates, three output coordinates (one is outside the code).
    i = sp.eye(2)
    embed = sp.eye(3)[:, :2]
    z0, z1 = i[:, 0], i[:, 1]

    def matched(v):
        w = sp.kronecker_product(v, embed * v)
        return w * w.T

    pairs = [sp.kronecker_product(z, embed * z) for z in (z0, z1)]
    phi = sum(pairs, sp.zeros(6, 1)) / sp.sqrt(2)
    x_match = sum((matched((z0 + s*z1) / sp.sqrt(2)) for s in (-1, 1)), sp.zeros(6))
    residual = sp.eye(6) + phi*phi.T - matched(z0) - matched(z1) - x_match
    assert residual == residual.T and residual**2 == residual and residual.trace() == 3


def test_changed_clicks_weaken_the_bound_and_malformed_records_fail():
    receipt = json.loads(ob.RECEIPT.read_bytes())
    banks = receipt["fresh_calibration_banks"] + [
        r for r in receipt["new_banks"] if r["probe"] == "minus" and r["instrument"] != "pair"]
    baseline = ob.infer(banks, ob.ATTEMPTS)
    assert np.isclose(baseline[0]["any_event_error_upper"], 0.09095837094712317)
    banks[0]["nonzero_joint_counts"] = [[0, 0, 8182], [0, 1, 10]]
    assert ob.infer(banks, ob.ATTEMPTS)[0]["any_event_error_upper"] > baseline[0]["any_event_error_upper"]
    banks[0]["nonzero_joint_counts"].pop()
    with pytest.raises(ValueError, match="incomplete raw cohort"):
        ob.infer(banks, ob.ATTEMPTS)


def test_native_replay_and_all_held_decisions():
    fresh = ob.report()
    ob.compare(json.loads(ob.RECEIPT.read_bytes()), fresh)
    for mode, passed in (("pair", True), ("address_then_forget", False)):
        panels = [b for b in fresh["held"] if b["instrument"] == mode]
        assert len(panels) == 12 and all(b["within_coarse_bound"] is passed for b in panels)


def test_every_held_interval_separates_the_readers():
    receipt = json.loads(ob.RECEIPT.read_bytes())
    radius = receipt["held_sampling_allowance"]
    for panel in receipt["held"]:
        limit = receipt["bounds"][panel["pair_start"] // 2]["any_event_error_upper"]
        error = abs(panel["rate"] - panel["forecast"])
        if panel["instrument"] == "pair":
            assert error + radius < limit
        else:
            assert error - radius > limit


def test_sampled_grid_matches_analytic_raw_law_in_all_72_settings():
    from decimal import Decimal, localcontext

    from ..sequential import GRID, grid_weights

    calls = ob.pair.address.schedule()
    writers = {
        "pair": ob.pair.read_branches([c for c in calls if c["bit"] != 0], 0),
        "address_then_forget": ob.pair.read_branches(calls, 0),
    }
    terminal = {r: ob.pair.read_branches(calls, r) for r in range(12)}
    # At the held angles even the irrational bright probabilities pick a
    # definite nearest grid integer, far from a half-integer tie.
    with localcontext() as context:
        context.prec = 50
        coarse_bright = {
            "plus": GRID,
            "minus": 0,
            "quarter": round(GRID * (2 + Decimal(2).sqrt()) / 4),
            "three_quarters": round(GRID * (2 - Decimal(2).sqrt()) / 4),
        }

    for start in range(0, 12, 2):
        for writer, head in writers.items():
            for probe in ("port0", "port1", *coarse_bright):
                if probe in ("port0", "port1", "plus"):
                    joint, _ = ob.pair.read_probe(head, terminal, start, probe)
                else:
                    joint = ob.native_joint(head, terminal, start, probe)
                sampled = grid_weights((joint / joint.sum()).ravel()).reshape(17, 17)
                expected = np.zeros((17, 17), dtype=np.int64)
                if probe in ("port0", "port1"):
                    port = int(probe == "port1")
                    first = start + port if writer == "address_then_forget" else start
                    expected[first, first ^ port] = GRID
                elif writer == "address_then_forget":
                    # Both first coordinates are equally probable. Each then
                    # gives either decoded terminal address with probability 1/2.
                    for first in (start, start + 1):
                        expected[first, first] = GRID // 4
                        expected[first, first ^ 1] = GRID // 4
                else:
                    expected[start, start] = coarse_bright[probe]
                    expected[start, start ^ 1] = GRID - coarse_bright[probe]
                np.testing.assert_array_equal(
                    sampled, expected, err_msg=f"pair {start}, {writer}, {probe}"
                )
