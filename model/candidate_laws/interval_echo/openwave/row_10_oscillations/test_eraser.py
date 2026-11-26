"""Check the ideal source eraser against passive remote counts and controls."""

import json
from pathlib import Path

import numpy as np
from scipy.stats import binom

from model.candidate_laws.interval_echo.openwave.row_10_oscillations.eraser import (
    forecast,
    read,
)


def test_sources_with_the_same_sharp_populations_can_give_opposite_eraser_exits():
    ket = np.array([0, 1, 1j, 0]) / np.sqrt(2)
    matrix = np.outer(ket, ket.conj())
    np.testing.assert_allclose(read(matrix)["Y"], [0, 1, 0, 0], atol=1e-13)
    np.testing.assert_allclose(read(matrix.conj())["Y"], [0, 0, 1, 0], atol=1e-13)
    np.testing.assert_allclose(
        read(np.diag(matrix.diagonal()))["Y"], [0, 0.5, 0.5, 0], atol=1e-13
    )


def test_forecast_preserves_passive_receivers_and_original_start_denominator():
    records = json.loads(Path(__file__).with_name("records.json").read_text())
    banks = {tuple(r["key"]): r["counts"] for r in records["banks"]}
    result = forecast()
    # Familywise allowance over both tails of 54 independent passive banks.
    per_tail = 0.01 / (2 * 54)
    comparisons = 0
    for held in result["held"]:
        for outcome in held["outcomes"]:
            pair = sum(outcome["granted_joint"]["X"])
            np.testing.assert_allclose(sum(outcome["granted_joint"]["Y"]), pair, atol=0)
            total = (
                pair
                + outcome["far_pair_missed"]
                + outcome["initial_arrival_failed"]
                + outcome["original_refused"]
            )
            assert abs(total - 1) < 1e-14
            prefix = (
                8,
                False,
                held["owner"],
                held["time"],
                "character",
                outcome["recipe"],
            )
            for panel in ("sharp", "left", "right"):
                counts = banks[(*prefix, panel)]
                detected, starts = sum(counts[:6]), sum(counts)
                assert binom.cdf(detected, starts, pair) > per_tail
                assert binom.sf(detected - 1, starts, pair) > per_tail
                comparisons += 1
    assert comparisons == 54
