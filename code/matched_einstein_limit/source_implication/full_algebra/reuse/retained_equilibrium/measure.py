"""Finite stationary return records, terminal feedback and a held area prediction."""

import argparse
import json
from pathlib import Path
import numpy as np

from . import control, gates, held

HERE = Path(__file__).parent


def successive(channels, source):
    # Preselect the same reversal for an independent sham pair; otherwise
    # conservation of the hidden parity alone would produce a correlation.
    weights = np.trace(source, axis1=1, axis2=2).real
    ch = channels.reshape(2, 24, 2, 16, 16)
    fast = ch[:, :, 1].sum(axis=1)
    p2 = 0.0
    sham = 0.0
    for z, f, w in zip(source, fast, weights, strict=True):
        first = f @ z.ravel()
        p2 += np.trace((f @ first).reshape(4, 4)).real
        sham += np.trace(first.reshape(4, 4)).real ** 2 / w
    rng = np.random.default_rng(8616229)
    shots = 2**42
    radius = np.sqrt(np.log(400) / (2 * shots)) + 2**-35
    # The two independent one-bit counts use the evaluated physical probabilities.
    draws = [int(rng.binomial(shots, round(p * 2**34) / 2**34)) for p in (p2, sham)]
    return dict(
        two_consecutive=p2,
        two_independent_at_matched_parity=sham,
        contrast=p2 - sham,
        shots_each=shots,
        simulated_early_early=draws,
        contrast_band_99=[
            (draws[0] - draws[1]) / shots + x * 2 * radius for x in (-1, 1)
        ],
    )


def run():
    x, ch, rows = control.compare()
    # The input uncertainty and both grouped integer channels are covered by
    # this total margin, including certified eigenvalue residuals (README).
    assert min(min(row["with_bit"]) for row in rows) > 0.00639
    assert max(row["without_bit"][1] for row in rows) < -0.00465
    return dict(
        schedule="12 directed routes and two charges with equal probability",
        flag="completed before tick 12, including every other outcome in the complement",
        certified_absolute_error_per_partial_sum_slack=1.3e-5,
        twirl=control.TWIRL,
        blocks=rows,
        physical_bus=gates.shuttle(),
        repeat=successive(ch, x),
        held=held.read(ch, x),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    result = run()
    if args.write:
        (HERE / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    else:
        print(json.dumps(result, indent=2))
