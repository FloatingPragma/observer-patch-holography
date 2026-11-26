"""Acquire sharp conversion-site and distal pair records on the native release."""

import argparse
import hashlib
import json
from itertools import product
from pathlib import Path

import numpy as np

from .. import eraser
from ..measure import OWNERS, TIMES, sector_fields
from .apparatus import Terminal, access, pointer, recipe

SEED = 26112641
STARTS = 2**52
BENCH = 2**22
MODES = ("native", "which_head", "skip_crossing", "contrast_four")


def bins(effect, density):
    rates = np.einsum("ijk,kj->i", effect, density).real
    if rates.min() < -2e-12 * np.trace(density).real:
        raise ValueError("negative material event")
    return np.maximum(rates, 0)


def sample(key, rates, n):
    p = np.asarray(rates, float)
    if p.min() < -2e-14 or abs(sum(p) - 1) > 2e-12:
        raise ValueError("incomplete instrument")
    p = np.maximum(p, 0)
    digest = hashlib.sha256(json.dumps([SEED, *key]).encode()).digest()
    rng = np.random.default_rng(int.from_bytes(digest[:16], "big"))
    order = np.argsort(p, kind="stable")
    return rng.multinomial(n, p[order] / sum(p))[np.argsort(order)].tolist()


def standards():
    lab = Terminal("tip")
    sharp = list(np.eye(4, dtype=complex))
    prepared = [lab.run(recipe(s))[0, :4] for s in (-1, 1)]
    rows = []
    for sense, omit, (j, col) in product(
        (-1, 1), (False, True), enumerate(sharp + prepared)
    ):
        # The charged converter and actual addressed stock preparation preserve
        # the two head indices on the selected branch (see access()).
        phase = complex(*access()["conversion_amplitude"])
        col = col * phase
        p = bins(pointer(sense, omit), np.outer(col, col.conj()))
        key = [sense, omit, j]
        rows.append(
            dict(key=key, counts=sample(["bench", *key], [*p, 1 - sum(p)], BENCH))
        )
    return rows


def acquire():
    apparatus = access()
    bench = standards()
    source = eraser.launch()
    held = []
    for k in (8, 4):
        waves = sector_fields(k, False)
        for owner, (t, time), sign, panel in product(
            OWNERS, enumerate(TIMES), (1, 2), (-1, 1)
        ):
            field = waves[owner][:, t, :]
            joint = source[sign] * (field @ field.conj().T)
            far = np.trace(joint).real
            arrived = np.trace(source[sign]).real
            for mode in MODES[:3] if k == 8 else MODES[3:]:
                # A local sharp first-bit read (vacancy flag test) has both
                # outcomes kept; its unread result removes the two cross blocks.
                if mode == "which_head":
                    same = (np.arange(4)[:, None] % 2) == (np.arange(4)[None, :] % 2)
                    measured = joint * same
                else:
                    measured = joint
                p = bins(pointer(panel, mode == "skip_crossing"), measured)
                np.testing.assert_allclose(sum(p), far, atol=5e-24, rtol=1e-12)
                # Full donor conversion is exact. No additional failure branch
                # has support, but its count is explicitly retained.
                key = [owner, time, sign, panel, mode]
                counts = sample(
                    key, [*p, arrived - far, 0, 1 / 3 - arrived, 2 / 3], STARTS
                )
                held.append(dict(key=key, counts=counts))
    return dict(
        seed=SEED,
        starts=STARTS,
        bench_starts=BENCH,
        bench=bench,
        held=held,
        source_apparatus=apparatus,
        bench_order=["success", "wrong flags", "away from hub", "source refusal"],
        held_order=[
            "success",
            "wrong flags",
            "away from hub",
            "pair missed",
            "donor failed",
            "initial charged arrival missed",
            "source refusal",
        ],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path, required=True)
    path = parser.parse_args().acquire
    path.write_text(json.dumps(acquire(), indent=1, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
