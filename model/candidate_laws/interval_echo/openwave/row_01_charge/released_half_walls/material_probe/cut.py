"""Whole-seam shutoff versus moving the two probes to the distant zero stock."""

import argparse
import json
from itertools import combinations, product
from pathlib import Path

from . import live
from .measure import HERE

MODES = ("live", "detached", "cut")
SEED = 26011264


def decode(raw):
    rows = live.load_banks(raw, MODES, SEED)
    report = {}
    for k in (8, 4):
        tally = {}
        for first, second in combinations(MODES, 2):
            wall = [0, 0]
            probes = [0, 0]
            for src, s, i, r in product(
                (-1, 1), (-1, 1), range(len(live.HOLDS)), live.ANALYZERS
            ):
                a, b = (rows[k, src, s, mode, i, r][1] for mode in (first, second))
                if r == "rail":
                    wall[0] += live.separate(
                        live.bound(a.sum(), live.LIVE_SHOTS),
                        live.bound(b.sum(), live.LIVE_SHOTS),
                    )
                    wall[1] += 1
                for lhits, rhits in (
                    (a[1].sum(), b[1].sum()),
                    (a[:, 1].sum(), b[:, 1].sum()),
                ):
                    probes[0] += live.separate(
                        live.bound(lhits, a.sum()), live.bound(rhits, b.sum())
                    )
                    probes[1] += 1
            tally[first + "_" + second] = dict(wall_rate=wall, conditional_probe=probes)
        example = {mode: int(rows[k, 1, 1, mode, 1, "rail"][1].sum()) for mode in MODES}
        report[k] = dict(comparisons=tally, wall_count_positive_at_quarter=example)
    return report


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    if args.output:
        raw = live.acquire(MODES, SEED)
        args.output.write_text(json.dumps(raw, separators=(",", ":")) + "\n")
    else:
        raw = json.loads((HERE / "cut_counts.json").read_text())
    print(json.dumps(decode(raw), indent=2))
