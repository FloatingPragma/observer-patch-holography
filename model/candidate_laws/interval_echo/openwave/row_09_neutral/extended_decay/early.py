"""Validate the early time-of-flight window using the independent source."""

import argparse
import json
from pathlib import Path
import numpy as np
from . import cold_ports
from .measure import HERE
from ..passive_decay.measure import compare_records

WAITS = tuple(3 + (j + 0.5) / 4 for j in range(8))


def verdict(data, require_each=False):
    settings = data["mounts"]["complete"]
    rows = settings["native"]["fixed"]
    p = settings["native"]["cold_preparation"]["admission"]

    def diff(arms):
        return (
            np.asarray(arms["Z"]["panels"])[:, 2]
            - np.asarray(arms["idle"]["panels"])[:, 2]
        )

    cut = diff(settings["exit_cut"]["fixed"])
    np.testing.assert_allclose(cut, 0, atol=3e-12)
    delta = p * diff(rows)
    low, heavy = rows["Z"]["unresolved_loss_light_gain_heavy_gain"][2:]
    margin = (min(delta) if require_each else np.mean(delta)) - p * data[
        "clock_false_light"
    ]
    print(
        "original-launch excesses",
        delta,
        "false-light adjusted qualification margin",
        margin,
    )
    return bool(margin > 0.01 and rows["Z"]["work"] > 0 and low > heavy)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    parser.add_argument("--ticks", action="store_true")
    args = parser.parse_args()
    if args.ticks:
        spacing = 2 * np.pi / cold_ports.outlet.PERIOD
        center = round(4 / spacing)
        waits = tuple(float(i * spacing) for i in (center - 1, center, center + 1))
    else:
        waits = WAITS
    result = cold_ports.acquire(
        "TICKS.md" if args.ticks else "EARLY.md", waits, ("complete",)
    )
    print("qualified", verdict(result, args.ticks))
    if args.acquire:
        with args.acquire.open("x") as f:
            json.dump(result, f, allow_nan=False, separators=(",", ":"))
            f.write("\n")
    else:
        compare_records(
            json.loads(
                (HERE / ("ticks.json" if args.ticks else "early.json")).read_text()
            ),
            result,
        )
