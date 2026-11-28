"""Bound ordinary-calibrated source transfer for every common storage fraction."""

import argparse
import json

from mpmath import iv

from . import (
    bridge_radar as bridge,
    paired_relay as pair,
    paired_survey as radar,
    stored_source as store,
)

RECEIPT = bridge.HERE / "continuous_storage.json"
BINS = 32


def blend(start, end, fraction):
    """Enclose the affine interpolation by interpolation along its lower/upper bounds."""
    lows = [start.a + f * (end.a - start.a) for f in (fraction.a, fraction.b)]
    highs = [start.b + f * (end.b - start.b) for f in (fraction.a, fraction.b)]
    return iv.mpf([min(x.a for x in lows), max(x.b for x in highs)])


def factors(sides):
    a, b, c = sides
    return (a + b + c, -a + b + c, a - b + c, a + b - c)


def triangle(bank):
    lengths = bridge.ranges(bank, store.L)
    return [(lengths[i, j] + lengths[j, i]) / 2 for i, j in ((0, 1), (1, 2), (0, 2))]


def difference_area(plus, minus, delta):
    # The supplied four factor differences may all be divided by the same f.
    # Telescoping the Heron product takes that f outside the sum exactly.
    terms = iv.mpf(0)
    xp, xm = iv.mpf(1), iv.mpf(1)
    for j, d in enumerate(delta):
        term = d
        for x in plus[:j] + minus[j + 1 :]:
            term *= x
        terms += term
        xp *= plus[j]
        xm *= minus[j]
    return terms / (4 * (iv.sqrt(xp) + iv.sqrt(xm)))


def analyze(data, arm):
    steps = data["steps"]
    caps = []
    for key in ("all_accept", "all_fail"):
        banks = pair.checked(data["panels"][key]["radar"][arm], (4, 6, 11), store.SHOTS)
        caps.append([triangle(b) for b in banks])
    # Products of the two held stocks coincide exactly in the complete-storage
    # arm. Both samples bound the same side means; either arm suffices here.
    caps[0][3] = caps[0][2]
    src0 = data["measurements"]["all_accept"]
    src1 = data["measurements"]["all_fail"]
    dh = iv.mpf(src1["held"])
    held_delta = factors([a - b for a, b in zip(caps[1][3], caps[1][2], strict=True)])
    ordinary_delta = [
        factors([a - b for a, b in zip(z[1], z[0], strict=True)]) for z in caps
    ]
    residuals = []
    for i in range(BINS):
        box = iv.mpf([i / BINS, (i + 1) / BINS])
        stocks = [
            factors([blend(x, y, box) for x, y in zip(a, b, strict=True)])
            for a, b in zip(*caps, strict=True)
        ]
        dm = [blend(x, y, box) for x, y in zip(*ordinary_delta, strict=True)]
        geom_o = difference_area(stocks[1], stocks[0], dm) / (2 * steps["ordinary"])
        geom_h_over_f = difference_area(stocks[3], stocks[2], held_delta) / (
            2 * steps["held"]
        )
        src_o = blend(iv.mpf(src0["ordinary"]), iv.mpf(src1["ordinary"]), box)
        residuals.append(radar.limits(geom_h_over_f * src_o - geom_o * dh))
    lo = min(x[0] for x in residuals)
    hi = max(x[1] for x in residuals)
    return dict(
        boxes=residuals, lower=lo, upper=hi, all_fractions_excluded=(lo > 0 or hi < 0)
    )


def run():
    data = json.loads(store.RECEIPT.read_text())
    if store.decode(data) != data["measurements"]:
        raise ValueError("Inconsistent source or radar bank")
    return dict(
        source=bridge.digest(store.RECEIPT),
        bins=BINS,
        results={
            arm: analyze(data, arm) for arm in data["panels"]["all_fail"]["radar"]
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = run()
    if args.write:
        with RECEIPT.open("x") as f:
            json.dump(data, f, separators=(",", ":"), allow_nan=False)
            f.write("\n")
    elif data != json.loads(RECEIPT.read_text()):
        raise ValueError("Continuous bounds differ")
    for arm, r in data["results"].items():
        print(arm, r["lower"], r["upper"], r["all_fractions_excluded"])
