"""Bound the effect of stellar calibration and arbitrary baryonic-force errors."""

import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
KPC_M = F("3.0856775814913673e19")
DEEP = F("1.2e-11")
DISK = (F(".4"), F(".6"))
BULGE = (F(".56"), F(".84"))
MATCH = F(11, 10)


def measured(row):
    """Return radius, orbital speed, speed error and signed force components."""
    radius, speed, error, gas, disk, bulge = [
        F(row[a:b])
        for a, b in [(19, 25), (26, 32), (33, 38), (39, 45), (46, 52), (53, 59)]
    ]
    components = [v * abs(v) * 10**6 / (radius * KPC_M) for v in (gas, disk, bulge)]
    return radius, speed, error, components


def force(components, disk, bulge):
    gas, star, centre = components
    return gas + disk * star + bulge * centre


def audit_pair(pair, lines):
    faint, bright = [
        measured(lines[pair[k]["table2_line"] - 1]) for k in ("faint", "bright")
    ]
    rf, vf, ef, cf = faint
    rb, vb, eb, cb = bright
    assert vb > 2 * eb
    high = (rb / rf) * ((vf + 2 * ef) / (vb - 2 * eb)) ** 2
    # For positive forces, min(t,sqrt(t)) <= H iff t <= max(H,H^2).
    limit = max(high, high**2)
    nominal = force(cf, F(".5"), F(".7")) / force(cb, F(".5"), F(".7"))
    assert nominal > limit
    calibrations = list(product(DISK, BULGE))
    corners = [(force(cf, d, b), force(cb, d, b)) for d, b in calibrations]
    assert all(min(f, h) > 0 for f, h in corners)
    ratios = [f / h for f, h in corners]
    return dict(
        galaxy=pair["galaxy"],
        faint_table_line=pair["faint"]["table2_line"],
        bright_table_line=pair["bright"]["table2_line"],
        observed_upper=float(high),
        nominal_baryonic_ratio=float(nominal),
        largest_permitted_baryonic_ratio=float(limit),
        minimum_independent_symmetric_fraction=float(
            (nominal - limit) / (nominal + limit)
        ),
        stellar_rectangle_baryonic_ratio=[float(min(ratios)), float(max(ratios))],
        remains_below_through_stellar_rectangle=min(ratios) > limit,
        remains_below_for_any_ten_percent_force_match=1 / MATCH > limit,
        eligible_corners_reaching_shape_bound=[
            [float(d), float(b)]
            for (d, b), (f, h) in zip(calibrations, corners)
            if 1 / MATCH <= f / h <= MATCH and f / h <= limit and max(f, h) < DEEP
        ],
        always_satisfies_pair_force_match=min(ratios) >= 1 / MATCH
        and max(ratios) <= MATCH,
        always_satisfies_deep_cut=all(max(f, h) < DEEP for f, h in corners),
    )


def run():
    pairs = json.loads((HERE / "results.json").read_text())["cuts"]["deep"]["pairs"]
    lines = (HERE.parent / "data/table2.dat").read_text().splitlines()
    return [
        audit_pair(p, lines)
        for p in pairs
        if p["speed_error_boxes"]["2"]["sign"] == "smaller_faint"
    ]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    receipt = HERE / "sensitivity.json"
    result = run()
    if parser.parse_args().write:
        receipt.write_text(json.dumps(result, indent=2) + "\n")
    else:
        assert json.loads(receipt.read_text()) == result
        print("Replayed the stellar mass and force-ratio sensitivity.")
