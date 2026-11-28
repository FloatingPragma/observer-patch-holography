"""Recompute the row 30 dispersion and material-choice contrasts from records."""

import json
from itertools import product
from math import isclose, isfinite, sin, sqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
EVIDENCE = ROOT / "model/candidate_laws/interval_echo/openwave"


def read(name):
    return json.loads((EVIDENCE / name).read_text(encoding="utf-8"))


def held_dispersion():
    base = "submission/row_30_confirmation"
    v3 = read(f"{base}/validation.json")
    v2 = read(f"{base}/v2/validation.json")
    if v3["status"] != "passed" or v2["status"] != "rejected":
        raise ValueError("Archived decisions differ")
    errors = []
    for board in v3["boards"]:
        spatial = board["spatial"]
        gate = spatial["quadratic_cutoff"]
        if not (
            spatial["records_complete"]
            and spatial["candidate_q"] < gate < spatial["launch_label_q"]
        ):
            raise ValueError(
                "Held spatial profiles do not discriminate launch wavelengths"
            )
        meter = board["dispersion"]["meter"]
        h, j = meter["h"], meter["j"]
        for mode, k in zip(
            board["dispersion"]["modes"], spatial["even_trained_k"], strict=True
        ):
            if k != mode["k"] or not all(isfinite(v) for v in (k, mode["omega"])):
                raise ValueError("Dispersion wavelength differs from held spatial fit")
            lattice = sqrt(4 * (j - h) ** 2 + 16 * j * h * sin(k / 2) ** 2)
            continuum = sqrt(4 * (j - h) ** 2 + 4 * j * h * k**2)
            if not isclose(lattice, mode["lattice"], rel_tol=0, abs_tol=1e-12):
                raise ValueError("Independent meter does not reproduce the forecast")
            errors.append(
                (k, abs(mode["omega"] - lattice), abs(mode["omega"] - continuum))
            )
    maximum = max(a for _, a, _ in errors)
    low_count = sum(k <= 0.5 for k, _, _ in errors)
    high_error = max(c for k, _, c in errors if k > 0.5)
    if len(errors) != 22 or maximum >= 0.000575 or high_error <= 0.0048:
        raise ValueError("Held dispersion comparison differs")
    print(
        f"V3: {len(errors)} modes, lattice max error {maximum:.9f}, "
        f"{low_count} long waves, shorter-wave continuum max error {high_error:.6f}"
    )


def material_choice():
    profiles = read("row_30_collective/sharp_growth/conditional_profiles.json")
    rows = profiles["rows"]
    links = len(profiles["route"]) - 1
    patterns = {"".join(bits) for bits in product("01", repeat=links)}
    if {row["record"] for row in rows} != patterns or len(rows) != len(patterns):
        raise ValueError("Incomplete set of sharp write histories")
    passed = set()
    for row in rows:
        statistic, gate = row["odd_profile_q"], row["quadratic_cutoff"]
        if not (isfinite(statistic) and isfinite(gate)):
            raise ValueError("Nonfinite first-write comparison")
        accepted = statistic < gate
        if accepted != row["profile_pass"]:
            raise ValueError("Archived profile decision differs from its gate")
        if accepted:
            passed.add(row["record"])
    if passed != {"0" * links, "1" * links}:
        raise ValueError("Material-choice contrast differs")
    print(f"Complete first-write census: {len(passed)}/{len(rows)} pass")
    for dimers in (14, 20, 36):
        print(
            f"{dimers} dimers, one fair write per connector: "
            f"uniform fraction 2^-{dimers - 2}"
        )


def polar_curvature():
    panels = read("row_30_collective/beat_speed/results.json")["polar"]
    if len(panels) != 8:
        raise ValueError("Need both ring sizes, strengths and directions")
    increase, decrease = [], []
    for row in panels.values():
        speed = row["speed"]
        if not all(0 <= lo <= hi and isfinite(hi) for lo, hi in speed):
            raise ValueError("Invalid polar speed interval")
        increase.append(speed[1][0] - speed[0][1])
        decrease.append(speed[-2][0] - speed[-1][1])
    if min(increase) <= 0 or min(decrease) <= 0:
        raise ValueError("No curvature-sign reversal")
    print(
        f"Polar rings: all {len(panels)} momentum rays reverse secant order; "
        f"minimum increase {min(increase):.6f}, minimum decrease {min(decrease):.6f}"
    )


if __name__ == "__main__":
    held_dispersion()
    material_choice()
    polar_curvature()
