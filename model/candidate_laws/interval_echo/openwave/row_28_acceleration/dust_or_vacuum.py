"""Separate a missed dust contribution from a vacuum term using moving observers."""

import json
from fractions import Fraction
from math import pi
from pathlib import Path

from model.candidate_laws.interval_echo_published_records.spherical_dust.experiment import (
    measured_frame,
)

HERE = Path(__file__).resolve().parent
SURVEY = HERE.parents[2] / "interval_echo_published_records/spherical_dust/results.json"


def compare(row, gain, include_collar):
    count = [Fraction(64 * (hi - lo), 7 * n) for n, lo, hi in row["mass_tallies"]]
    metered = sum(count) if include_collar else count[0]
    level = row["inverse_radius"]
    frames = []
    for obs in row["observers"]:
        gamma = float(measured_frame(obs["frame_phase_stamps"], level)[0])
        curvature = -15 * obs["area_over_r4"] / (4 * pi)
        frames.append((gamma**2, curvature - 6 * gain * float(metered) * gamma**2))
    # Two contractions suffice to distinguish constant vacuum from comoving dust.
    # The remaining eight frames check the interpolation on this spherical family.
    at_rest = frames[0]
    fast = max(enumerate(frames), key=lambda f: f[1][0])[0]
    g, d = frames[fast]
    slope = (d - at_rest[1]) / (g - at_rest[0])
    vacuum = at_rest[1] - slope * at_rest[0]
    return dict(
        inverse_radius=level,
        dust_coefficient=slope,
        vacuum_intercept=vacuum,
        residual_range=[min(d for g, d in frames), max(d for g, d in frames)],
        max_other_frame_error=max(
            abs(d - (vacuum + slope * g))
            for i, (g, d) in enumerate(frames)
            if i not in (0, fast)
        ),
    )


def run():
    survey = json.loads(SURVEY.read_text())
    cases = {}
    for name, stock, include_collar in (
        ("ordinary", "ordinary", True),
        ("held", "held", True),
        ("omit_counted_collar", "held", False),
        ("half_collar_force", "half_collar_force", True),
    ):
        rows = survey["ball_surveys"][stock]
        cases[name] = [
            compare(
                row,
                survey["ordinary_calibration"][str(row["inverse_radius"])]["gain"],
                include_collar,
            )
            for row in rows
        ]
    return cases


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
