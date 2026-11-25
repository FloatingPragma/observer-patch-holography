"""Odd source changes on one set of stationary material guides."""

import argparse
import json
from pathlib import Path

import numpy as np

from . import measure

HERE = Path(__file__).parent
SITES = ((6, 6, 6), (5, 6, 6), (7, 7, 7), (4, 7, 7))


def source(b, s):
    load = np.zeros((10, 10, 10))
    for v in measure.CORNERS:
        ix, iy, iz = (2 * x - 11 for x in v)
        load[tuple(x - 1 for x in v)] = (1 + b * ix + s * ix * iy * iz) / 8
    if abs(b) + abs(s) > 1:
        raise ValueError("negative source probability")
    return load


def field(q, b, s):
    # In the isotropic case Phi=q*expected occupation and has unit source.
    return float(q[0]) * measure.poisson(source(b, s), q)


def run():
    steps, q = measure.measured_steps(None)
    central = field(q, 0, 0)
    dipole = field(q, 1, 0) - central
    octupole = field(q, 0, 1) - central
    # Every axis reflection changes the parity sign. Each spherical radial
    # bin can be checked without pairing any two noisy observations.
    for axis in range(3):
        np.testing.assert_allclose(
            np.flip(octupole, axis), -octupole, atol=1e-15, rtol=0
        )
        sign = -1 if axis == 0 else 1
        np.testing.assert_allclose(
            np.flip(dipole, axis), sign * dipole, atol=1e-15, rtol=0
        )
    single = []
    for site in SITES:
        point = np.zeros((10, 10, 10))
        point[tuple(v - 1 for v in site)] = 1
        single.append(float(q[0] * measure.poisson(point, q)[site]))
    panels = []
    for bias, seed in ((-0.5, 930825), (0.5, 930826)):
        phi = field(q, 0, bias)
        forecast = [float(phi[site] / single[j]) for j, site in enumerate(SITES)]
        counts, ticks = measure.acquire(steps, seed, SITES, bias)
        panels.append(
            dict(
                bias=bias,
                seed=seed,
                hits=counts[5:],
                forecasts=forecast,
                exit_axis_counts=counts[:3],
                last_arrival=ticks,
            )
        )
    return dict(
        shots_per_panel=measure.SHOTS,
        sites=[list(s) for s in SITES],
        dipole_potential_slope=[float(dipole[s]) for s in SITES],
        octupole_potential_slope=[float(octupole[s]) for s in SITES],
        simultaneous_99_radius=float(
            np.sqrt(np.log(2 * 14 / 0.01) / (2 * measure.SHOTS))
        ),
        panels=panels,
    )


def verify(actual, saved):
    for key in ("shots_per_panel", "sites"):
        if actual[key] != saved[key]:
            raise ValueError("parity observation metadata differs")
    for key in (
        "dipole_potential_slope",
        "octupole_potential_slope",
        "simultaneous_99_radius",
    ):
        np.testing.assert_allclose(actual[key], saved[key], atol=1e-14, rtol=0)
    for fresh, old in zip(actual["panels"], saved["panels"], strict=True):
        for key in ("bias", "seed", "hits", "exit_axis_counts", "last_arrival"):
            if fresh[key] != old[key]:
                raise ValueError(f"parity detector field {key} differs")
        np.testing.assert_allclose(
            fresh["forecasts"], old["forecasts"], atol=2e-14, rtol=0
        )
        expected = np.array([*fresh["forecasts"], *([1 / 3] * 3)])
        observed = np.array([*old["hits"], *old["exit_axis_counts"]]) / measure.SHOTS
        if max(abs(observed - expected)) > actual["simultaneous_99_radius"]:
            raise ValueError("held parity forecast lies outside the simultaneous band")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "multipoles.json"
    if args.write:
        verify(result, result)
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        verify(result, json.loads(path.read_text()))
        print("Replayed both parity-source acquisitions and their angular forecasts.")
