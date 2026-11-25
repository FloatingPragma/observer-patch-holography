"""Predict the terminal records of the three material guides from their step law."""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.fft import dstn, idstn

from .. import measure as rails

HERE = Path(__file__).parent
SHOTS = 65536
SITES = ((8, 5, 5), (5, 8, 5))
DIRECTIONS = ((-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1))
CORNERS = tuple(product((5, 6), repeat=3))


def measured_steps(fast_axis):
    _, paths, _ = rails.layout()
    extra = (
        ((0, 0, 0),) * 2
        if fast_axis is None
        else DIRECTIONS[2 * fast_axis : 2 * fast_axis + 2]
    )
    measured = []
    for delta in (*DIRECTIONS, *extra):
        observations = []
        for axis, change in enumerate(delta):
            if change == 0:
                observations.append(0)
                continue
            path = paths[axis]
            # The complete single-guide chart uses the same pulse at every
            # active link. Test both directions on all eleven links.
            for j in range(11):
                source, dest = (j, j + 1) if change > 0 else (j + 1, j)
                observed = rails.pulse((path[j], path[j + 1]), path[source])
                if observed != path[dest]:
                    raise ValueError("step address did not transfer")
            observations.append(change)
        measured.append(observations)
    tag_steps = np.array(measured, dtype=int)
    if np.any(tag_steps.sum(axis=0)) or np.any(abs(tag_steps).sum(axis=1) > 1):
        raise ValueError("Green reader requires a centered, one-guide tag")
    covariance = tag_steps.T @ tag_steps / len(tag_steps)
    if not np.array_equal(covariance, np.diag(np.diag(covariance))):
        raise ValueError("guide steps have an off-axis covariance")
    return tag_steps, np.diag(covariance) / 2


def poisson(load, conductance):
    """Expected visits per release before a boundary stop; one unit per appointment."""
    n = load.shape[0]
    lap = 2 * (1 - np.cos(np.pi * np.arange(1, n + 1) / (n + 1)))
    eigen = (
        conductance[0] * lap[:, None, None]
        + conductance[1] * lap[None, :, None]
        + conductance[2] * lap[None, None, :]
    )
    visits = idstn(dstn(load, type=1, norm="ortho") / eigen, type=1, norm="ortho")
    padded = np.pad(visits, 1)
    rhs = np.zeros_like(load)
    for axis, q in enumerate(conductance):
        rhs -= (
            q
            * np.diff(padded, n=2, axis=axis)[
                tuple(slice(None) if j == axis else slice(1, -1) for j in range(3))
            ]
        )
    if np.max(np.abs(rhs - load)) > 2e-15:
        raise ValueError("Dirichlet Green function failed a vertex balance")
    return padded


def predictions(conductance):
    injection = np.zeros((10, 10, 10))
    for vertex in CORNERS:
        injection[tuple(v - 1 for v in vertex)] = 1 / 8
    visits = poisson(injection, conductance)
    predictions = []
    for axis in (0, 1, 2):
        slices = np.moveaxis(visits, axis, 0)
        predictions.append(
            float(conductance[axis] * (slices[1].sum() + slices[-2].sum()))
        )
    # A release always starts below both planes. The final-side indicator
    # telescopes the signed crossings on the corresponding full section.
    for axis in (0, 1):
        slices = np.moveaxis(visits, axis, 0)
        predictions.append(float(conductance[axis] * (slices[7] - slices[8]).sum()))
    for site in SITES:
        point = np.zeros_like(injection)
        point[tuple(v - 1 for v in site)] = 1
        predictions.append(float(visits[site] / poisson(point, conductance)[site]))
    return predictions


def acquire(steps, seed, sites=SITES, parity_bias=0):
    rng = np.random.default_rng(seed)
    corners = np.array([rails.release_corner(v) for v in CORNERS])
    if parity_bias:
        # First two guide bits are fair. The third fixes a parity drawn
        # with probability 1/4 or 3/4, using two further fair bits.
        if parity_bias not in (-0.5, 0.5):
            raise ValueError("the instrument has two calibrated parity settings")
        pair = rng.integers(4, size=SHOTS)
        positive = rng.integers(4, size=SHOTS) < (3 if parity_bias > 0 else 1)
        ix, iy = pair // 2, pair % 2
        iz = positive.astype(int) ^ ix ^ iy
        positions = corners[4 * ix + 2 * iy + iz]
    else:
        positions = corners[rng.integers(8, size=SHOTS)]
    active = np.arange(SHOTS)
    hits = np.column_stack([np.all(positions == site, axis=1) for site in sites])
    ticks = 0
    while len(active):
        ticks += 1
        positions[active] += steps[rng.integers(8, size=len(active))]
        observed = positions[active]
        for j, site in enumerate(sites):
            hits[active, j] |= np.all(observed == site, axis=1)
        active = active[np.all((observed > 0) & (observed < 11), axis=1)]
    drain = (positions == 0) | (positions == 11)
    if not np.all(drain.sum(axis=1) == 1):
        raise ValueError("loss or ambiguity at the outer material stop")
    return [
        *drain.sum(axis=0).tolist(),
        *(positions[:, :2] >= 8).sum(axis=0).tolist(),
        *hits.sum(axis=0).tolist(),
    ], ticks


def run():
    panels = []
    for fast_axis, seed in ((None, 930823), (0, 930824)):
        steps, q = measured_steps(fast_axis)
        forecast = predictions(q)
        counts, ticks = acquire(steps, seed)
        panels.append(
            dict(
                extra_tags="idle" if fast_axis is None else "x forward and backward",
                seed=seed,
                squared_steps_per_appointment=(2 * q).tolist(),
                conductance=q.tolist(),
                forecast=forecast,
                counts=counts,
                last_arrival=ticks,
            )
        )
    return dict(
        shots_per_panel=SHOTS,
        events=[
            "x drain",
            "y drain",
            "z drain",
            "exit x>=8",
            "exit y>=8",
            "visit (8,5,5)",
            "visit (5,8,5)",
        ],
        family_size=14,
        simultaneous_99_frequency_radius=float(
            np.sqrt(np.log(2 * 14 / 0.01) / (2 * SHOTS))
        ),
        panels=panels,
    )


def check(actual, saved):
    if actual.keys() != saved.keys() or len(actual["panels"]) != len(saved["panels"]):
        raise ValueError("record fields or panel count differ")
    for key in actual.keys() - {"panels"}:
        if key == "simultaneous_99_frequency_radius":
            np.testing.assert_allclose(actual[key], saved[key], atol=1e-15, rtol=0)
        elif actual[key] != saved[key]:
            raise ValueError(f"acquisition metadata changed at {key}")
    for panel, old in zip(actual["panels"], saved["panels"], strict=True):
        if panel.keys() != old.keys():
            raise ValueError("panel fields differ")
        for key in panel:
            if key == "forecast":
                np.testing.assert_allclose(panel[key], old[key], atol=3e-14, rtol=0)
            elif panel[key] != old[key]:
                raise ValueError(f"simulation record changed at {key}")
        if (
            np.max(abs(np.array(old["counts"]) / SHOTS - panel["forecast"]))
            > actual["simultaneous_99_frequency_radius"]
        ):
            raise ValueError("frozen Green forecast misses the sensor records")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "results.json"
    if args.write:
        check(result, result)
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        check(result, json.loads(path.read_text()))
        print("Verified both complete material-guide banks and their held forecasts.")
