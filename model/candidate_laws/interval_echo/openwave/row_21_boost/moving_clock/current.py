"""Independent route current, spatial fringe and stationary line clock."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from . import apparatus as lab
from . import readback as data_only
from .measure import HERE, SHOTS, TARGETS, histories


def route_meter(psi, rail, exchange):
    # Restore two adjacent cells to a local pair: swap the intervening rail
    # with the second selected input, then mix. Dwell is calibrated to the
    # exchange strength of the reader, independently of the flight.
    n = len(psi)
    ports = psi.reshape(-1).copy()
    a = rail
    middle = (a + 1) % (2 * n)
    b = (a + 2) % (2 * n)
    h = np.array(
        [[32 + exchange / 2, -exchange / 2], [-exchange / 2, 32 + exchange / 2]]
    )
    full = expm(-1j * np.pi * h / exchange)
    quarter = expm(-0.5j * np.pi * h / exchange)
    ports[[middle, b]] = full @ ports[[middle, b]]
    # Native unmatched-site phase, calibrated with a uniform two-rail source.
    # A full exchange's unequal-load phase is exp(-32 i pi/exchange).
    phase = np.exp(-32j * np.pi / exchange)
    corrected = np.array([ports[a], ports[middle] / phase])
    rows = []
    for tag in (1, 1j):
        out = abs(quarter @ (corrected * np.array([1, tag]))) ** 2
        rows.append([out[0], out[1], 1 - out.sum()])
    return rows


def acquisition():
    rng = np.random.default_rng(211943)
    result = []
    for level in (1, 2):
        n = len(lab.mount(level)) // 2
        for scale, tag in ((1.0, "native"), (0.8, "changed")):
            layers = lab.gates(scale, 8.0)
            for target in TARGETS:
                h = round(target * n / (2 * np.pi))
                z = lab.arms(n, h, layers, scale, 8.0, plane=True)
                branches = []
                for psi in z:
                    first = psi @ layers[0].T
                    last = lab.step(psi, layers)
                    rails = [
                        lab.draw((abs(s) ** 2).sum(axis=0), SHOTS, rng).tolist()
                        for s in (psi, first, last)
                    ]
                    bins = [
                        lab.draw(p, SHOTS, rng).tolist()
                        for p in route_meter(psi, int(h < 0), scale)
                    ]
                    branches.append(dict(rails=rails, fringes=bins))
                hist = histories(z, layers)
                phase = []
                for a, b in hist:
                    phase.append(
                        [
                            int(
                                rng.binomial(
                                    SHOTS,
                                    abs(a[0, int(h < 0)] + p * b[0, int(h < 0)]) ** 2
                                    / 2,
                                )
                            )
                            for p in (1, 1j, -1, -1j)
                        ]
                    )
                result.append(
                    dict(
                        key=f"{level}_{tag}_{h}",
                        level=level,
                        n=n,
                        h=h,
                        exchange=scale,
                        branches=branches,
                        origin_phase=phase,
                    )
                )
    return dict(shots=SHOTS, seed=211943, rows=result)


def spatial(counts, shots):
    c = np.array(counts, dtype=float)
    if np.any(c.sum(axis=1) != shots):
        raise ValueError("Missing spatial comparison exits")
    # Difference after the iX quarter rotation is -2 Im(conj(a)b).
    # Second setting rotates the second amplitude by i. Its contrast reads the real part.
    im, re = (c[:, 1] - c[:, 0]) / shots
    err = data_only.uncertainty(c[:, :2], shots).sum(axis=1)
    z = re + 1j * im
    bound = np.arcsin(np.linalg.norm(err) / abs(z))
    return float(np.angle(z)), float(bound)


def read(data):
    nshot = data["shots"]
    refs = {}
    out = []
    for row in data["rows"]:
        readings = []
        for arm in row["branches"]:
            c = np.array(arm["rails"], dtype=float)
            if np.any(c.sum(axis=1) != nshot):
                raise ValueError("Missing rail read")
            v = (c[2, 0] - c[1, 0]) / nshot
            dv = data_only.uncertainty(c[[2, 1], 0], nshot).sum()
            k, dk = spatial(arm["fringes"], nshot)
            readings.append([v, dv, k, dk])
        v, dv = (np.array(readings[0])[:2] + np.array(readings[1])[:2]) / 2
        k = (readings[0][2] - readings[1][2]) / 2
        dk = (readings[0][3] + readings[1][3]) / 2
        sig = np.array(row["origin_phase"], dtype=float)
        z = (sig[:, 0] - sig[:, 2] + 1j * (sig[:, 3] - sig[:, 1])) / nshot
        err = data_only.uncertainty(sig, nshot)
        du = np.arcsin(np.hypot(err[:, 0] + err[:, 2], err[:, 1] + err[:, 3]) / abs(z))
        increments = np.angle(z[1:] / z[:-1])
        if np.any(abs(increments) + du[:-1] + du[1:] >= np.pi):
            raise ValueError("Stationary phase branch unresolved")
        e = increments.sum() / (2 * lab.T)
        de = (du[0] + du[-1]) / (2 * lab.T)
        key = (row["n"], row["exchange"])
        if row["h"] == 0:
            refs[key] = e, de
        m, dm = refs[key]
        ideal = (e - k * v) / m
        box = [
            (a - b * c) / f
            for a, b, c, f in product(
                [e - de, e + de], [k - dk, k + dk], [v - dv, v + dv], [m - dm, m + dm]
            )
        ]
        edges = [np.sqrt(1 - q * q) for q in (v - dv, v + dv)]
        if -dv <= v <= dv:
            edges.append(1.0)
        out.append(
            dict(
                key=row["key"],
                branches=readings,
                velocity=float(v),
                velocity_error=float(dv),
                winding=k,
                winding_error=dk,
                energy=e,
                energy_error=de,
                ideal=float(ideal),
                ideal_interval=[min(box), max(box)],
                lorentz_interval=[min(edges), max(edges)],
            )
        )
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path)
    dest = p.parse_args().output
    if dest:
        dest.write_text(json.dumps(acquisition(), indent=1, allow_nan=False) + "\n")
    else:
        print(
            json.dumps(read(json.loads((HERE / "current.json").read_text())), indent=1)
        )
