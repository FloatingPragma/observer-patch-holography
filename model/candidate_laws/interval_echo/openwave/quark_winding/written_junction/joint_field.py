"""Inclusive signed field read on the same stationary eigenstate as the triple-wall record."""

import argparse
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import currents
from . import backgrounds, field_probe


def cores(n):
    sites = np.flatnonzero(n)
    for a, b in zip(sites, np.roll(sites, -1)):
        gap = (b - a) % len(n)
        if gap not in (1, 4):
            continue
        region = (a, b) if gap == 1 else ((a + 2) % len(n),)
        f = np.zeros(len(n))
        for p in region:
            for step in (-1, 0, 1):
                f[(p + step) % len(n)] += 1 / 3
        q = (2 if gap == 1 else -1) / 3
        assert abs(f @ (n - 1 / 3) - q) < 1e-12
        yield ("positive" if gap == 1 else "negative"), f, q


def instrument(path, reference, heads, f, sign, theta=np.pi / 2):
    # Convert forward ring bits into the absolute address convention of the
    # microscopic probe. Repeat the first address only to read the closing link.
    cycle = list(path) + [path[0]]
    forward = np.diff(cycle) > 0
    before = np.where(forward, reference, 1 - reference)
    after = np.where(forward, heads, 1 - heads)
    return field_probe.interferometer(before, after, cycle, np.r_[f, f[0]], theta, sign)


def unit_control(path):
    # A separate ordinary unit moves from site 0 onto a flat electrode at site 1.
    f = np.zeros(len(path))
    f[1] = 1
    old = np.zeros(len(path), int)
    old[0] = 1
    return {str(s): instrument(path, old, np.zeros_like(old), f, s) for s in (-1, 1)}


def panel(row, site, ext, geometry, rng):
    length = row["length"]
    path = geometry["ports"]
    g = backgrounds.bits(row["gauss"], length)
    data = backgrounds.block(g, np.array(geometry["reference"]))
    w, heads, loads, a, tight = data[:5]
    refs = {"zero_head": np.zeros(length, int)}
    for phase in range(3):
        bulk = (np.arange(length) % 3 == phase).astype(int)
        if row["n"] > length // 3:
            bulk[site] += 1
        ids = np.flatnonzero((loads == bulk).all(axis=1))
        if len(ids):
            refs["phase_" + str(phase)] = heads[ids[0]]
    result = dict(
        length=length,
        n=row["n"],
        background=row["gauss"],
        references=list(refs),
        tight_probability=row["tight"],
        reads={},
    )
    if not tight.any():
        return result
    feat, state = backgrounds.terminal(data, ext)
    prob = state**2
    for name, ref in refs.items():
        for core in ("positive", "negative"):
            for sign in (1, -1):
                y = np.zeros((2, 4))
                total = 0
                field = 0
                equals = 0
                for i in np.flatnonzero(tight):
                    readings = [(f, q) for tag, f, q in cores(loads[i]) if tag == core]
                    for f, q in readings:
                        wt = prob[i] / len(readings)
                        delta = f @ (heads[i] - ref - np.roll(heads[i] - ref, 1))
                        assert (
                            abs(delta - f @ (loads[i] - (g + ref - np.roll(ref, 1))))
                            < 1e-12
                        )
                        record = instrument(path, ref, heads[i], f, sign)
                        expected = np.angle(np.exp(1j * np.pi / 2 * sign * delta))
                        assert (
                            abs(np.angle(np.exp(1j * (record["phase"] - expected))))
                            < 1e-11
                        )
                        y += wt * np.array(
                            [record["quadratures"][k] for k in ("x", "y")]
                        )
                        total += wt
                        field += wt * delta
                        equals += wt * (abs(delta - q) < 1e-10)
                np.testing.assert_allclose(y.sum(axis=1), total, atol=1e-10)
                inclusive = np.column_stack((y, np.full(2, max(0.0, 1 - total))))
                shots = 2**19
                counts = [
                    rng.multinomial(shots, p / p.sum()).tolist() for p in inclusive
                ]
                fringe = (y[:, 0] - y[:, 1]) / total
                result["reads"][f"{name}/{core}/{sign}"] = dict(
                    mean_field_at_core=field / total,
                    fraction_with_target=equals / total,
                    phase_units=np.angle(fringe[0] + 1j * fringe[1])
                    / (sign * np.pi / 2),
                    visibility=np.linalg.norm(fringe),
                    joint_xy=inclusive,
                    counts=counts,
                )
    return result


def run(scans):
    rng = np.random.default_rng(12013002)
    out = {}
    for name, pair in scans["paired"].items():
        out[name] = [
            panel(
                r,
                pair["additional_site"],
                name == "extension",
                scans["geometry"][str(r["length"])],
                rng,
            )
            for r in pair["selected"]
        ]
    return currents._portable(
        dict(
            unit=unit_control(scans["geometry"]["12"]["ports"]),
            shots=2**19,
            outcome_order=["bright", "dark", "loss_a", "loss_b", "no_tight_triple"],
            laws=out,
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    s = json.loads(Path(__file__).with_name("backgrounds.json").read_text())
    d = run(s)
    with a.output.open("x") as f:
        f.write(json.dumps(d, indent=2, allow_nan=False) + "\n")
