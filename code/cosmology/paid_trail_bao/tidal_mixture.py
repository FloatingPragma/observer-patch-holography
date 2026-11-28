"""Find verifiable joint-radius obstructions for independent annular scale mixtures."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import math

import numpy as np
from scipy.optimize import linprog

from . import tidal_scale_kernel as kernel
from . import tidal_scales as scales

HERE = scales.HERE


def support_witness(angles_lo, angles_hi, gains, lower, upper):
    # gains are exact positive surface-density and ordinary-conversion factors.
    # Angles are signed integers in units of 1/ANGLE. The remaining divisor is pi.
    n, annuli, bins = angles_lo.shape
    norm = np.maximum(
        1, np.maximum(abs(np.array(lower, float)), abs(np.array(upper, float)))
    )
    gain = np.asarray(gains, float) / (kernel.ANGLE * np.pi)
    low = angles_lo * gain[None, :, None] / norm[:, None, None]
    high = angles_hi * gain[None, :, None] / norm[:, None, None]

    def flatten(cube):
        return cube.reshape(n, annuli * bins).T

    # s and t weight lower and upper requirements at the selected observing
    # radii. z bounds the signed sum at every length in a given source annulus.
    A = np.column_stack(
        (flatten(high), -flatten(low), -np.repeat(np.eye(annuli), bins, axis=0))
    )
    cost = np.r_[
        -np.array(lower, float) / norm, np.array(upper, float) / norm, np.ones(annuli)
    ]
    eq = np.r_[np.ones(2 * n), np.zeros(annuli)][None, :]
    sol = linprog(
        cost,
        A_ub=A,
        b_ub=np.zeros(len(A)),
        A_eq=eq,
        b_eq=[1.0],
        bounds=[(0, None)] * (2 * n) + [(None, None)] * annuli,
        method="highs",
    )
    assert sol.success, sol.message
    if sol.fun >= -1e-8:
        return dict(certified_exclusion=False)
    scaled = np.maximum(0, sol.x[: 2 * n]) / np.tile(norm, 2)
    s, t = np.rint(scaled * (10**7 / scaled.max())).astype(np.int64).reshape(2, n)
    return verify(s, t, angles_lo, angles_hi, gains, lower, upper)


def verify(s, t, angles_lo, angles_hi, gains, lower, upper):
    assert np.all(s >= 0) and np.all(t >= 0) and all(g >= 0 for g in gains)
    assert (int(s.sum()) + int(t.sum())) * 2 * kernel.ANGLE < 2**63
    # Replay the proposed separating weights with integer angles, rational
    # decimal data and an outward bound on pi. The solver gives no verdict.
    ticks = np.tensordot(s, angles_hi, axes=(0, 0)) - np.tensordot(
        t, angles_lo, axes=(0, 0)
    )
    peak = sum(
        g * int(p) / kernel.ANGLE for g, p in zip(gains, ticks.max(1), strict=True)
    )
    pi_low = Q("3.14159265358979323846264338327950288")
    pi_high = pi_low + Q(1, 10**35)
    capacity = peak / (pi_low if peak >= 0 else pi_high)
    need = sum(
        int(a) * x - int(b) * y for a, b, x, y in zip(s, t, lower, upper, strict=True)
    )
    assert capacity < need
    normalization = int(max(s.max(), t.max()))
    return dict(
        certified_exclusion=True,
        lower_weights={str(i): int(x) for i, x in enumerate(s) if x},
        upper_weights={str(i): int(x) for i, x in enumerate(t) if x},
        required_sum=str(need),
        outward_capacity=str(capacity),
        normalized_gap=float((need - capacity) / normalization),
    )


def run(stored=None):
    path = HERE / "tidal_scales.json"
    envelope = json.loads(path.read_bytes())
    inputs, *_ = scales.budget.inventory.rp.galaxy.transfer.budget.load_inputs()
    profiles, mu = scales.source_data()
    G = Q(
        str(scales.budget.inventory.rp.galaxy.transfer.budget.local.profiles.galaxies.G)
    )
    knots, K, refs = profiles["source"]
    moments = kernel.kernel_bins(knots, K)
    maps = scales.pixel_maps()
    out = {}
    for name, panel in envelope["disks"].items():
        pending = [
            (int(cal), tag)
            for cal, table in panel["source"].items()
            for tag, stats in table.items()
            if stats["separate_annulus_ring_upper"] == stats["rings"]
        ]
        if not pending:
            continue
        r, lo, hi = scales.target(inputs[name])
        ticks = math.lcm(*(x.denominator for x in r))
        edges = np.array([int(x * ticks) for x in r])
        lengths = kernel.length_edges(float(r[-1]))
        ones = np.ones((1, len(r)), dtype=np.int64)
        boxes = [
            kernel.disk_boxes(
                int(x), edges, ones, ones, lengths, ticks, moments, keep_annuli=True
            )
            for x in edges
        ]
        al, ah = map(np.asarray, zip(*boxes, strict=True))
        area = [b * b - a * a for a, b in zip([0] + r[:-1], r, strict=True)]
        out[name] = {}
        for cal, tag in pending:
            stock = [Q(str(x)) for x in maps[tag][name]]
            factor = G * mu * moments[-1] / refs[cal]
            gains = [
                (b - a) / A * factor
                for a, b, A in zip([0] + stock[:-1], stock, area, strict=True)
            ]
            key = f"{cal}/{tag}"
            saved = stored["joint"][name][key] if stored is not None else {}
            if saved.get("certified_exclusion"):

                def weights(field):
                    return np.array(
                        [saved[field].get(str(i), 0) for i in range(len(r))], np.int64
                    )

                out[name][key] = verify(
                    weights("lower_weights"),
                    weights("upper_weights"),
                    al,
                    ah,
                    gains,
                    lo,
                    hi,
                )
            else:
                out[name][key] = support_witness(al, ah, gains, lo, hi)
    return dict(sha256={path.name: sha256(path.read_bytes()).hexdigest()}, joint=out)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    dst = HERE / "tidal_mixture.json"
    old = None if args.write else json.loads(dst.read_bytes())
    answer = run(old)
    if args.write:
        dst.write_text(json.dumps(answer, indent=2) + "\n")
    else:
        assert answer == old
    print(
        {
            n: sum(v["certified_exclusion"] for v in row.values())
            for n, row in answer["joint"].items()
        }
    )
