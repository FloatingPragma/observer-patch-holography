"""Source tomography by uniform exterior illumination on the protected two-cell rig."""

import argparse
from fractions import Fraction
from math import factorial, log, sqrt
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix

from .. import exterior

HERE = Path(__file__).resolve().parent
TIMES = (1, 2, 4, 8, 16)
SHOTS = 2**24
EPS = sqrt(log(2 * 2400 / 0.01) / (2 * SHOTS))


def trace_sites(columns):
    """Trace over the mediator and equally weighted, independent launch columns."""
    waves = columns.reshape(10, 10, -1)
    return np.einsum("ari,asi->rs", waves, waves.conj()) / waves.shape[-1]


def generator(k, occupation, tag, blocks, cost, parked=False):
    a = exterior.belt(blocks[k, 1])
    a[5, 5] += cost[occupation, 1]
    b = exterior.belt(blocks[k, tag])
    h = np.kron(np.zeros_like(a) if parked else a, np.eye(10)) + np.kron(np.eye(10), b)
    h[2, 2] += cost[1, tag]
    return a, h


def distributions(a, h):
    """One conditional exterior channel per mediator launch, at each held wait."""
    med = exterior.unitary(a, 12)
    inputs = np.kron(med, np.eye(10))
    e, v = eigh(h)
    spectral = v.T @ inputs
    out = {}
    for time in TIMES:
        waves = (v @ (np.exp(-1j * time * e)[:, None] * spectral)).reshape(100, 10, 10)
        probabilities = np.stack(
            [trace_sites(waves[:, j]).diagonal().real for j in range(10)]
        )
        np.testing.assert_allclose(probabilities.sum(axis=0), 1.0, rtol=0, atol=1e-13)
        np.testing.assert_allclose(probabilities.sum(axis=1), 1.0, rtol=0, atol=1e-13)
        out[time] = probabilities
    return out


def acquisition():
    blocks, _ = exterior.base.calibration()
    original = np.random.default_rng(202611259)
    held = np.random.default_rng(202611260)
    rows = []
    for k in (8, 0):
        cost = exterior.local_occupation_cost(k)
        for pole in (0, 1):
            for tag in (1, 2):
                for parked in (False, True):
                    a, h = generator(k, pole, tag, blocks, cost, parked)
                    for time, probs in distributions(a, h).items():
                        # Separate fresh preparations choose the launch uniformly inside its group.
                        curve = probs[[0, 5]]
                        others = probs[[1, 2, 3, 4, 6, 7, 8, 9]].mean(axis=0)
                        counts = [
                            original.multinomial(SHOTS, p).tolist() for p in curve
                        ]
                        fresh = held.multinomial(SHOTS, others)
                        forecast = (1 - np.sum(counts, axis=0) / SHOTS) / 8
                        err = np.max(abs(forecast - fresh / SHOTS))
                        rows.append(
                            dict(
                                k=k,
                                pole=pole,
                                tag=tag,
                                parked=parked,
                                time=time,
                                primary_site_counts=counts,
                                held_site_counts=fresh.tolist(),
                                held_max_error=float(err),
                            )
                        )
    return rows


def evolve_integer(twice, re, im, time, steps=128, order=16, bits=44):
    """Round every complex Taylor term on a fixed dyadic mesh; bound errors in 2-norm."""
    # For a symmetric generator, ||H||_2 <= the absolute row-sum norm.
    norm = int(np.abs(twice).sum(axis=1).max())
    if norm > 2 * steps:
        raise ValueError("increase the step count to keep the scaled norm at most one")
    scale = 2**bits
    H = csr_matrix(twice)
    # Round each real and imaginary component down by at most one mesh unit.
    # In the Taylor recurrence a term error feeds all later terms with total
    # multiplier at most e. The tail on an input of norm at most two is < 6/(m+1)!.
    component_norm = int(np.ceil(np.sqrt(2 * len(twice))))
    local = Fraction(3 * order * component_norm, scale) + Fraction(
        6, factorial(order + 1)
    )
    accumulated = steps * time * local
    if accumulated >= Fraction(1, 4):
        raise ValueError("error bound would not preserve the induction norm")
    for _ in range(steps * time):
        dr, di = re, im
        re, im = re.copy(), im.copy()
        for j in range(1, order + 1):
            # A term has norm at most 2/j! up to accumulated rounding. The
            # 3*scale component check also guards the packed integer multiply.
            if norm * 3 * scale >= 2**63:
                raise OverflowError("packed integer matrix action would overflow")
            dr, di = (H @ di) // (2 * steps * j), (-H @ dr) // (2 * steps * j)
            re += dr
            im += di
        if max(abs(re).max(), abs(im).max()) > 2 * scale:
            raise OverflowError("bounded propagation exceeded its component range")
    return re, im, accumulated


def certified_gate():
    """Enclose the marked receiver at t=4 on both occupied/empty sharp boundary launches."""
    blocks, _ = exterior.base.calibration()
    cost = exterior.local_occupation_cost(8)
    scale = 2**44
    result = []
    for pole in (0, 1):
        a, h = generator(8, pole, 1, blocks, cost)
        # In a blank single-load belt every hopping entry is exactly -1/2.
        # In the protected even memory the coincident units shift energy by -32.
        cycle = np.eye(10, dtype=np.int64)
        ring = -np.roll(cycle, 1, axis=1) - np.roll(cycle, -1, axis=1)
        med_twice = ring.copy()
        med_twice[5, 5] -= 64 * pole
        joint_twice = np.kron(med_twice, cycle) + np.kron(cycle, ring)
        joint_twice[2, 2] -= 64
        np.testing.assert_allclose(2 * a, med_twice, atol=1e-13, rtol=0)
        np.testing.assert_allclose(2 * h, joint_twice, atol=1e-13, rtol=0)
        u = cycle[:, [0]] * scale
        ar, ai, err1 = evolve_integer(med_twice, u, np.zeros_like(u), 12)
        ident = np.eye(10, dtype=np.int64)
        r, i, err2 = evolve_integer(
            joint_twice, np.kron(ar, ident), np.kron(ai, ident), 4
        )
        prob = Fraction(
            sum(int(x) ** 2 for x in np.r_[r[2::10].ravel(), i[2::10].ravel()]),
            10 * scale**2,
        )
        # A projection has unit norm. Each corresponding output vector differs by
        # at most err1+err2 and both vector norms sum to less than three.
        error = 3 * (err1 + err2)
        result.append(
            [
                float(np.nextafter(float(prob - error), -np.inf)),
                float(np.nextafter(float(prob + error), np.inf)),
            ]
        )
    return result


def run():
    rows = acquisition()
    return dict(
        chart=exterior.chart(),
        shots=SHOTS,
        times=list(TIMES),
        simultaneous_probability_radius=EPS,
        held_forecast_radius=1.25 * EPS,
        certified_gate_probability_empty_then_occupied=certified_gate(),
        rows=rows,
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    record = run()
    msg = json.dumps(record, indent=2, allow_nan=False) + "\n"
    if ap.parse_args().write:
        (HERE / "results.json").write_text(msg)
    else:
        print(msg, end="")
