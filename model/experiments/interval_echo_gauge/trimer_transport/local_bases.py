"""Certify trimer singlet fidelity using four independent position-readout bases."""

import argparse
import json
from itertools import combinations
from math import pi

import numpy as np
from scipy.linalg import expm
from scipy.stats import beta

from . import measure

CASES = ("idle", "common_AB", "different_order", "dephased")


def bases():
    w = np.exp(2j * pi / 3)
    d = np.arange(3)
    fourier = w ** np.outer(d, d) / np.sqrt(3)
    return [np.eye(3)] + [fourier @ np.diag(w ** (k * d * d)) for k in range(3)]


def compile_basis(target):
    contacts, private = measure.operators(measure.KAPPA)

    def phases(angles):
        ans = np.eye(3, dtype=complex)
        for p, angle in enumerate(angles):
            wait = (-angle) % (2 * pi) / measure.KAPPA
            ans = expm(-1j * wait * private[p]) @ ans
        return ans

    def planar(a, b, gate):
        t = 2 * np.arctan2(abs(gate[0, 1]), abs(gate[0, 0]))
        p = np.angle(gate[0, 0])
        q = np.angle(gate[0, 1]) - pi / 2
        left = np.zeros(3)
        right = np.zeros(3)
        offset = np.zeros(3)
        left[[a, b]] = [(p + q) / 2, -(p + q) / 2]
        right[[a, b]] = [(p - q) / 2, -(p - q) / 2]
        offset[[a, b]] = t * (measure.KAPPA + 0.5)
        return (
            phases(left)
            @ phases(offset)
            @ expm(-1j * t * contacts[a, b])
            @ phases(right)
        )

    remaining = np.array(target, dtype=complex)
    inverses = []
    for col in range(2):
        for b in range(2, col, -1):
            a = b - 1
            x, y = remaining[[a, b], col]
            length = np.hypot(abs(x), abs(y))
            if abs(y) < 1e-14:
                continue
            gate = np.array([[x.conjugate(), y.conjugate()], [-y, x]]) / length
            remaining[[a, b]] = gate @ remaining[[a, b]]
            inverses.append((a, b, gate.conj().T))
    result = phases(np.angle(remaining.diagonal()))
    for a, b, gate in reversed(inverses):
        result = planar(a, b, gate) @ result
    np.testing.assert_allclose(result, target, rtol=0, atol=7e-14)
    return result


def binomial_interval(k, n, comparisons):
    tail = 0.01 / (2 * comparisons)
    lower = 0.0 if k == 0 else float(beta.ppf(tail, k, n - k + 1))
    upper = 1.0 if k == n else float(beta.ppf(1 - tail, k + 1, n - k))
    return [lower, upper]


def pair_intervals(bins, comparisons=48):
    n = sum(bins[2:])
    if sum(bins) != measure.STARTS or len(bins) != 29 or min(bins) < 0 or n == 0:
        raise ValueError("Missing source or terminal records")
    out = []
    for left, right in combinations(range(3), 2):
        k = sum(
            v
            for word, v in zip(measure.TUPLES, bins[2:], strict=True)
            if word[left] == word[right]
        )
        out.append(binomial_interval(k, n, comparisons))
    return out


def fidelity(intervals):
    lo, hi = np.asarray(intervals).sum(axis=(0, 1))
    # The sampler rounds each launch-bin probability to 14 decimal places.
    # Admission exceeds 3e-5. Rounding 27 bins changes each conditional event
    # by < 1e-8, so 1e-7 also covers the native basis matrix tolerance.
    return [max(0, 1 - float(hi) / 3 - 1e-7), min(1, 1 - float(lo) / 6 + 1e-7)]


def launch_yield(rows):
    bins = [row["counts"] for row in rows]
    total = len(bins) * measure.STARTS
    rejected = sum(row[0] for row in bins)
    excluded = sum(row[1] for row in bins)
    admitted = total - rejected - excluded
    # Use only the first frame for the acceptance interval. This
    # permits every frame to have a slightly different rounding error.
    p_lo, p_hi = binomial_interval(sum(bins[0][2:]), measure.STARTS, 52)
    f_lo, f_hi = fidelity([pair_intervals(row, 52) for row in bins])
    return dict(
        total_starts=total,
        filter_refused=rejected,
        orbital_refused=excluded,
        admitted=admitted,
        success_per_million_interval=[
            1e6 * max(0, p_lo - 5e-13) * f_lo,
            1e6 * min(1, p_hi + 5e-13) * f_hi,
        ],
    )


def run():
    sec, arms, _, _, _, _, finite, meta = measure.trimer.source()
    state = measure.tensor(sec, np.delete(arms, 3, axis=1), finite)
    source = meta["filter_acceptances"][-1]
    admit = float(np.vdot(state, state).real)
    direct = measure.native_word("AB")
    changed = measure.native_word("BA")
    frames = [compile_basis(b) for b in bases()]
    rng = np.random.default_rng(12113585)
    info = {}
    for case in CASES:
        rows = []
        for frame in frames:
            first = frame if case == "idle" else frame @ direct
            other = frame @ changed if case == "different_order" else first
            signal = measure.detector(state, first, other, other, case == "dephased")
            counts = measure.count(
                rng,
                measure.STARTS,
                np.r_[1 - source, source * (1 - admit), source * signal.ravel()],
            )
            rows.append(
                dict(counts=counts, equal_pair_intervals=pair_intervals(counts))
            )
        info[case] = dict(
            bases=rows,
            fidelity_interval=fidelity([row["equal_pair_intervals"] for row in rows]),
            launch_yield=launch_yield(rows),
        )
    return dict(starts=measure.STARTS, cases=info)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--write", action="store_true")
    data = run()
    dest = measure.HERE / "local.json"
    if p.parse_args().write:
        dest.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(data, json.loads(dest.read_bytes()), "local bases")
    for k, v in data["cases"].items():
        print(k, v["fidelity_interval"])
