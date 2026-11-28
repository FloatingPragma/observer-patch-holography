"""Resolve energy on the original 99-state written seam without translating its records."""

import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path

import numpy as np

from model.experiments.interval_echo_einstein_station.mobile_collar_source import (
    measure as formation,
)

HERE = Path(__file__).parent
LAW = formation.emission.LAW


def four_product(twice_rho):
    t = twice_rho.reshape(11, 9, 11, 9)
    left = t.trace(axis1=1, axis2=3)
    right = t.trace(axis1=0, axis2=2)
    state = np.kron(left, right)
    assert np.trace(state) == 4
    return state


def energy_projectors(k):
    h = formation.emission.four_h(k)
    eye = np.eye(LAW.DIM, dtype=np.int64)
    levels = sorted({0, k, k + 1, 4 * k, 4 * k + 1})
    operators = []
    for e in levels:
        p = eye.copy()
        den = 1
        for other in levels:
            if other != e:
                p = p @ (h - 4 * other * eye)
                den *= 4 * (e - other)
        assert np.array_equal(h @ p, 4 * e * p)
        operators.append((e, p, den))
    np.testing.assert_allclose(sum(p / n for _, p, n in operators), eye, atol=1e-14)
    return h, operators


def spectrum(four, operators):
    probs = {e: F(int((four * p.T).sum()), 4 * den) for e, p, den in operators}
    assert sum(probs.values()) == 1 and min(probs.values()) >= 0
    return probs


def exchange_memories(h, branches):
    """Physically swap the two sharp memories, retaining both entire load pairs."""
    before = after = restored = 0
    for r, s in product(branches, repeat=2):
        m, loads = divmod(int(r), 9)
        n, other = divmod(int(s), 9)
        u, v = 9 * n + loads, 9 * m + other
        before += int(h[r, r] + h[s, s])
        after += int(h[u, u] + h[v, v])
        back = (9 * (v // 9) + u % 9, 9 * (u // 9) + v % 9)
        assert back == (r, s)
        restored += int(h[back[0], back[0]] + h[back[1], back[1]])
    return {
        "before": F(before, 16),
        "after": F(after, 16),
        "restored": F(restored, 16),
        "work_on_two_contacts": F(after - before, 16),
    }


def run_case(k):
    h, projectors = energy_projectors(k)
    echo = LAW.echo()
    pattern = LAW._write_pattern()
    results = []
    for i, (a, b) in enumerate(LAW.SRC):
        v = np.r_[np.zeros(LAW.N_SRC, dtype=int), pattern[:, i]]
        np.testing.assert_allclose(echo[:, i], v / np.sqrt(2), atol=3e-15)
        coherent = np.outer(v, v)
        sharp = np.diag(v * v)
        blank = np.zeros_like(coherent)
        blank[i, i] = 2
        row = dict(loads=[a, b], ordinary_incident=F(int(h[i, i]), 4))
        for name, state in (
            ("sharp", sharp),
            ("coherent", coherent),
            ("unwritten", blank),
        ):
            ref = four_product(state)
            delta = ref - 2 * state
            t = delta.reshape(11, 9, 11, 9)
            assert (
                not t.trace(axis1=1, axis2=3).any()
                and not t.trace(axis1=0, axis2=2).any()
            )
            formed, erased = spectrum(2 * state, projectors), spectrum(ref, projectors)
            gain = sum(e * (erased[e] - formed[e]) for e in formed)
            assert gain == F(int((h * delta.T).sum()), 16)
            if name == "unwritten":
                assert gain == 0
            row[name] = dict(
                retained_law=formed,
                independent_law=erased,
                retained_energy=sum(e * p for e, p in formed.items()),
                independent_energy=sum(e * p for e, p in erased.items()),
                gain=gain,
            )
        row["memory_swap"] = exchange_memories(h, np.flatnonzero(v))
        assert row["memory_swap"]["work_on_two_contacts"] == 2 * row["sharp"]["gain"]
        results.append(row)
    return results


def calibrate():
    # Full native spectral register from a finite Fourier clock. A controlled
    # dwell and readable clock are supplied. Each level is an integer below 64.
    h, operators = energy_projectors(8)
    evec = np.linalg.eigh(h / 4)
    e, v = evec
    bins = np.arange(64)
    amplitudes = np.fft.fft(np.exp(2j * np.pi * np.outer(bins, e) / 64), axis=0) / 64
    source = {}
    for tag, op, den in operators:
        inferred = (v * abs(amplitudes[tag]) ** 2) @ v.T
        np.testing.assert_allclose(inferred, op / den, atol=3e-13)
        source[tag] = round(float(np.trace(inferred)))
    assert sum(source.values()) == 99
    # One fair proposal per eigenword and one 7/8 shutter per energy unit.
    s = F(7, 8)
    normal = sum(n * s**val for val, n in source.items())
    ref = sum(float(s**val / normal) * (p / den) for val, p, den in operators)
    levels, axes = np.linalg.eigh(ref)
    K = (axes * (-np.log(levels))) @ axes.T
    beta = np.log(8 / 7)
    error = np.max(abs(K - np.log(float(normal)) * np.eye(99) - beta * h / 4))
    assert error < 2e-12
    return dict(
        level_multiplicities=source,
        survival_fraction=float(normal / 99),
        modular_slope=beta,
        maximum_operator_residual=error,
    )


def encode(obj):
    if isinstance(obj, F):
        return str(obj)
    raise TypeError(type(obj).__name__)


def run():
    return dict(
        stiffness_eight=run_case(8), exchange_only=run_case(0), thermometer=calibrate()
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    outcome = json.loads(json.dumps(run(), default=encode))
    path = HERE / "source.json"
    if args.write:
        path.write_text(json.dumps(outcome, indent=2) + "\n")
    else:
        saved = json.loads(path.read_text())
        assert (
            outcome["stiffness_eight"] == saved["stiffness_eight"]
            and outcome["exchange_only"] == saved["exchange_only"]
        )
        assert (
            outcome["thermometer"]["level_multiplicities"]
            == saved["thermometer"]["level_multiplicities"]
        )
        for key in ("survival_fraction", "modular_slope"):
            np.testing.assert_allclose(
                outcome["thermometer"][key], saved["thermometer"][key]
            )
