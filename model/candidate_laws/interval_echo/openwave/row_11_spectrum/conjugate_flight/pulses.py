"""Rebuild or verify nearest-neighbour preparation words on a closed material route."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares

from .measure import HERE, LAW, seam_data


def contacts(n):
    if n == 3:
        return [(0, 1), (1, 2)] * 6
    return [(p, (p + 1) % n) for p in range(n)] * 4


def native_word(state, word, kappa, sign, *, inverse=False):
    h, _ = seam_data(kappa)
    labels = [LAW.SRC.index((1 + sign, 1)), LAW.SRC.index((1, 1 + sign))]
    block = h[np.ix_(labels, labels)]
    pulses = list(zip(contacts(len(state)), word, strict=True))
    v = state.copy()
    if inverse:
        pulses = [(pair, (-t) % (2 * np.pi)) for pair, t in reversed(pulses)]
    for (a, b), time in pulses:
        v[[a, b]] = expm(-1j * time * block) @ v[[a, b]]
    return v


def search(target, seed):
    """Discover exposure times with a bounded, reproducible least-squares search."""
    n = len(target)
    bonds = contacts(n)

    def propagate(times, derivative):
        v = np.zeros(n, complex)
        v[1] = 1
        jac = np.zeros((n, len(bonds)), complex)
        for j, ((a, b), t) in enumerate(zip(bonds, times, strict=True)):
            z = np.exp(-8.5j * t)
            c, s = z * np.cos(t / 2), 1j * z * np.sin(t / 2)
            gate = np.array([[c, s], [s, c]])
            v[[a, b]] = gate @ v[[a, b]]
            if derivative:
                jac[[a, b]] = gate @ jac[[a, b]]
                jac[[a, b], j] = -1j * (8.5 * v[[a, b]] - 0.5 * v[[b, a]])
        if derivative:
            return np.vstack((jac.real, jac.imag))
        residual = v - target
        return np.r_[residual.real, residual.imag]

    rng = np.random.default_rng(seed)
    for _ in range(4):
        result = least_squares(
            lambda t: propagate(t, False),
            rng.uniform(0, 2 * np.pi, len(bonds)),
            jac=lambda t: propagate(t, True),
            gtol=2e-13,
            ftol=2e-13,
            xtol=2e-13,
            max_nfev=900,
        )
        if np.linalg.norm(result.fun) < 1e-10:
            return (result.x % (2 * np.pi)).round(14).tolist()
    raise ArithmeticError("Target was not synthesized within the search budget")


def compile_words():
    n = 15
    words = {
        str(m): search(
            np.exp(2j * np.pi * m * np.arange(n) / n) / np.sqrt(n), 47133 + m
        )
        for m in (0, 1, 2)
    }
    words["bond_minus"] = search(np.array([1, -1, 0], complex) / np.sqrt(2), 8122)
    return words


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    data = json.dumps(compile_words(), indent=2) + "\n"
    if p.parse_args().write:
        Path(HERE / "pulses.json").write_text(data)
    else:
        print(data, end="")
