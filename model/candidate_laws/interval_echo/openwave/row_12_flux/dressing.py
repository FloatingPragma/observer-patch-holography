"""Separate direct wall motion from paths through the complementary load words."""

import argparse
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare
from .measure import HERE, LENGTH, Ring, interpolate, ramsey


def energies(ring, m, flux, contact_cost):
    basis, source = ring.momenta[m]
    h = basis.conj().T @ ring.generator(flux) @ basis
    potential = abs(basis.T) ** 2 @ contact_cost
    keep = np.arange(len(potential)) != source
    gap = potential[source] - potential[keep]
    if np.any(gap < 15.99):
        raise ValueError("The zero-order wall has a degenerate off-wall channel")
    direct = h[source, source].real
    virtual = np.sum(abs(h[keep, source]) ** 2 / gap)
    return direct, direct + virtual


def dressed_bounds(ring, m, phi):
    basis, source = ring.momenta[m]
    h = basis.conj().T @ ring.generator(phi) @ basis
    keep = np.arange(len(h)) != source
    energy, vectors = np.linalg.eigh(h[np.ix_(keep, keep)])
    a = h[source, source].real
    weights = abs(vectors.conj().T @ h[keep, source]) ** 2
    if energy[-1] >= a:
        raise ArithmeticError("Wall level is not separated from its complement")

    def update(x):
        return a + np.sum(weights / (x - energy))

    lo = a
    for _ in range(2):
        hi = update(lo)
        lo = update(hi)
    return lo, hi


def fresh_acquisition(ring, rng):
    rows = []
    for m in (1, 2, 3, 4):
        for phi in (np.pi / 3, 3 * np.pi / 5, 3 * np.pi / 4):
            left = dressed_bounds(ring, m, phi)
            right = dressed_bounds(ring, m, -phi)
            predicted = [left[0] - right[1], left[1] - right[0]]
            ep, _, plus, _ = ring.line(phi, m)
            em, _, minus, _ = ring.line(-phi, m)
            if not predicted[0] - 1e-12 <= ep - em <= predicted[1] + 1e-12:
                raise ArithmeticError("Full eigenline outside its Schur enclosure")
            read = ramsey(np.sum(minus.conj() * plus, axis=1), rng)
            rows.append(
                dict(
                    m=m,
                    flux=phi,
                    predicted_bounds=predicted,
                    spectral=ep - em,
                    read=read,
                )
            )
    return rows


def run():
    initial = json.loads((HERE / "results.json").read_text())
    result = []
    rng = np.random.default_rng(120438)
    for n, charge in ((3, -1 / 3), (4, 2 / 3)):
        ring = Ring(n, candidate.KAPPA)
        control = Ring(n, 0)
        cost = (ring.h - control.h).diagonal()
        measured = next(
            s
            for s in initial["studies"]
            if s["load"] == n and s["contrast"] == candidate.KAPPA
        )
        zero = [energies(ring, m, 0, cost)[0] for m in range(LENGTH)]
        rows = []
        for point in measured["held"]:
            m, phi = point["m"], point["flux"]
            plus = np.array(energies(ring, m, phi, cost))
            minus = np.array(energies(ring, m, -phi, cost))
            direct, dressed = plus - minus
            k = 6 * np.pi * m / LENGTH
            shift = 3 * charge * phi / LENGTH
            rigid_direct = interpolate(zero, 3, k - shift) - interpolate(
                zero, 3, k + shift
            )
            rows.append(
                dict(
                    m=m,
                    flux=phi,
                    direct=direct,
                    rigid_direct=rigid_direct,
                    second_order=dressed,
                    spectral=point["dominant_splitting"],
                    observed=point["read"]["frequency"],
                )
            )
        result.append(
            dict(load=n, charge=charge, rows=rows, fresh=fresh_acquisition(ring, rng))
        )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write("\n")
        return
    stored = json.loads((HERE / "dressing.json").read_text())
    compare(json.loads(json.dumps(result, allow_nan=False)), stored, "dressed_flux")
    print("Reproduced both wall dispersions and all 24 held dressed-wall splittings.")


if __name__ == "__main__":
    main()
