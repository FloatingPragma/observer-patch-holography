"""Follow a wall around one applied flux cycle and compare equal-load controls."""

import argparse
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare
from .measure import HERE, LENGTH, Ring, SHOTS


def certify(ring):
    # An absolute row-sum bound holds for the Hermitian flux derivative at any twist.
    speed = float(max(np.sum(abs(ring.h * ring.d), axis=1)) / LENGTH)
    radius = np.pi / 64
    min_upper, max_lower, min_wall_gap, max_coupling = np.inf, -np.inf, np.inf, 0.0
    for phi in 2 * np.pi * np.arange(64) / 64:
        h = ring.generator(phi)
        for basis, source in ring.momenta:
            block = basis.conj().T @ h @ basis
            vals = np.linalg.eigvalsh(block)
            keep = np.arange(len(block)) != source
            complement = np.linalg.eigvalsh(block[np.ix_(keep, keep)])
            min_upper = min(min_upper, vals[-1])
            max_lower = max(max_lower, vals[-2])
            min_wall_gap = min(
                min_wall_gap, block[source, source].real - complement[-1]
            )
            max_coupling = max(max_coupling, np.linalg.norm(block[keep, source]))
    # For the upper vector, ||Qv||/|Pv| <= ||QHP||/(a - max spectrum(QHQ)).
    gap = min_wall_gap - 2 * radius * speed
    weight = 1 / (1 + ((max_coupling + radius * speed) / gap) ** 2) if gap > 0 else 0.0
    return dict(
        speed_bound=speed,
        uniform_band_gap=min_upper - max_lower - 2 * radius * speed,
        uniform_wall_weight_lower=weight,
    )


def ramp(ring, m, load, rng):
    basis, source = ring.momenta[m]
    # Word zero is the lexicographically first contiguous configuration.
    compact = np.flatnonzero(abs(basis[0]) > 0).item()
    x = np.eye(basis.shape[1], dtype=complex)[:, [source, compact]]
    for j in range(128):
        h = basis.conj().T @ ring.generator(2 * np.pi * (j + 0.5) / 128) @ basis
        e, u = np.linalg.eigh(h)
        x = u @ (np.exp(-1j * e * 16 / 128)[:, None] * (u.conj().T @ x))
    final = ring.unitary.conj()[:, None] * (basis @ x)
    out = {}
    for column, name in enumerate(("wall", "compact")):
        probs = []
        for other, wall in ring.momenta:
            amp = other.conj().T @ final[:, column]
            seen = abs(amp) ** 2
            probs.extend((seen[wall], sum(np.delete(seen, wall))))
        probs = np.array(probs).reshape(LENGTH, 2)
        expected = (m - load) % LENGTH
        if np.sum(probs[expected]) < 1 - 1e-10:
            raise ValueError("Final momentum violates the large-gauge identity")
        counts = rng.multinomial(SHOTS, (probs / probs.sum()).ravel().round(12))
        out[name] = dict(
            momentum_wall_other_counts=counts.reshape(LENGTH, 2).tolist(),
            wall_probability=float(probs[:, 0].sum()),
        )
    return dict(
        initial_momentum=m,
        predicted_final_momentum=expected,
        wall_translation_shift=((3 * (expected - m) + 5) % 10) - 5,
        inputs=out,
    )


def run():
    rng = np.random.default_rng(120439)
    out = []
    for load in (3, 4):
        for k in (candidate.KAPPA, 0):
            ring = Ring(load, k)
            out.append(
                dict(
                    load=load,
                    contrast=k,
                    certificate=certify(ring),
                    ramps=[ramp(ring, m, load, rng) for m in (1, 3)],
                )
            )
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = run()
    if args.output:
        with args.output.open("x") as file:
            json.dump(data, file, indent=2, allow_nan=False)
            file.write("\n")
    else:
        compare(
            json.loads(json.dumps(data)),
            json.loads((HERE / "cycle.json").read_text()),
            "cycle",
        )
    print(
        "Verified the flux-cycle band estimates and both final momentum and wall flags."
    )


if __name__ == "__main__":
    main()
