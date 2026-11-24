"""Neutral-subtracted pair energies from a complete static endpoint menu."""

import argparse
import json

import numpy as np
from scipy.linalg import svdvals

from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import federation
from . import experiment as ex


def minimum(n, left, right, j, h):
    diagonal = np.full(n + 1, j)
    diagonal[0] *= abs(left)
    diagonal[-1] *= abs(right)
    b = np.diag(diagonal) - h * np.eye(n + 1, k=-1)
    values = svdvals(b)
    parity = (-1) ** n * left * right
    return float(-sum(values) + (2 * values[-1] if parity < 0 else 0))


def run():
    rng = np.random.default_rng(201127)
    j, h = ex.couplings(candidate.KAPPA)
    _, path, _ = ex.guide(federation(2))
    rows = []
    for n in range(2, (len(path) - 2) // 2 + 1):
        signs = (-1, 0, 1)
        e = np.array([[minimum(n, a, b, j, h) for b in signs] for a in signs])
        u = e + e[1, 1] - e[:, 1, None] - e[None, 1, :]
        reads = [[ex.phase_reader(u[a, b], rng) for b in range(3)] for a in range(3)]
        measured = np.array([[x["energy"] for x in row] for row in reads])
        corner = measured[::2, ::2]
        rows.append(
            dict(
                n=n,
                exact=u.tolist(),
                counts=reads,
                even=corner.mean(),
                charge_product=np.mean(corner * np.outer([-1, 1], [-1, 1])),
            )
        )
    return dict(shots=ex.SHOTS, signs=[-1, 0, 1], rows=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = ex.portable(run())
    path = ex.HERE / "neutral.json"
    if args.write:
        path.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    elif json.loads(path.read_text()) != result:
        raise ValueError("Neutral-subtracted energies differ")
    print("Verified nine signed boundary pairs on all 44 tested lengths.")


if __name__ == "__main__":
    main()
