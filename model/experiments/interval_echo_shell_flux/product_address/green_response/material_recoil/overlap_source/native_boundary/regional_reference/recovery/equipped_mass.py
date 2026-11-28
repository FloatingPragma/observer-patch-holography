"""Transport the blank kinetic calibration to actual equipped loading-standard contacts."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply

from . import mobile_standard as blank

source = blank.source
HERE = blank.HERE


def run():
    rng = np.random.default_rng(437111)
    old = blank.run()
    board, _, _, _ = source.matched.fp.make_board()
    lap = np.diag(board.degree) - board.adjacency().toarray()
    ev, basis = eigh(lap)
    ids = np.flatnonzero(np.diff(ev) > 1e-8)[:3] + 1
    probs = np.zeros((3, 2, 2))
    discrepancy = 0.0
    for parity in range(8):
        four, _ = source.timing.model(1, parity)
        h = four / 4
        en, modes = eigh(h)
        for i, idx in enumerate(ids):
            a, b = basis[:, 0], basis[:, idx]
            launch = np.kron(((a + b) / np.sqrt(2))[:, None], np.eye(8))
            readers = [(a + b) / np.sqrt(2), (a - 1j * b) / np.sqrt(2)]
            for j, t in enumerate((1.0, 2.0)):
                wave = modes @ (np.exp(-1j * t * en)[:, None] * (modes.T @ launch))
                audit = expm_multiply(-1j * t * sparse.csr_matrix(h), launch)
                discrepancy = max(discrepancy, float(abs(wave - audit).max()))
                for k, an in enumerate(readers):
                    projected = np.einsum(
                        "p,pij->ij", an.conj(), wave.reshape(board.n, 8, 8)
                    )
                    probs[i, j, k] += (abs(projected) ** 2).sum() / 64
    assert discrepancy < 1e-10
    hits = rng.binomial(blank.SHOTS, probs)
    clock_hits = rng.binomial(blank.SHOTS, (1 + np.cos([1.0, 2.0])) / 2)
    # Four blank calibration readouts suffice (two fringes and original
    # clock at both appointments); two new clocks and twelve held outcomes.
    eps = np.sqrt(np.log(2 * 18 / 0.01) / (2 * blank.SHOTS))

    def clock_bounds(n):
        m = 2 * n / blank.SHOTS - 1
        return np.arccos(np.clip(m + 2 * eps, -1, 1)), np.arccos(
            np.clip(m - 2 * eps, -1, 1)
        )

    x, y = 2 * np.array(old["laws"]["native"]["yes_counts"][0][0]) / blank.SHOTS - 1
    phi = np.arctan2(y, x)
    ephi = np.arcsin(2 * np.sqrt(2) * eps / np.hypot(x, y))
    tlo, thi = clock_bounds(np.array(old["clock_yes"]))
    amin = (phi - ephi) / (ev[ids[0]] * thi[0])
    amax = (phi + ephi) / (ev[ids[0]] * tlo[0])
    tlo, thi = clock_bounds(clock_hits)
    rows = []
    for i, idx in enumerate(ids):
        for j in range(2):
            low, high = ev[idx] * amin * tlo[j], ev[idx] * amax * thi[j]
            mid = (low + high) / 2
            forecast = (1 + np.array([np.cos(mid), np.sin(mid)])) / 2
            measured = hits[i, j] / blank.SHOTS
            # Both Cartesian fringe probabilities have slope at most 1/2.
            width = eps + (high - low) / 4
            fringe = 2 * measured - 1
            floor = float(np.linalg.norm(np.maximum(abs(fringe) - 2 * eps, 0)))
            cap = float(np.linalg.norm(abs(fringe) + 2 * eps))
            rows.append(
                dict(
                    shell=i + 1,
                    appointment=j + 1,
                    forecast=forecast.tolist(),
                    measured=measured.tolist(),
                    comparison_half_width=float(width),
                    excluded_components=(abs(forecast - measured) > width).tolist(),
                    signed_fringe_length=float(np.linalg.norm(fringe)),
                    length_interval=[floor, cap],
                    pure_phase_excluded=cap < 1 or floor > 1,
                )
            )
    return dict(
        protocol_sha256=sha256((HERE / "EQUIPPED_MASS.md").read_bytes()).hexdigest(),
        seed=437111,
        trials=blank.SHOTS,
        clock_yes=clock_hits.tolist(),
        yes_counts=hits.tolist(),
        probabilities=probs.tolist(),
        source_free_kinetic_interval=[amin, amax],
        joint_probability_half_width=float(eps),
        direct_amplitude_residual=discrepancy,
        rows=rows,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = json.loads(json.dumps(run(), allow_nan=False))
    path = HERE / "equipped_mass.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2) + "\n")
    else:
        source.region.check(json.loads(path.read_text()), data)
    print(
        [
            (r["excluded_components"], r["length_interval"])
            for r in data["rows"]
        ]
    )
