"""Measure the blank-background dispersion relative to an isolated sign clock."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import eigvalsh, expm

from . import measure as source

HERE = source.HERE
SHOTS = 2**26
SEED = 704189


def fringe_modes(lap):
    """Read the prepared amplitudes in contact-map port order."""
    preparation = (HERE / "fringe_modes.json").read_bytes()
    modes = np.array(json.loads(preparation))
    ev = eigvalsh(lap)
    ids = np.r_[0, np.flatnonzero(np.diff(ev) > 1e-8)[:3] + 1]
    spatial = ev[ids]
    assert modes.shape == (4, lap.shape[0])
    np.testing.assert_allclose(modes @ modes.T, np.eye(4), atol=1e-13, rtol=0)
    np.testing.assert_allclose(lap @ modes.T, modes.T * spatial, atol=1e-13, rtol=0)
    return spatial[1:], modes, sha256(preparation).hexdigest()


def run():
    rng = np.random.default_rng(SEED)
    board, seam, _, _ = source.matched.fp.make_board()
    lap = np.diag(board.degree) - board.adjacency().toarray()
    h, _, _ = source.matched.fp.generator(board, seam, 1, 8, False)
    native = h.toarray()
    np.testing.assert_array_equal(2 * (native - 48 * np.eye(board.n)), lap)
    spatial, modes, preparation_sha256 = fringe_modes(lap)
    times = np.array([1.0, 2.0])
    # Empty magnitude-one receiver, with both source parity generators checked.
    vacuum_even, _ = source.timing.model(0, 0)
    vacuum_odd, _ = source.timing.model(0, 2)
    np.testing.assert_array_equal(
        vacuum_odd - vacuum_even, 4 * np.eye(len(vacuum_even))
    )
    clock_prob = (1 + np.cos(times)) / 2
    clock_hits = rng.binomial(SHOTS, clock_prob)
    # 2 clock settings, 2 laws * 3 modes * 2 appointments * 2 phase settings.
    eps = np.sqrt(np.log(2 * 26 / 0.01) / (2 * SHOTS))
    clock = np.arccos(2 * clock_hits / SHOTS - 1)
    lo_t = np.arccos(np.clip(2 * (clock_hits / SHOTS + eps) - 1, -1, 1))
    hi_t = np.arccos(np.clip(2 * (clock_hits / SHOTS - eps) - 1, -1, 1))
    out = {}
    for label, op in (
        ("native", native),
        ("half_exchange", 48 * np.eye(board.n) + lap / 4),
    ):
        probs = np.empty((3, 2, 2))
        for i, b in enumerate(modes[1:]):
            a = modes[0]
            launch = (a + b) / np.sqrt(2)
            analysers = [(a + b) / np.sqrt(2), (a - 1j * b) / np.sqrt(2)]
            for j, t in enumerate(times):
                wave = expm(-1j * t * op) @ launch
                probs[i, j] = [abs(np.vdot(v, wave)) ** 2 for v in analysers]
        counts = rng.binomial(SHOTS, probs)
        x, y = 2 * counts[0, 0] / SHOTS - 1
        phase = np.arctan2(y, x)
        # A disk containing the entire simultaneous two-coordinate rectangle.
        phase_width = np.arcsin(2 * np.sqrt(2) * eps / np.hypot(x, y))
        A = phase / (spatial[0] * clock[0])
        lo = (phase - phase_width) / (spatial[0] * hi_t[0])
        hi = (phase + phase_width) / (spatial[0] * lo_t[0])
        forecast = []
        for i in range(3):
            for j in range(2):
                if i == 0 and j == 0:
                    continue
                arg = A * spatial[i] * clock[j]
                pred = (1 + np.array([np.cos(arg), np.sin(arg)])) / 2
                delta = max(
                    abs(arg - lo * spatial[i] * lo_t[j]),
                    abs(hi * spatial[i] * hi_t[j] - arg),
                )
                observed = counts[i, j] / SHOTS
                forecast.append(
                    dict(
                        mode=i + 1,
                        appointment=j + 1,
                        probabilities=pred.tolist(),
                        measured=observed.tolist(),
                        half_width=float(eps + delta / 2),
                        passed=bool(np.max(abs(pred - observed)) <= eps + delta / 2),
                    )
                )
        out[label] = dict(
            coefficient=A,
            coefficient_interval=[lo, hi],
            kinetic_mass=1 / (2 * A),
            kinetic_mass_interval=[1 / (2 * hi), 1 / (2 * lo)],
            yes_counts=counts.tolist(),
            forward_probabilities=probs.tolist(),
            held=forecast,
        )
    return dict(
        protocol_sha256=sha256((HERE / "MOBILE_STANDARD.md").read_bytes()).hexdigest(),
        fringe_modes_sha256=preparation_sha256,
        seed=SEED,
        shots=SHOTS,
        graph_eigenvalues=spatial.tolist(),
        clock_yes=clock_hits.tolist(),
        inferred_clock=clock.tolist(),
        clock_intervals=np.stack((lo_t, hi_t), axis=1).tolist(),
        binary_half_width=float(eps),
        laws=out,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = json.loads(json.dumps(run(), allow_nan=False))
    path = HERE / "mobile_standard.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2) + "\n")
    else:
        source.region.check(json.loads(path.read_text()), data)
    for name, result in data["laws"].items():
        print(name, result["coefficient_interval"], result["kinetic_mass_interval"])
