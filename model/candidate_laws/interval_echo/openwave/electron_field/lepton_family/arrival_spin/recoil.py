"""Read offset-independent dispersion masses from a single cyclic address and a port counter."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm, hankel

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (  # noqa: E402
    anchor_spectrum as spectrum,
    rest_partner as rest,
)

SHOTS = 2**24
SURVEY = np.arange(640)
BLOCK = np.arange(256)
DELAYS = 2 ** np.arange(10, 30)
PHASES = np.array([1, 1j, -1, -1j])


def apparatus(board, route, j, effect):
    b = rest.guide.path_incidence(board, route)
    n = len(route)
    zero = np.ones(n) / np.sqrt(n)
    spatial = np.sqrt(2 / n) * np.cos(j * np.pi * (np.arange(n) + 0.5) / n)
    np.testing.assert_allclose(b @ zero, 0, atol=1e-14)
    kinetic = b @ spatial
    rate = np.linalg.norm(kinetic)
    np.testing.assert_allclose(b.T @ kinetic, rate**2 * spatial, atol=1e-13)
    # Vertex and edge columns span the excited wave under the full generator.
    basis = np.zeros((n + len(b), 3))
    basis[:n, 0] = zero
    basis[:n, 1] = spatial
    basis[n:, 2] = kinetic / rate
    masses, chars = np.linalg.eigh(effect)
    for m in masses:
        h = rest.guide.generator(b, np.array([[m]]))
        reduced = np.array([[m, 0, 0], [0, m, rate], [0, rate, -m]])
        np.testing.assert_allclose(h @ basis, basis @ reduced, atol=2e-13)
        np.testing.assert_allclose(
            expm(-0.375j * h) @ basis, basis @ expm(-0.375j * reduced), atol=2e-13
        )
    # The uniform source is on vertices and the first edge is anchored at port zero.
    weights = abs(chars[0]) ** 2
    a, b0 = zero[0], spatial[0]

    def acquire(times, rng):
        probs = np.zeros((len(times), 4))
        for m, w in zip(masses, weights, strict=True):
            # Canonical clock rates make the long-wait seeded replay portable.
            energy = np.round(np.hypot(rate, m), 14)
            m = np.round(m, 14)
            c, s = np.cos(energy * times), np.sin(energy * times) / energy
            first = a * np.exp(-1j * m * times)
            excited = b0 * (c - 1j * m * s)
            edge = kinetic[0] * s
            probs += (
                w
                * (
                    abs(first[:, None] + excited[:, None] * PHASES) ** 2
                    + edge[:, None] ** 2
                )
                / 2
            )
        return rng.binomial(SHOTS, np.clip(probs.round(11), 0, 1))

    return acquire, 2 * a * b0


def fringe(counts, normalization):
    return ((counts[:, 0] - counts[:, 2]) + 1j * (counts[:, 1] - counts[:, 3])) / (
        SHOTS * normalization
    )


def read(survey, blocks, norm):
    z = fringe(survey, norm)
    matrix = hankel(z[:192], z[191:])
    u, s, _ = np.linalg.svd(matrix, full_matrices=False)
    # Union bound on four binomial port probabilities at every survey wait, then Frobenius.
    eps = np.sqrt(np.log(8 * len(z) / 0.001) / (2 * SHOTS))
    cutoff = np.sqrt(8 * matrix.size) * eps / norm
    rank = int(np.count_nonzero(s > cutoff))
    if not rank:
        raise ValueError("No resolved recoil line")
    basis = u[:, :rank]
    roots = np.linalg.eigvals(np.linalg.lstsq(basis[:-1], basis[1:], rcond=None)[0])
    omega = np.sort(np.angle(roots))
    for delay, counts in zip(DELAYS, blocks, strict=True):
        arg = np.angle(spectrum.amplitudes(BLOCK, fringe(counts, norm), omega))
        omega = (arg + 2 * np.pi * np.rint((delay * omega - arg) / (2 * np.pi))) / delay
    return omega, {
        "cutoff": float(cutoff),
        "first_ten_singular_values": s[:10].tolist(),
    }


def fit_first_two(n, first, second):
    shape = 4 * np.sin(np.pi * np.arange(1, 4) / (2 * n)) ** 2
    ratio = shape[1] / shape[0]
    mass = (second**2 - ratio * first**2) / (2 * (ratio * first - second))
    stiffness = (first**2 + 2 * mass * first) / shape[0]
    if min(mass) <= 0 or min(stiffness) <= 0:
        raise ValueError("The two recoil records give no positive dispersion mass")
    prediction = np.sqrt(mass**2 + stiffness * shape[2]) - mass
    return mass, stiffness, prediction


def acquire_study(board, route, effect, seed):
    rng = np.random.default_rng(seed)
    rows = []
    positive = []
    for j in (1, 2, 3):
        acquire, norm = apparatus(board, route, j, effect)
        survey = acquire(SURVEY, rng)
        blocks = [acquire(delay + BLOCK, rng) for delay in DELAYS]
        freq, diag = read(survey, blocks, norm)
        # For a common kinetic coefficient the heaviest mass takes the least energy for the same kick.
        positive.append(freq[freq > 0][::-1])
        rows.append(
            {
                "j": j,
                "normalization": float(norm),
                "all_frequencies": freq.tolist(),
                "rank": diag,
                "survey_counts": spectrum.pack(survey),
                "refinement_counts": spectrum.pack(blocks),
            }
        )
        if j == 2:
            mass, stiffness, third = fit_first_two(len(route), *positive)
    return {
        "masses_from_positive_gaps": mass.tolist(),
        "stiffness_from_positive_gaps": stiffness.tolist(),
        "ratios": (mass / mass[0]).tolist(),
        "predicted_third_positive_recoil": third.tolist(),
        "measured_third_positive_recoil": positive[2].tolist(),
        "heldout_residual": (positive[2] - third).tolist(),
        "reads": rows,
    }


def end_cost_check(board, route, effect, source):
    # A local energy on the two terminal incidence grades preserves the guide reversal,
    # all transitions and the uniform positive rest states. Test the held spatial mode.
    b = rest.guide.path_incidence(board, route)
    n = len(route)
    uniform = np.r_[np.ones(n) / np.sqrt(n), np.zeros(len(b))]
    spatial = np.r_[
        np.sqrt(2 / n) * np.cos(3 * np.pi * (np.arange(n) + 0.5) / n), np.zeros(len(b))
    ]
    inputs = (uniform[:, None] + spatial[:, None] * PHASES) / np.sqrt(2)
    mass = np.array(source["masses_from_positive_gaps"])
    pos = np.array(source["predicted_third_positive_recoil"])
    predicted = np.r_[pos, -pos - 2 * mass]
    design = np.exp(1j * np.outer(SURVEY, predicted))
    norm = 2 * uniform[0] * spatial[0]
    bars = {}
    rng = np.random.default_rng(178035)
    eigen, chars = np.linalg.eigh(effect)
    for cost in (0.0, 1.0):
        probs = np.zeros((len(SURVEY), 4))
        for m, w in zip(eigen, abs(chars[0]) ** 2, strict=True):
            h = rest.guide.generator(b, np.array([[m]]))
            for index in (n, len(h) - 1):
                h[index, index] += cost
            np.testing.assert_allclose(h @ uniform, m * uniform, atol=2e-13)
            e, v = np.linalg.eigh(h)
            projection = v.T @ inputs
            phase = np.exp(-1j * np.outer(SURVEY, e))
            for detector in (0, n):
                probs += w * abs((phase * v[detector]) @ projection) ** 2
        counts = rng.binomial(SHOTS, np.clip(probs.round(11), 0, 1))
        z = fringe(counts, norm)
        half = len(z) // 2
        coeff = np.linalg.lstsq(design[:half], z[:half], rcond=None)[0]
        residue = z[half:] - design[half:] @ coeff
        bars[str(cost)] = {
            "terminal_cost": cost,
            "counts": spectrum.pack(counts),
            "held_half_rms_complex_fringe_error": float(
                np.sqrt(np.mean(abs(residue) ** 2))
            ),
        }
    return bars


def run():
    board, route = rest.guide.longest_unpinned_route()
    effect, _ = rest.guide.register_effect()
    family = acquire_study(board, route, effect, 178033)
    return {
        "shots_per_phase": SHOTS,
        "route": route,
        "native_exposure_delays": DELAYS.tolist(),
        "seed": 178033,
        "family": family,
        "terminal_cost_control": end_cost_check(board, route, effect, family),
        "scalar_control": acquire_study(
            board, route, np.eye(3) * np.trace(effect).real / 3, 178034
        ),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("x") as f:
        json.dump(run(), f, indent=2, allow_nan=False)
        f.write("\n")
