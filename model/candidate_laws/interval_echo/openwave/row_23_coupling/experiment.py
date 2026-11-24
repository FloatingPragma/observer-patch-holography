"""Return interferometry of a native empty-write scatterer on complete federations."""

import argparse
import json
from itertools import combinations
from math import factorial
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import shortest_path
from scipy.stats import beta

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_30_formed import (
    experiment as formed,
)
from model.executable_model import carrier

HERE = Path(__file__).resolve().parent
N = 2**36
RADIAL_N = 2**23
IDENTITY = np.eye(2)
PAIR = list(combinations(range(3), 2))
PARITY = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
STRENGTHS = {"blank": None, "written": candidate.KAPPA, "cut_only": 0}
NODES, WEIGHTS = np.polynomial.legendre.leggauss(6)
TIMES = (np.arange(160)[:, None] * 0.5 + 0.25 * (NODES + 1)).ravel()
WEIGHTS = np.tile(WEIGHTS * 0.25, 160)


def encoded(a):
    """Keep the complex observable as separately named real and imaginary parts."""
    arr = np.asarray(a)
    return {"real": arr.real.round(8).tolist(), "imag": arr.imag.round(8).tolist()}


def geometry(board, seam):
    lap = board.laplacian().toarray() / 2
    energy, vec = eigh(lap)
    rc = min(set(board.neighbours[seam[0]]) & set(board.neighbours[seam[1]]))
    ports = [*seam, rc]
    targets = [0, 0.5, 1, 1.5, 2, 2.5]
    indices = sorted({int(np.argmin(abs(energy - k * k / 2))) for k in targets})
    shells = np.maximum(energy[indices], 0)
    probes = list(np.eye(board.n)[:, ports].T)
    probes.extend((probes[i] + probes[j]) / np.sqrt(2) for i, j in PAIR)
    cuts = []
    for k in (0.75, 2.25):
        keep = energy <= k * k / 2 + 1e-12
        source = vec[:, keep] @ vec[rc, keep]
        weight = np.vdot(source, source).real
        probes.append(source / np.sqrt(weight))
        cuts.append((k, keep, weight))
    return lap, energy, vec, ports, shells, np.asarray(probes), cuts


def returns(energy, vectors, sources, rng):
    amps = np.exp(-1j * np.outer(TIMES, energy)) @ ((sources @ vectors) ** 2).T
    # An ancillary two-rail comparison with a spatially unresolved outcome.
    p = np.stack([(1 + amps.real) / 2, (1 + amps.imag) / 2], axis=0)
    return rng.binomial(N, np.clip(p.round(11), 0, 1)).astype(np.int64)


def radial(board, energy, vec, ports, rng):
    delta = energy[:, None] / 32 - np.fft.fftfreq(256)[None, :]
    # Normalized finite Fourier clock amplitude at each bin.
    kernel = np.exp(-1j * np.pi * 255 * delta) * np.sinc(256 * delta) / np.sinc(delta)
    accepted = (32 * np.fft.fftfreq(256) < -2) & (32 * np.fft.fftfreq(256) > -12)
    dist = (
        shortest_path(
            board.adjacency(), directed=False, indices=ports[:2], unweighted=True
        )
        .min(axis=0)
        .astype(int)
    )
    banks = []
    for i in ports[:2]:
        wave = vec @ (vec[i, :, None] * kernel)
        weights = (abs(wave) ** 2)[:, accepted].sum(axis=1)
        bins = np.bincount(dist, weights=weights)
        p = np.r_[bins, max(0, 1 - bins.sum())]
        banks.append(rng.multinomial(RADIAL_N, p / p.sum()).tolist())
    counts = np.sum(banks, axis=0)
    radius = np.sqrt(
        np.dot(counts[:-1], (np.arange(len(counts) - 1) + 0.5) ** 2) / counts[:-1].sum()
    )
    tag = float(counts[:-1].sum() / (2 * RADIAL_N))
    eps = np.sqrt(np.log(2 * 20000 / 0.01) / (4 * RADIAL_N))
    # Layer-cake decomposition of the nonnegative squared radius; use
    # one-sided exact binomial bounds for each pooled tail on each source.
    # Bounding the joint moment directly covers the differing source laws.
    radii_sq = (np.arange(len(counts) - 1) + 0.5) ** 2
    increments = np.diff(radii_sq, prepend=0)
    upper_numerators = []
    for bank in banks:
        tail = np.cumsum(np.array(bank[:-1])[::-1])[::-1]
        up = beta.ppf(1 - 0.01 / 20000, tail + 1, RADIAL_N - tail)
        upper_numerators.append(increments @ np.where(tail == RADIAL_N, 1, up))
    upper_radius = np.sqrt(np.mean(upper_numerators) / (tag - eps))
    return {
        "radius_upper_simultaneous_99": float(upper_radius),
        "counts_by_distance_then_failed_tag": banks,
        "radius_links": float(radius),
        "accepted_fraction": float(counts[:-1].sum() / (2 * RADIAL_N)),
    }


def measure(raw, z):
    real, imag = 2 * raw / N - 1
    weights = -1j * np.exp(1j * z * TIMES) * WEIGHTS
    values = weights @ (real + 1j * imag)
    matrix = np.diag(values[:3])
    for p, (i, j) in enumerate(PAIR):
        matrix[i, j] = matrix[j, i] = values[p + 3] - (values[i] + values[j]) / 2
    return matrix, values[6:]


def allowance(z):
    # Each component is a sum of independent bounded quadrature means.
    # 20000 covers all decoded scalar quadratures in both experiments.
    err = np.sqrt(
        2
        * np.log(2 * 20000 / 0.01)
        * np.sum((WEIGHTS * np.exp(-z.imag * TIMES)) ** 2)
        / N
    )
    # Gaussian remainder on each half-unit interval, then the missing tail.
    bound = max(abs(z + 9), abs(z - 6))
    coeff = factorial(6) ** 4 / (13 * factorial(12) ** 3)
    quad = coeff * 0.5**13 * bound**12 / (1 - np.exp(-z.imag / 2))
    tail = np.exp(-80 * z.imag) / z.imag
    return float(np.sqrt(2) * (err + quad + tail))


def norm(x):
    return float(np.linalg.norm(x, 2))


def inverse_error(matrix, error):
    q = norm(np.linalg.inv(matrix))
    if q * error >= 1:
        raise ValueError("Unresolved inverse in this experiment")
    return q * q * error / (1 - q * error)


def interaction(blank, target, z):
    b, w = blank[:2, :2], target[:2, :2]
    ib = np.linalg.inv(b)
    vertex = ib @ (w - b) @ ib
    inferred = ib - np.linalg.inv(w)
    err = 3 * allowance(z)  # two row entries bounded by epsilon and 2 epsilon
    dinv = inverse_error(b, err)
    vertex_err = (norm(ib) + dinv) ** 2 * 2 * err + dinv * (2 * norm(ib) + dinv) * norm(
        w - b
    )
    return vertex, inferred, inverse_error(b, err) + inverse_error(w, err), vertex_err


def dressed(v, g):
    return np.linalg.solve(IDENTITY - v @ g, v)


def held_prediction(blank, potential, potential_err, eps):
    b, x = blank[:2, :2], blank[2:, :2]
    vertex = dressed(potential, b)
    resolvent = np.linalg.inv(IDENTITY - potential @ b)
    matrix_err = potential_err * (norm(b) + 3 * eps) + norm(potential) * 3 * eps
    inverse_bound = norm(resolvent) / (1 - norm(resolvent) * matrix_err)
    if inverse_bound < 0:
        raise ValueError("Unresolved held-response interval")
    vertex_err = inverse_bound * (potential_err + matrix_err * norm(vertex))
    xerr = 2 * np.sqrt(2) * eps
    pred_err = (
        eps
        + (norm(x) + xerr) ** 2 * vertex_err
        + xerr * (2 * norm(x) + xerr) * norm(vertex)
    )
    return (
        blank[2, 2] + (x @ vertex @ x.T).item(),
        blank[2, 2] + (x @ potential @ x.T).item(),
        pred_err,
    )


def shell_read(raw, control, shells, radius):
    panels = []
    for width in (0.5, 0.25):
        rows = []
        for energy in shells:
            z = energy + 1j * width
            b = measure(raw, z)[0]
            w = measure(control, z)[0]
            t, v, err, terr = interaction(b, w, z)
            if not rows:
                v0, err0 = v, err
            pred, born, e_pred = held_prediction(b, v0, err0, allowance(z))
            delta = float(abs(w[2, 2] - pred))
            q = float(np.sqrt(2 * energy))
            rows.append(
                {
                    "k": q,
                    "k_radius": q * radius,
                    "vertex_parity": encoded(PARITY.T @ t @ PARITY),
                    "vertex_operator_error_bound": terr,
                    "inferred_parity": encoded(PARITY.T @ v @ PARITY),
                    "inferred_operator_error_bound": err,
                    "held_return": encoded(w[2, 2]),
                    "single_pass_error": float(abs(born - w[2, 2])),
                    "held_prediction_error": delta,
                    "held_joint_error_bound": e_pred + allowance(z),
                    "constant_interaction_difference": norm(v - v0),
                    "constant_interaction_joint_bound": err + err0,
                }
            )
        panels.append({"width": width, "shells": rows})
    return panels


def band_match(raw, control, energy, vec, ports, cuts):
    z = -1 + 0.5j
    _, held_blank = measure(raw, z)
    _, held_write = measure(control, z)
    z0 = 0.5j
    _, v, v_err, _ = interaction(measure(raw, z0)[0], measure(control, z0)[0], z0)
    local = vec[ports[:2], :]
    studies = []
    for j, (k, keep, weight) in enumerate(cuts):
        g_high = (local[:, ~keep] / (z - energy[~keep])) @ local[:, ~keep].T
        g_low = (local[:, keep] / (z - energy[keep])) @ local[:, keep].T
        x = (
            (vec[ports[2], keep] / (z - energy[keep]))
            @ local[:, keep].T
            / np.sqrt(weight)
        )
        matched = dressed(v, g_high)
        if not j:
            coarse = matched
            coarse_high = g_high

        def predicted(coeff):
            return held_blank[j] + x @ dressed(coeff, g_low) @ x

        # Propagate the calibration uncertainty before decimation. The exact
        # block identity combines high and low returns into the full g.
        def budget(g):
            t = dressed(v, g)
            inv = norm(np.linalg.inv(IDENTITY - v @ g))
            margin = 1 - inv * v_err * norm(g)
            if margin <= 0:
                raise ValueError("Unresolved effective coefficient")
            err = inv * v_err * (1 + norm(g @ t)) / margin
            return 2 * allowance(z) + float(np.linalg.norm(x) ** 2) * err

        studies.append(
            {
                "cutoff_k": k,
                "retained_modes": int(keep.sum()),
                "source_projection_weight": float(weight),
                "matched_parity": encoded(PARITY.T @ matched @ PARITY),
                "held_filtered_return": encoded(held_write[j]),
                "matched_error": float(abs(held_write[j] - predicted(matched))),
                "matched_joint_bound": budget(g_low + g_high),
                "coarse_joint_bound": budget(g_low + coarse_high),
                "bare_joint_bound": budget(g_low),
                "bare_in_reduced_space_error": float(abs(held_write[j] - predicted(v))),
                "coarse_in_reduced_space_error": float(
                    abs(held_write[j] - predicted(coarse))
                ),
            }
        )
    return studies


def compute():
    formed.native()
    rng = np.random.default_rng(230471)
    results, archive = [], {}
    for lvl in (0, 1):
        board = carrier.federation(lvl)
        seam = carrier.EMITTER_SEAM[lvl]
        lap, e0, v0, ports, shells, sources, cuts = geometry(board, seam)
        raw = {}
        for name, k in STRENGTHS.items():
            if k is None:
                energy, vec = e0, v0
            else:
                h = lap.copy()
                h[np.ix_(seam, seam)] -= k * IDENTITY + np.outer(
                    PARITY[:, 1], PARITY[:, 1]
                )
                energy, vec = eigh(h)
            raw[name] = returns(energy, vec, sources, rng)
            archive[f"level{lvl}_{name}"] = raw[name]
            if name == "written":
                rad = radial(board, energy, vec, ports, rng)
        z_probes = [measure(raw["blank"], p + 0.5j)[0] for p in shells]
        # Audit the bare inverse on inputs independent from the derivation.
        exact = [(v0[ports, :] / (p + 0.5j - e0)) @ v0[ports, :].T for p in shells]
        assert all(
            np.max(abs(a - b)) < 2 * allowance(p + 0.5j)
            for a, b, p in zip(z_probes, exact, shells, strict=True)
        )
        results.append(
            {
                "cells": board.n // 12,
                "ports": ports,
                "core": rad,
                "frequency_shells": shells.tolist(),
                "blank_band_gaps": [
                    [float(a), float(b)]
                    for a, b in zip(e0[:-1], e0[1:], strict=True)
                    if b - a > 0.15
                ],
                "written": shell_read(
                    raw["blank"], raw["written"], shells, rad["radius_links"]
                ),
                "cut_only": shell_read(
                    raw["blank"], raw["cut_only"], shells, rad["radius_links"]
                ),
                "two_resolutions": band_match(
                    raw["blank"], raw["written"], e0, v0, ports, cuts
                ),
            }
        )
    return {
        "ramsey_enrollments_per_quadrature_and_hold": N,
        "studies": results,
    }, archive


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, help="write JSON and raw NPZ with the same stem"
    )
    args = parser.parse_args()
    result, counts = compute()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        np.savez_compressed(args.output.with_suffix(".npz"), **counts)
    else:
        stored = json.loads((HERE / "results.json").read_text())
        with np.load(HERE / "results.npz") as raw:
            for key in counts:
                np.testing.assert_allclose(counts[key], raw[key], atol=2, rtol=0)

        def compare(a, b):
            if isinstance(a, dict):
                assert a.keys() == b.keys()
                for key in a:
                    compare(a[key], b[key])
            elif isinstance(a, list):
                for x, y in zip(a, b, strict=True):
                    compare(x, y)
            else:
                np.testing.assert_allclose(a, b, atol=2e-7, rtol=1e-8)

        compare(result, stored)
        print(
            "Reproduced both federations, source radii, shell spectra and held filtered probes."
        )


if __name__ == "__main__":
    main()
