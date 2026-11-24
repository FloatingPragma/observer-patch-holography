"""A mobile load binding two stationary native writes at surveyed separations."""

import argparse
import json
from pathlib import Path

import numpy as np

from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import eigsh

from model.candidate_laws.interval_echo.openwave.row_30_formed import (
    experiment as formed,
)
from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare
from model.executable_model import EMITTER_SEAM, federation
from .measure import HERE, K, SHOTS


def survey(board, first):
    degrees = sorted(board.degree[list(first)])
    d = shortest_path(
        board.adjacency(), directed=False, unweighted=True, indices=list(first)
    ).min(axis=0)
    out = {}
    for edge in sorted(board.seams):
        if sorted(board.degree[list(edge)]) != degrees:
            continue
        r = int(min(d[list(edge)])) + 1
        if 2 <= r <= 12:
            out.setdefault(r, edge)
    return out


def ground(lap, seams, strength):
    # Empty writes remove the native blank kinetic and contrast terms at both endpoints.
    h = lap.tolil(copy=True)
    for a, b in seams:
        h[a, a] -= strength + 0.5
        h[b, b] -= strength + 0.5
        h[a, b] += 0.5
        h[b, a] += 0.5
    h = h.tocsr()
    source = np.ones(h.shape[0]) / np.sqrt(h.shape[0])
    e, v = eigsh(h, k=1, which="SA", v0=source, tol=2e-14)
    return float(e[0]), float(np.linalg.norm(h @ v[:, 0] - e[0] * v[:, 0]))


def measured_calibration(rows, board):
    rng = np.random.default_rng(311943 + board.cells)
    gain = rows[3]["gain"]
    bank = []
    phase_error = np.arcsin(np.sqrt(2) * np.sqrt(2 * np.log(16 / 0.01) / SHOTS))
    estimates = []
    for wait in (0.125, 8):
        phase = np.exp(1j * gain * wait)
        quadratures = np.array([phase.real, phase.imag])
        counts = rng.binomial(SHOTS, (1 + quadratures) / 2)
        real, imag = 2 * counts / SHOTS - 1
        angle = np.arctan2(imag, real)
        if estimates:
            angle += 2 * np.pi * round((wait * estimates[0] - angle) / (2 * np.pi))
        estimates.append(float(angle / wait))
        bank.append(dict(wait=wait, real_imag_plus=counts.tolist()))
    if (phase_error / 0.125) * 8 + phase_error >= np.pi:
        raise ArithmeticError("Calibration cannot identify its winding")
    estimate = estimates[-1]
    error = float(phase_error / 8)
    minimum_dark = float(np.cos(np.pi * error / (2 * estimate)) ** 2)
    for r, sample in rows.items():
        if r <= 3:
            continue
        wait = np.pi * r / (3 * estimate)
        p = np.sin(wait * sample["gain"] / 2) ** 2
        sample["fresh_meter_calibration"] = dict(
            wait=wait, dark=int(rng.binomial(SHOTS, round(p, 12)))
        )
    return dict(
        seed=311943 + board.cells,
        quadratures=bank,
        inferred_gain=estimate,
        gain_sampling_allowance=error,
        coulomb_dark_probability_floor=minimum_dark,
    )


def scan(level, rng):
    board = federation(level)
    first = EMITTER_SEAM[level]
    lap = board.laplacian().asfptype() / 2
    available = survey(board, first)
    e0, err0 = ground(lap, [first], K)
    rows = {}
    for r, second in available.items():
        e1, err1 = ground(lap, [second], K)
        e12, err12 = ground(lap, [first, second], K)
        # With contrast removed the one-load operator is a graph Laplacian;
        # its uniform mode is the lowest mode for all three preparations.
        cut = lap.copy().tolil()
        for a, b in (first, second):
            cut[a, a] -= 0.5
            cut[b, b] -= 0.5
            cut[a, b] += 0.5
            cut[b, a] += 0.5
        np.testing.assert_allclose(cut @ np.ones(board.n), 0, atol=1e-14)
        rows[r] = dict(
            distance=r,
            second=list(map(int, second)),
            reference=min(e0, e1),
            ground_pair=e12,
            single_core_energies=[e0, e1],
            gain=min(e0, e1) - e12,
            largest_residual=max(err0, err1, err12),
        )
    fitted = 3 * rows[3]["gain"]
    if fitted <= 0:
        raise ValueError("No attractive calibration")
    for r, v in rows.items():
        forecast = fitted / r
        wait = np.pi / forecast
        dark = np.sin(wait * v["gain"] / 2) ** 2
        # Bound the error on the Ramsey amplitude jointly from the two eigen-residuals.
        v.update(
            coulomb_gain=forecast,
            wait=wait,
            dark_probability=float(dark),
            dark_count=int(rng.binomial(SHOTS, round(dark, 12))),
            cut_only_dark_count=0,
            free_amplitude_residual_bound=2 * wait * v["largest_residual"],
        )
    return dict(
        cells=board.cells,
        first=list(map(int, first)),
        fitted_coefficient=fitted,
        independent_three_link_calibration=measured_calibration(rows, board),
        samples=[rows[k] for k in sorted(rows)],
    )


def run():
    formed.native()
    rng = np.random.default_rng(311942)
    return dict(
        seed=311942, shots=SHOTS, studies=[scan(level, rng) for level in (0, 1)]
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        compare(result, json.loads((HERE / "mediator.json").read_text()), "mediator")
        print("Reproduced the indirectly mediated attraction and held distance scan.")


if __name__ == "__main__":
    main()
