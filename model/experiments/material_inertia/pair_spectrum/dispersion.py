"""Read the energy transfer of graph-calibrated centre-density gratings."""

import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np

from . import measure as flight

HERE = Path(__file__).resolve().parent
COUNT = 128
STEP = 131072
SHOTS = 2**36
BIN_WIDTH = 2 * np.pi / (COUNT * STEP)


def geometry(board):
    graph = nx.Graph((a // 12, b // 12) for a, b in board.external)
    lap = nx.laplacian_matrix(graph, nodelist=range(board.cells)).toarray()
    vals, v = np.linalg.eigh(lap)
    shells = np.unique(vals.round(9))[1:4]
    result = []
    for lam in shells:
        block = v[:, abs(vals - lam) < 2e-8]
        p = block @ block.T
        np.testing.assert_allclose(lap @ p, lam * p, atol=1e-9)
        result.append((float(vals[abs(vals - lam) < 2e-8].mean()), p))
    return result


def timing_histogram(gaps, compressed, g):
    # Each timing bin adds amplitudes only for the same final material state.
    times = np.arange(COUNT) * STEP
    last = times[-1]
    before = g[:, None] * np.exp(-1j * gaps[: len(g), None] * (last - times))
    states = (compressed @ before) * np.exp(-1j * gaps[:, None] * times)
    bins = np.fft.ifft(states, axis=1)
    probabilities = abs(bins) ** 2
    np.testing.assert_allclose(
        probabilities.sum(), np.mean(np.linalg.norm(states, axis=0) ** 2), atol=1e-14
    )
    return probabilities


def sample(n, p, rng):
    p = np.asarray(p, float)
    np.testing.assert_allclose(p.sum(), 1, atol=3e-12, rtol=0)
    p = np.maximum(p, 0).round(13)
    return rng.multinomial(n, p / p.sum())


def read(counts):
    counts = np.array(counts)
    peak = int(np.argmax(counts))
    indexes = np.mod(np.arange(peak - 1, peak + 2), COUNT)
    share = float(counts[indexes].sum() / counts.sum())
    signed = (peak + COUNT // 2) % COUNT - COUNT // 2
    return {
        "signed_peak_bin": int(signed),
        "peak_gap": float(signed * BIN_WIDTH),
        "counts_near_peak_fraction": share,
        "has_majority_peak": share > 0.5,
    }


def row(spatial, lam, eig, basis, g, admission, projector_error, residual, x, y, rng):
    amp = (spatial[x // 12] + spatial[y // 12]) / 2
    compressed = basis.T @ (amp[:, None] * basis[:, :2])
    probs = timing_histogram(eig - eig[0], compressed, g)
    same, new = probs[:2].sum(axis=0), probs[2:].sum(axis=0)
    fail, passed = sample(SHOTS, [1 - admission, admission], rng)
    bins = sample(passed, np.r_[same, new, 1 - same.sum() - new.sum()], rng)
    # Same source and gratings, changed local law with all native generator
    # rates halved in flight. Prep and the final spectral projectors stay fixed.
    slower = timing_histogram((eig - eig[0]) / 2, compressed, g)
    # Compare before conditioning. Each free history has the same total time.
    # With an exact 2-line initial projector, the low (24) final cluster has a
    # separate gap-error bound after projecting the grating back into the band.
    duration = (COUNT - 1) * STEP
    ground_projector_error = residual / min(
        (eig[2] - eig[1]) / 2 - residual, flight.CUT - eig[-1]
    )
    amplitude_error = (
        2 * ground_projector_error / np.sqrt(admission)
        + 2 * projector_error
        + ground_projector_error
        + duration * residual
    )
    cond_error = 2 * amplitude_error / np.sqrt(new.sum())
    inelastic = bins[COUNT : 2 * COUNT]
    max_mass = float(
        np.max(inelastic + np.roll(inelastic, 1) + np.roll(inelastic, -1))
        / inelastic.sum()
    )
    sampling_error = float(
        np.sqrt(np.log(2 * 9 * COUNT / 0.01) / (2 * inelastic.sum()))
    )
    return {
        "any_three_bins_max_fraction": max_mass,
        "any_three_bins_upper_99_with_numerical_error": max_mass
        + sampling_error
        + cond_error,
        "spatial_laplacian": lam,
        "owner_profile": spatial.tolist(),
        "original_enrollments": SHOTS,
        "rest_preparation_failed": int(fail),
        "passed_rest": int(passed),
        "elastic_counts": bins[:COUNT].tolist(),
        "inelastic_counts": bins[COUNT : 2 * COUNT].tolist(),
        "grating_or_final_energy_filter_failed": int(bins[-1]),
        "inelastic_probability_given_rest": float(new.sum()),
        "numerical_conditional_probability_bound": float(cond_error),
        "read": read(bins[COUNT : 2 * COUNT]),
        "half_rate_expectation_conditional": (
            slower[2:].sum(axis=0) / slower[2:].sum()
        ).tolist(),
    }


def fit(rows):
    first = rows[0]
    b = first["read"]["signed_peak_bin"]
    a = b * BIN_WIDTH / first["spatial_laplacian"]
    held = []
    for r in rows[1:]:
        ratio = r["spatial_laplacian"] / first["spatial_laplacian"]
        pred = b * ratio
        seen = r["read"]["signed_peak_bin"]
        allowed = ratio + 1
        distance = abs((seen - pred + COUNT / 2) % COUNT - COUNT / 2)
        held.append(
            {
                "spatial_laplacian": r["spatial_laplacian"],
                "predicted_bin": float(pred),
                "observed_bin": seen,
                "tolerance_bins": float(allowed),
                "passes": bool(
                    distance <= allowed
                    and r["read"]["has_majority_peak"]
                    and first["read"]["has_majority_peak"]
                ),
            }
        )
    return {"first_A": float(a), "held": held}


def scalar_reader_control(board, shells):
    graph = nx.Graph((a // 12, b // 12) for a, b in board.external)
    lap = nx.laplacian_matrix(graph, nodelist=range(board.cells)).toarray()
    vals, v = np.linalg.eigh(lap)
    coefficient = 4 * BIN_WIDTH
    rows = []
    rng = np.random.default_rng(381042)
    for lam, projector in shells:
        f = projector[:, 0] / max(abs(projector[:, 0]))
        compressed = v.T @ (f[:, None] * v[:, :1])
        probs = timing_histogram(coefficient * (vals - vals[0]), compressed, np.ones(1))
        bins = sample(SHOTS, np.r_[probs[1:].sum(axis=0), 1 - probs[1:].sum()], rng)
        rows.append(
            {
                "spatial_laplacian": lam,
                "inelastic_bin_counts": bins[:-1].tolist(),
                "other": int(bins[-1]),
                "read": read(bins[:-1]),
            }
        )
    return {"declared_coefficient": coefficient, "shells": rows, "fit": fit(rows)}


def run():
    b = flight.federation(0)
    x, y, _, _, h, cost = flight.defect.generators(b)
    reduced = flight.PairReduction(h, cost)
    eig, basis, residual, _ = reduced.lowest(h, 24)
    file = json.loads((HERE / "results.json").read_text())
    np.testing.assert_allclose(eig, file["lines"], rtol=0, atol=2e-12)
    cut = (eig[1] + eig[2]) / 2
    # An independent inertia count verifies the subspace used for the rest filter.
    certificate = reduced.inertia(cut)
    if certificate["upper_lines_below"] != 2:
        raise ValueError("Unresolved rest doublet")
    stationary = []
    shells = geometry(b)
    rng = np.random.default_rng(381041)
    for launch in file["flights"]:
        if launch["admission"] < 0.001:
            continue
        src = flight.defect.ordered_index(*launch["ports"], b.n)
        g = basis[src, :2]
        admission = g @ g
        rows = []
        for lam, p in shells:
            f = p[:, launch["owner"]]
            rows.append(
                row(
                    f / abs(f).max(),
                    lam,
                    eig,
                    basis,
                    g / np.sqrt(admission),
                    admission,
                    max(
                        file["cluster_projector_error"],
                        residual / (flight.CUT - eig[-1]),
                    ),
                    residual,
                    x,
                    y,
                    rng,
                )
            )
        stationary.append(
            {
                "launch": launch["owner"],
                "ports": launch["ports"],
                "rest_doublet_admission": float(admission),
                "shells": rows,
                "fit": fit(rows),
            }
        )
    return {
        "seed": 381041,
        "timing_bins": COUNT,
        "tick_spacing": STEP,
        "gap_per_bin": BIN_WIDTH,
        "energy_range": [float(eig[0]), float(eig[-1])],
        "rest_cut_inertia": certificate,
        "generator_residual": float(residual),
        "sources": stationary,
        "scalar_reader_control": scalar_reader_control(b, shells),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(
            json.loads(json.dumps(result)),
            json.loads((HERE / "dispersion.json").read_text()),
            "pair dispersion",
        )
        print(
            "Reproduced the geometry-only calibration and three held dispersion studies."
        )
