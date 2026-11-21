"""Cross-cell detector tomography for the free one-load sector of Interval Echo.

python -m model.experiments.interval_echo_space [--output NEW_DIRECTORY]
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate as echo
from model.executable_model import federation
from model.executable_model.carrier import seam_digest

HERE = Path(__file__).resolve().parent
SHOTS = 2**22
WORDS = np.array(list(product((-1, 1), repeat=3))) / np.sqrt(3)
TRAINING = np.vstack((np.zeros(3), np.eye(3), -np.eye(3)))


def slow_frame(home_clicks):
    """Infer the cell incidence and slow space from the short-exposure home table."""
    off_diagonal = np.sort(home_clicks[np.triu_indices(len(home_clicks), 1)])
    split = np.argmax(np.diff(off_diagonal))
    adjacency = home_clicks > off_diagonal[split:split + 2].mean()
    np.fill_diagonal(adjacency, False)
    if not np.array_equal(adjacency, adjacency.T):
        raise ValueError("Detected home incidence is asymmetric")
    spectrum, vectors = np.linalg.eigh(np.diag(adjacency.sum(axis=1)) - adjacency)
    slow = spectrum[spectrum > 1e-8][0]
    selected = vectors[:, np.abs(spectrum - slow) < 1e-8]
    projector = selected @ selected.T
    basis = []
    for vector in projector.T:
        residual = vector.copy()
        for v in basis:
            residual -= (v @ residual) * v
        if np.linalg.norm(residual) > 1e-8:
            basis.append(residual / np.linalg.norm(residual))
    return np.array(basis).T, int(adjacency.sum() // 2)


def local_tables(states, cells):
    clicks = np.abs(states[cells])**2
    quiet = 1 - clicks.sum(axis=1, keepdims=True)
    if np.min(quiet) < -1e-12:
        raise ArithmeticError("Local click outcomes exceed unit probability")
    return np.concatenate((clicks, np.maximum(quiet, 0)), axis=1)


def acquisition(table, frame, directions):
    # A source draws a sharp port from the indicated mixture, then an isolated
    # trial in each receiver cell ends in one of its twelve ports or quiet.
    weights = 1 / 12 + directions @ frame.T / 24
    if np.min(weights) < 0 or not np.allclose(weights.sum(axis=1), 1, rtol=0, atol=1e-12):
        raise ValueError("Source preparation is outside the probability simplex")
    return np.einsum("ws,ros->wro", weights, table)


def draw_panels(rng, probabilities, shots=SHOTS):
    series = [[rng.multinomial(shots, p / p.sum()) for p in receivers] for receivers in probabilities]
    return np.array(series, dtype=np.int64)


def fit(train):
    baseline = train[0]
    response = (train[1:4] - train[4:7]) / 2
    return baseline, response


def click_matrix(tables):
    return tables[:, :, :12].reshape(len(tables), -1)


def relative_rms(residual, centered_signal):
    return float(np.linalg.norm(residual) / np.linalg.norm(centered_signal))


def decode(observed, candidates):
    distances = np.linalg.norm(observed[:, None] - candidates, axis=2)
    labels = distances.argmin(axis=1)
    return int(np.sum(labels == np.arange(len(labels))))


def panel(table, frame, count_rng):
    train = acquisition(table, frame, TRAINING)
    held = acquisition(table, frame, WORDS)
    baseline, responses = fit(train)
    pred = baseline + np.einsum("wi,iro->wro", WORDS, responses)
    full_error = np.max(np.abs(pred - held))
    contrasts = click_matrix(responses)
    u, singular, vh = np.linalg.svd(contrasts, full_matrices=False)
    rank_two = (u[:, :2] * singular[:2]) @ vh[:2]
    centered = click_matrix(held - baseline)
    error_two = relative_rms(centered - WORDS @ rank_two, centered)
    sigma_ratio = singular[-1] / singular[0]
    singles = [np.linalg.svd(responses[:, r, :12], compute_uv=False) for r in range(table.shape[0])]
    # Record a physically unmodulated source for each intended held word.
    sham = acquisition(table, frame, np.zeros_like(WORDS))
    sham_contrast = float(np.max(np.abs(sham - baseline)))
    train_counts = draw_panels(count_rng, train)
    held_counts = draw_panels(count_rng, held)
    sham_counts = draw_panels(count_rng, sham)
    sampled_base, sampled_responses = fit(train_counts / SHOTS)
    sampled_contrasts = click_matrix(sampled_responses)
    su, ss, svh = np.linalg.svd(sampled_contrasts, full_matrices=False)
    sampled_two = (su[:, :2] * ss[:2]) @ svh[:2]
    predictions = sampled_base[:, :12].ravel() + WORDS @ sampled_contrasts
    reduced = sampled_base[:, :12].ravel() + WORDS @ sampled_two
    recorded = click_matrix(held_counts / SHOTS)
    decoded = decode(recorded, predictions)
    metrics = {
        "singular_values": singular.tolist(), "smallest_to_largest": float(sigma_ratio),
        "each_recipient_smallest_to_largest": [float(s[-1] / s[0]) for s in singles],
        "held_exact_max_error": float(full_error), "rank_two_relative_rms": error_two,
        "held_decoded_of_eight": decoded, "sham_decoded_of_eight": decode(click_matrix(sham_counts / SHOTS), predictions),
        "sham_exact_max_contrast": sham_contrast,
        "sampled_full_click_rms": float(np.sqrt(np.mean((recorded - predictions)**2))),
        "sampled_rank_two_click_rms": float(np.sqrt(np.mean((recorded - reduced)**2))),
        "held_per_recipient_min_max_detection": [float((1 - held[:, :, -1]).min()), float((1 - held[:, :, -1]).max())],
        "passes_exact": bool(full_error < 1e-10 and sigma_ratio >= 0.05 and error_two > 0.05 and sham_contrast < 1e-14),
        "passes_sampled": decoded >= 7,
    }
    arrays = {"train": train, "held": held, "train_counts": train_counts, "held_counts": held_counts, "sham_counts": sham_counts}
    return metrics, arrays


def measure():
    names_rng = np.random.default_rng(73109)
    count_rng = np.random.default_rng(92017)
    arrays, sites = {}, []
    for level in range(3):
        board = federation(level)
        h = echo.generator(board)
        for sender in (0, board.cells // 3, 2 * board.cells // 3):
            neighbors = sorted({b // 12 for a in range(12 * sender, 12 * (sender + 1))
                                for b in board.neighbours[a] if b // 12 != sender})
            if len(neighbors) != 3:
                raise ValueError("Three distinct neighbouring patches required")
            cells = np.array([12 * c + names_rng.permutation(12) for c in [sender, *neighbors]])
            launch = np.zeros((board.n, 12), dtype=complex)
            launch[cells[0], np.arange(12)] = 1
            short = expm_multiply(-1j * h / 8, launch)
            home = local_tables(short, cells[:1])[0]
            frame, edges = slow_frame(home[:12])
            prefix = f"n{board.cells}_s{sender}"
            arrays[prefix + "_home"] = home
            arrays[prefix + "_frame"] = frame
            site = {"cells": board.cells, "seams_sha256": seam_digest(board), "sender": sender, "receivers": neighbors,
                    "anonymous_local_port_order": (cells % 12).tolist(), "local_rank": frame.shape[1],
                    "inferred_edges": edges, "dwells": {}}
            if frame.shape[1] == 3:
                for dwell in (2, 4):
                    evolved = expm_multiply(-1j * dwell * h, launch)
                    table = local_tables(evolved, cells[1:])
                    metrics, records = panel(table, frame, count_rng)
                    site["dwells"][str(dwell)] = metrics
                    for label, value in records.items():
                        arrays[f"{prefix}_t{dwell}_{label}"] = value
            sites.append(site)
    result = {"candidate": echo.NAME, "kappa": echo.KAPPA, "trials_per_local_setting": SHOTS,
              "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
              "sites": sites}
    return result, arrays


def normalized(value):
    if isinstance(value, dict):
        return {k: normalized(v) for k, v in value.items()}
    if isinstance(value, list):
        return [normalized(v) for v in value]
    if isinstance(value, float):
        return round(value, 12)
    return value


def record_or_verify(measurement, description, summary_file, record_file):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--output", type=Path, help="write receipt and local records to a fresh destination")
    args = parser.parse_args()
    result, arrays = measurement()
    result = normalized(result)
    if args.output:
        args.output.mkdir(exist_ok=True, parents=True)
        paths = [args.output / summary_file, args.output / record_file]
        if any(p.exists() for p in paths):
            parser.error("the output measurement exists")
        np.savez_compressed(paths[1], **arrays)
        paths[0].write_text(json.dumps(result, indent=2) + "\n")
    else:
        if result != json.loads((HERE / summary_file).read_text()):
            raise SystemExit("Published summary differs from measurement")
        with np.load(HERE / record_file, allow_pickle=False) as stored:
            if set(stored.files) != set(arrays):
                raise SystemExit("Raw-record keys differ")
            for k, values in arrays.items():
                saved = stored[k]
                if values.dtype.kind in "iu":
                    np.testing.assert_array_equal(values, saved)
                else:
                    np.testing.assert_allclose(values, saved, rtol=0, atol=1e-12)
    return result


def main():
    result = record_or_verify(measure, __doc__, "result.json", "records.npz")
    panels = [p for site in result["sites"] for p in site["dwells"].values()]
    print(f"{len(panels)} of 18 panels acquired, {sum(p['passes_exact'] for p in panels)} exact gates, "
          f"{sum(p['passes_sampled'] for p in panels)} finite-count gates.")


if __name__ == "__main__":
    main()
