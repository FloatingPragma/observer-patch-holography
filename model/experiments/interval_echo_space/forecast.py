"""Transfer a frozen hop law from short local surveys to held receiver records.

python -m model.experiments.interval_echo_space.forecast [--output NEW_DIRECTORY]
"""
from hashlib import sha256

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate as echo
from model.executable_model import federation
from model.executable_model.carrier import seam_digest
from .instrument import (HERE, WORDS, acquisition, click_matrix, decode, draw_panels, local_tables,
                         record_or_verify, relative_rms, slow_frame)

DWELL = 1 / 32
HOP = 0.5
SHOTS = 2**20


def enroll(board, sender):
    def expand(owners):
        return {endpoint // 12 for owner in owners for port in range(12 * owner, 12 * (owner + 1))
                for endpoint in board.neighbours[port]}

    immediate = sorted(expand([sender]) - {sender})
    outer = sorted(expand(immediate) - {*immediate, sender})
    if len(immediate) != 3:
        raise ValueError("The sender needs three distinct neighbours")
    return [sender, *immediate, *outer]


def infer_operator(survey):
    clicks = survey[:-1]
    off_diagonal = np.sort(clicks[np.triu_indices_from(clicks, k=1)])
    split = np.argmax(np.diff(off_diagonal))
    adjacency = clicks > off_diagonal[split:split + 2].mean()
    np.fill_diagonal(adjacency, False)
    degree = np.rint((1 - clicks.diagonal()) / (HOP * DWELL)**2).astype(int)
    if not np.array_equal(adjacency, adjacency.T) or np.any(degree < adjacency.sum(axis=1)):
        raise ValueError("Survey produced inconsistent total degrees or asymmetric contacts")
    return HOP * (np.diag(degree) - adjacency), degree, adjacency


def analyze(seen, predicted, frame, rng):
    forecast = acquisition(predicted, frame, WORDS)
    held = acquisition(seen, frame, WORDS)
    baseline_seen = acquisition(seen, frame, np.zeros((1, 3)))
    baseline_forecast = acquisition(predicted, frame, np.zeros((1, 3)))
    centered = click_matrix(held - baseline_seen)
    centered_forecast = click_matrix(forecast - baseline_forecast)
    u, singular, vh = np.linalg.svd(centered_forecast, full_matrices=False)
    reduced = (u[:, :2] * singular[:2]) @ vh[:2]
    counts = draw_panels(rng, held, SHOTS)
    decoded = decode(click_matrix(counts / SHOTS), click_matrix(forecast))
    relative = relative_rms(centered - centered_forecast, centered)
    relative_two = relative_rms(centered - reduced, centered)
    max_error = float(np.max(np.abs(held - forecast)))
    metrics = {
        "max_held_probability_error": max_error, "centered_relative_rms": relative,
        "rank_two_centered_relative_rms": relative_two, "decoded_of_eight": decoded,
        "passes_forecast": bool(max_error < 0.001 and relative < 0.05),
        "passes_rank_two_discrimination": bool(relative_two > 0.05),
        "passes_sampled": decoded >= 7,
    }
    arrays = {"forecast": forecast, "held": held, "counts": counts, "uniform_seen": baseline_seen,
              "uniform_forecast": baseline_forecast}
    return metrics, arrays


def measure():
    layout_rng = np.random.default_rng(5376)
    count_rng = np.random.default_rng(912844)
    sites, arrays = [], {}
    for level in range(3):
        board = federation(level)
        sender = board.cells // 2
        owners = enroll(board, sender)
        labels = np.array([12 * owner + layout_rng.permutation(12) for owner in owners])
        flat = labels.ravel()
        launch = np.zeros((board.n, len(flat)), dtype=complex)
        launch[flat, np.arange(len(flat))] = 1
        h = echo.generator(board)
        survey_state = expm_multiply(-1j * DWELL * h, launch)
        survey = local_tables(survey_state, flat[None, :])[0]
        inferred_h, degree, adjacency = infer_operator(survey)
        frame, edges = slow_frame(survey[:12, :12])
        prefix = f"n{board.cells}"
        arrays[prefix + "_survey"] = survey
        arrays[prefix + "_frame"] = frame
        site = {"cells": board.cells, "seams_sha256": seam_digest(board), "sender": sender, "patch": owners,
                "receivers": owners[1:4], "anonymous_local_port_order": (labels % 12).tolist(),
                "survey_native_edges_detected": int(adjacency.sum() // 2),
                "inferred_degree_multiplicities": {str(v): int(np.sum(degree == v)) for v in np.unique(degree)},
                "sender_edges": edges, "sender_rank": frame.shape[1], "dwells": {}}
        if frame.shape[1] == 3:
            for dwell in (0.75, 1.5, 3):
                # The forecast uses only the patch survey. The acquisition retains all external links.
                forecast_state = expm(-1j * dwell * inferred_h)[:, :12]
                forecast = local_tables(forecast_state, np.arange(12, 48).reshape(3, 12))
                held_state = expm_multiply(-1j * dwell * h, launch[:, :12])
                seen = local_tables(held_state, labels[1:4])
                metrics, records = analyze(seen, forecast, frame, count_rng)
                site["dwells"][str(dwell)] = metrics
                for key, value in records.items():
                    arrays[f"{prefix}_t{dwell}_{key}"] = value
        sites.append(site)
    return {"candidate": echo.NAME, "kappa": echo.KAPPA, "trials_per_local_setting": SHOTS,
            "protocol_sha256": sha256((HERE / "FORECAST_PROTOCOL.md").read_bytes()).hexdigest(),
            "sites": sites}, arrays


if __name__ == "__main__":
    receipt = record_or_verify(measure, __doc__, "forecast.json", "forecast_records.npz")
    panels = [r for site in receipt["sites"] for r in site["dwells"].values()]
    for name in ("passes_forecast", "passes_rank_two_discrimination", "passes_sampled"):
        print(f"{name} {sum(p[name] for p in panels)} of 9")
