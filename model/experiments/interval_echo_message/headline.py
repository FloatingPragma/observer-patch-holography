"""Reconstruct the channel and score the retained counts without producer imports."""
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent


def rebuild(station, exposure):
    with np.load(HERE.parents[2] / "model/executable_model/carriers" /
                 f"FEDERATION_LEVEL{station['level']}.npz", allow_pickle=False) as data:
        edges = data["edges"]
    n = 12 * station["cells"]
    a, b = edges.T
    rows = np.arange(len(edges))
    incidence = sparse.coo_matrix((np.r_[np.ones(len(a)), -np.ones(len(b))],
                                   (np.r_[rows, rows], np.r_[a, b])), shape=(len(a), n)).tocsr()
    exchange = (incidence.T @ incidence) * .5
    degree = np.bincount(edges.ravel(), minlength=n)
    # Count the unit squared contrast on incident physical and reference seams.
    references = 6 - degree
    if np.any(references < 0) or not np.array_equal(exchange.diagonal() * 2 + references, np.full(n, 6)):
        raise ValueError("degree correction")
    state = np.zeros((n, 2), dtype=complex)
    state[station["ports"], [0, 1]] = 1
    # Drop the scalar phase. Taylor recurrence, ||L/2|| <= 6 by Gershgorin.
    steps = max(1, int(np.ceil(exposure * 8)))
    step = exchange * (-1j * exposure / steps)
    for _ in range(steps):
        term = state.copy()
        result = state.copy()
        for order in range(1, 19):
            term = (step @ term) / order
            result += term
        state = result
    mass = abs(state)**2
    return np.vstack((mass[:12], mass[12:].sum(axis=0))).T


def check(protocol, receipt):
    gains = []
    worst = 0.
    if len(receipt["stations"]) != len(protocol["stations"]):
        raise ValueError("station count")
    for station, panel in zip(protocol["stations"], receipt["stations"], strict=True):
        if panel["cells"] != station["cells"]:
            raise ValueError("wrong board")
        for s, row in zip(protocol["exposures"], panel["windows"], strict=True):
            if s != row["exposure"]:
                raise ValueError("wrong window")
            p = rebuild(station, s)
            actual = np.asarray(row["probabilities"])
            if actual.shape != (2, 13):
                raise ValueError("missing outcome")
            error = float(np.max(abs(p - actual)))
            np.testing.assert_allclose(actual, p, atol=protocol["numeric_tolerance"], rtol=0)
            worst = max(worst, error)
            oracle = .5 + abs(p[0] - p[1]).sum() / 4
            np.testing.assert_allclose(row["detection"], p[:, :12].sum(axis=1), rtol=0,
                                       atol=protocol["numeric_tolerance"])
            if (abs(oracle - row["optimal"]) > protocol["numeric_tolerance"]
                    or abs(2 * oracle - 1 - row["tv"]) > protocol["numeric_tolerance"]):
                raise ValueError("oracle")
            learn = np.array(row["calibration_counts"])
            if (learn.shape != (2, 13) or learn.dtype.kind not in "iu" or learn.min() < 0
                    or np.any(learn.sum(axis=1) != protocol["calibration_per_bit"])):
                raise ValueError("calibration enrolment")
            table = (learn[1] > learn[0]).astype(int)
            if table.tolist() != row["lookup"]:
                raise ValueError("calibration")
            for arm in ("native", "no_encoding"):
                result = row[arm]
                c = np.array(result["counts"])
                n = protocol["trials"]
                if (c.shape != (2, 13) or c.dtype.kind not in "iu" or c.min() < 0
                        or np.any(c.sum(axis=1) != n // 2)):
                    raise ValueError("enrolment")
                correct = int(c[table, np.arange(13)].sum())
                if (correct != result["correct"] or correct / n != result["accuracy"]
                        or result["total"] != n or n - correct != result["errors"]
                        or c[:, -1].sum() != result["quiet"]
                        or c[table[-1], -1] != result["quiet_correct"] or result["abstentions"]):
                    raise ValueError("retained decisions")
                radius = math.sqrt(math.log(2 * protocol["comparisons"] / protocol["alpha"]) / (2 * n))
                np.testing.assert_allclose(result["interval"],
                                           [max(0., correct / n - radius), min(1., correct / n + radius)],
                                           atol=1e-14, rtol=0)
            gains.append((station["cells"], s, row["native"]["correct"], round(float(oracle), 8)))
    return worst, gains


if __name__ == "__main__":
    error, rows = check(json.loads((HERE / "protocol.json").read_text()),
                        json.loads((HERE / "receipt.json").read_text()))
    print("cells, exposure, correct, classical optimum:", *rows, sep="\n")
    print("Maximum absolute probability discrepancy:", error)
