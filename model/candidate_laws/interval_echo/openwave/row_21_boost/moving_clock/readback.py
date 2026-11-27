"""Estimates from port counts and the reference hand's three addresses only."""

import json
from pathlib import Path

import numpy as np
from scipy.special import betaincinv

HERE = Path(__file__).resolve().parent
FAIL = 0.01 / (2 * 30000)


def uncertainty(k, n):
    k = np.asarray(k)
    lower = np.zeros(k.shape)
    upper = np.ones(k.shape)
    ix = k > 0
    lower[ix] = betaincinv(k[ix], n - k[ix] + 1, FAIL)
    ix = k < n
    upper[ix] = betaincinv(k[ix] + 1, n - k[ix], 1 - FAIL)
    return np.maximum(k / n - lower, upper - k / n)


def center(counts):
    n = len(counts)
    return int(
        np.rint(
            n * np.angle(counts @ np.exp(2j * np.pi * np.arange(n) / n)) / (2 * np.pi)
        )
    )


def phase_read(rows, stations, rail, total):
    if np.any(rows.sum(axis=(2, 3)) != total):
        raise ValueError("Incomplete local position records")
    k = rows[
        np.arange(len(stations)), :, np.asarray(stations) % rows.shape[2], rail
    ].astype(float)
    return phase_bins(k, total)


def phase_bins(counts, total):
    k = np.asarray(counts, dtype=float)
    if np.any((k < 0) | (k > total)):
        raise ValueError("Invalid terminal-bin total")
    sig = ((k[:, 0] - k[:, 2]) + 1j * (k[:, 3] - k[:, 1])) / total
    errors = uncertainty(k, total)
    radius = np.hypot(errors[:, 0] + errors[:, 2], errors[:, 1] + errors[:, 3])
    if np.any(abs(sig) <= radius):
        return dict(resolved=False, upper_contrast=float(np.max(abs(sig) + radius)))
    tol = np.arcsin(radius / abs(sig))
    inc = np.angle(sig[1:] / sig[:-1])
    margin = np.pi - abs(inc) - tol[1:] - tol[:-1]
    return dict(
        resolved=bool(np.all(margin > 0)),
        elapsed=float(inc.sum()),
        phase_error=float(tol[0] + tol[-1]),
        min_branch_margin=float(min(margin)),
        min_contrast=float(np.min(abs(sig) - radius)),
    )


def report(data, bank):
    out = []
    zero = {}
    for row in data["rows"]:
        n, ex, h = row["n"], row["exchange"], row["h"]
        key = row["key"]
        pilot = bank[key + "_pilot"]
        if np.any(pilot.sum(axis=(2, 3)) != data["shots"]):
            raise ValueError("Missing band-pilot exits")
        x0, x1 = [center(t.sum(axis=(0, 2))) for t in pilot]
        if abs(x1 - x0) > n // 2:
            x1 += n if x1 < 0 else -n
        ts = np.rint(np.linspace(x0, x1, len(row["stations"]))).astype(int)
        if ts.tolist() != row["stations"]:
            raise ValueError("Held reader misses pilot endpoint")
        e, t, r = row["hand_positions"]
        if (e, t, r) != (t - abs(x1 - x0), len(ts) - 1, t + abs(x1 - x0)) or min(
            row["radar"]["forward_then_return"][-1], e
        ) <= 0:
            raise ValueError("Two-way material echo or clock failed")
        estimates = {
            s: phase_read(bank[key + "_" + s], ts, row["rail"], data["shots"])
            for s in ["packet", "plane", "incoherent"]
        }
        if h == 0:
            zero[n, ex] = estimates
        forecast = np.sqrt(e * r) / t
        res = {}
        for s in ["packet", "plane"]:
            val, ref = estimates[s], zero[n, ex][s]
            if val["resolved"] and ref["resolved"]:
                x, dx = val["elapsed"], val["phase_error"]
                y, dy = ref["elapsed"], ref["phase_error"]
                if y <= dy or x <= dx:
                    raise ValueError("Unresolved positive clock interval")
                val["ratio"] = float(x / y)
                val["interval"] = [
                    float((x - dx) / (y + dy)),
                    float((x + dx) / (y - dy)),
                ]
                val["passes"] = bool(
                    val["interval"][0] <= forecast <= val["interval"][1]
                )
            res[s] = val
        out.append(
            dict(
                key=key,
                h=h,
                d=x1 - x0,
                forecast=forecast,
                packet=res["packet"],
                plane=res["plane"],
                mixture=estimates["incoherent"],
            )
        )
    return out


if __name__ == "__main__":
    with np.load(HERE / "records.npz") as bank:
        print(
            json.dumps(
                report(json.loads((HERE / "records.json").read_text()), bank), indent=1
            )
        )
