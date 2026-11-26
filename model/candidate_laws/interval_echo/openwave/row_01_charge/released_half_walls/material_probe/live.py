"""One freely evolving word amplitudes and two mobile meters, including backreaction."""

import argparse
import json
from functools import cache
from itertools import product
from math import log, sqrt
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from .. import ring
from .flight import ANALYZERS, mount, verify_plates
from .measure import HERE, draw, source

LIVE_SHOTS = 2**28
LIVE_SEED = 26011263
HOLDS = (1 / 16, 1 / 4)
MODES = ("live", "detached", "stopped")


@cache
def resident(k):
    h = ring.generator(14, 1, k)
    start = np.zeros(h.shape[0])
    start[np.searchsorted(ring.basis(14)[0], ring.launch(14))] = 1
    a = -3j * h
    state = expm_multiply(a, start, traceA=a.diagonal().sum())
    np.testing.assert_allclose(
        abs(state) ** 2, ring.released(14, 1, k)[4], atol=2e-12, rtol=0
    )
    return state


def combined(k, src, probe, mode):
    mounts, _ = mount(0, 14)
    bonds = tuple(m[0] for m in mounts)
    words, gate, _ = source(14, 1, src, k, bonds)
    sums = [words[:, j] + words[:, (j + 1) % 14] for j in bonds]
    h = ring.generator(14, 1, k)
    if mode == "stopped":
        h = sparse.diags(h.diagonal())
    exchange = np.array([[0, -0.5], [-0.5, 0]])
    two_meter = sparse.kron(sparse.csr_matrix(exchange), np.eye(2)) + sparse.kron(
        sparse.eye(2), exchange
    )
    joint = sparse.kron(h, sparse.eye(4)) + sparse.kron(
        sparse.eye(len(words)), two_meter
    )
    occup = np.array([1, 1 - probe])
    if mode in ("live", "stopped"):
        # For each actual source-to-tip contact the ring already carries κ n².
        # Probe-only squared loads sum to a constant on both equal-degree rails.
        # The remaining term is −2 κ n(load on sensing tip) in either sector.
        diag = (
            -2
            * k
            * (
                sums[0][:, None] * np.repeat(occup, 2)
                + sums[1][:, None] * np.tile(occup, 2)
            )
        )
        joint += sparse.diags(diag.ravel(), dtype=float)
    elif mode == "cut":
        # Remove both squared-load self terms, four full physical seams in all.
        ring_loss = sum(words[:, j] ** 2 + words[:, (j + 1) % 14] ** 2 for j in bonds)
        tip_loss = 2 * (np.repeat(occup**2, 2) + np.tile(occup**2, 2))
        joint -= sparse.diags(
            (k * (ring_loss[:, None] + tip_loss)).ravel(), dtype=float
        )
    elif mode != "detached":
        raise ValueError("Unexpected physical probe treatment")
    return joint.tocsr(), gate


def bins(after, gate, read):
    q = np.array([[1, 1j], [1j, 1]]) / np.sqrt(2)
    if read == "phase":
        after = after @ np.kron(q, q).T
    p = abs(after.reshape(-1, 2, 2)) ** 2
    # State index zero places the local defect in the first rail, click one.
    return np.stack(
        [p[~gate].sum(axis=0)[::-1, ::-1], p[gate].sum(axis=0)[::-1, ::-1]]
    ).ravel()


def acquire(modes=MODES, seed=LIVE_SEED):
    verify_plates()
    records = []
    for k, src, s in product((8, 4), (-1, 1), (-1, 1)):
        psi = np.kron(resident(k), np.ones(4) / 2)
        for mode in modes:
            h, gate = combined(k, src, s, mode)
            for i, t in enumerate(HOLDS):
                a = -1j * t * h
                after = expm_multiply(a, psi, traceA=a.diagonal().sum()).reshape(-1, 4)
                for read in ANALYZERS:
                    key = [k, src, s, mode, i, read]
                    p = bins(after, gate, read)
                    counts = draw([seed, *key], p, size=LIVE_SHOTS)
                    records.append([key, counts])
    return dict(seed_key_prefix=seed, shots=LIVE_SHOTS, banks=records)


def bound(hits, n):
    if n == 0:
        return np.array([0.0, 1.0])
    delta = sqrt(log(2 * 2048 / 0.01) / (2 * n))
    f = hits / n
    return np.array([f - delta, f + delta])


def separate(x, y):
    return bool(x[0] > y[1] or y[0] > x[1])


def load_banks(raw, modes, seed):
    if raw["shots"] != LIVE_SHOTS or raw["seed_key_prefix"] != seed:
        raise ValueError("Unexpected live seed or exposure")
    rows = {}
    for key, v in raw["banks"]:
        key = tuple(key)
        if (
            key in rows
            or len(v) != 8
            or sum(v) != LIVE_SHOTS
            or any(type(x) is not int or x < 0 for x in v)
        ):
            raise ValueError("Invalid probe cohort")
        rows[key] = np.array(v).reshape(2, 2, 2)
    required = set(
        product((8, 4), (-1, 1), (-1, 1), modes, range(len(HOLDS)), ANALYZERS)
    )
    if required != rows.keys():
        raise ValueError("Missing or unknown live setting")
    return rows


def decode(raw):
    rows = load_banks(raw, MODES, LIVE_SEED)
    report = {}
    for k in (8, 4):
        changed = 0
        wall_total = 0
        clicks_changed = 0
        clicks_total = 0
        examples = []
        for src, s, i, read in product((-1, 1), (-1, 1), range(len(HOLDS)), ANALYZERS):

            def accepted(mode):
                return rows[k, src, s, mode, i, read][1]

            live, detached, stopped = (accepted(mode) for mode in MODES)
            if read == "rail":
                changed += int(
                    separate(
                        bound(live.sum(), LIVE_SHOTS), bound(detached.sum(), LIVE_SHOTS)
                    )
                )
                wall_total += 1
                if k == 8 and src == s == 1:
                    examples.append(
                        dict(
                            dwell=HOLDS[i],
                            wall_live=int(live.sum()),
                            wall_detached=int(detached.sum()),
                            wall_stopped=int(stopped.sum()),
                            left_click_live=int(live[1].sum()),
                            left_click_stopped=int(stopped[1].sum()),
                        )
                    )
            for connected_hits, stopped_hits in (
                (live[1].sum(), stopped[1].sum()),
                (live[:, 1].sum(), stopped[:, 1].sum()),
            ):
                clicks_changed += int(
                    separate(
                        bound(connected_hits, live.sum()),
                        bound(stopped_hits, stopped.sum()),
                    )
                )
                clicks_total += 1
        report[k] = dict(
            direct_wall_rates_resolving_backaction=[changed, wall_total],
            held_single_probe_rates_resolving_live_from_stopped=[
                clicks_changed,
                clicks_total,
            ],
            examples=examples,
        )
    return report


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    if args.output:
        raw = acquire()
        args.output.write_text(json.dumps(raw, separators=(",", ":")) + "\n")
    else:
        raw = json.loads((HERE / "live_counts.json").read_text())
    print(json.dumps(decode(raw), indent=2))
