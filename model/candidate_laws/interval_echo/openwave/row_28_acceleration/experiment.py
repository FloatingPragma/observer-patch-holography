"""Routed material distances during local sharp-record settling."""

import argparse
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.interval_echo_consensus.asynchronous import priority

HERE = Path(__file__).resolve().parent
EPOCHS = (0, 2, 4, 6, 8, 12, 16, 24, 32, 40, 48)
SHOTS = 2**24
NSEEDS = 8


def calibrate():
    law = candidate.native_law
    h0, h1 = law.local_term_parts()
    h = h0 + candidate.KAPPA * h1
    flip = np.zeros((law.N_SRC, law.N_SRC))
    for i, (a, b) in enumerate(law.SRC):
        flip[law.SRC.index((b, a)), i] = 1
    eligible = np.diag([0, *[int(abs(d) == 2) for _, d in law.LABELS]])
    changed = h + np.kron(eligible, (np.eye(law.N_SRC) - flip) / 2)
    us = [expm(-1j * np.pi * op) for op in (h, changed)]
    clock = expm(-1j * h)
    basis = np.eye(law.DIM)
    np.testing.assert_allclose(
        abs(clock[law.SRC.index((0, 1)), law.SRC.index((1, 0))]) ** 2,
        np.sin(0.5) ** 2,
        atol=1e-12,
    )
    tests = []
    for sign in (-2, 0, 2):
        if sign:
            i = law.SRC.index((1 + sign // 2, 1 - sign // 2))
            mem = (law.echo() @ basis[:, i]).reshape(-1, law.N_SRC)[
                :, law.SRC.index((1, 1))
            ]
        else:
            mem = np.eye(1 + law.N_LAB)[0]
        for tag in (1, 2):
            for direction in (0, 1):
                start = (tag, 0) if direction == 0 else (0, tag)
                source = np.kron(mem, np.eye(law.N_SRC)[law.SRC.index(start)])
                reads = []
                for u in us:
                    out = (u @ source).reshape(-1, law.N_SRC)
                    prob = np.linalg.norm(out[:, law.SRC.index(start[::-1])]) ** 2
                    reads.append(round(float(prob), 12))
                assert reads == [float(sign == 0), 1.0]
                tests.append([sign, tag, direction, reads])
    return tests


def clocks():
    rng = np.random.default_rng(280101)
    counts = rng.binomial(SHOTS, np.sin(0.5) ** 2, 48)
    p = counts / SHOTS
    eps = np.sqrt(np.log(2 * len(p) / 0.01) / (2 * SHOTS))

    def inv(v):
        return 2 * np.arcsin(np.sqrt(v))

    return dict(
        shots=SHOTS,
        counts=counts.tolist(),
        waits=inv(p).tolist(),
        lower=inv(p - eps).tolist(),
        upper=inv(p + eps).tolist(),
    )


def geometry(n, edges, source, active=None):
    a, b = edges.T if active is None else edges[active].T
    graph = sparse.csr_matrix(
        (np.ones(2 * len(a)), (np.r_[a, b], np.r_[b, a])), shape=(n, n)
    )
    ds = shortest_path(graph, directed=False, unweighted=True, indices=source)
    return np.where(np.isfinite(ds), ds, -1).astype(np.int32)


@lru_cache(maxsize=3)
def layout(level):
    board = federation(level)
    edges = np.array(board.seams)
    stations = [0]
    raw = [geometry(board.n, edges, 0)]
    while len(stations) < 8:
        stations.append(int(np.argmax(np.min(raw, axis=0))))
        raw.append(geometry(board.n, edges, stations[-1]))
    return (
        board.n,
        edges,
        np.array([priority(e) for e in board.seams]),
        stations,
        np.array(raw),
    )


def history(n, edges, ranks, initial, probability, rng):
    word = initial.copy()
    a, b = edges.T
    commit = np.zeros(len(a), dtype=np.int16)
    sign = np.zeros(len(a), dtype=np.int8)
    epoch = 0
    while True:
        ready = (commit == 0) & (abs(word[a] - word[b]) == 2)
        if not ready.any():
            break
        epoch += 1
        rank = np.full(n, ranks.max() + 1, dtype=int)
        np.minimum.at(rank, a[ready], ranks[ready])
        np.minimum.at(rank, b[ready], ranks[ready])
        enabled = ready & (ranks == rank[a]) & (ranks == rank[b])
        accepted = enabled & (rng.random(len(a)) < probability)
        assert len(np.unique(edges[accepted])) == 2 * accepted.sum()
        commit[accepted] = epoch
        sign[accepted] = word[a[accepted]] - word[b[accepted]]
        word[a[accepted]] = word[b[accepted]] = 1
    return commit, sign, epoch


def trial(level, seed, probability, clock):
    n, edges, ranks, stations, d0 = layout(level)
    rng = np.random.default_rng(seed)
    init = 2 * rng.integers(0, 2, n, dtype=np.int8)
    commit, sign, settled = history(n, edges, ranks, init, probability, rng)
    df = geometry(n, edges, stations, commit == 0)
    masks = [
        (d0 > 0)
        & (3 * d0 <= (i + 1) * d0.max(axis=1)[:, None])
        & (3 * d0 > i * d0.max(axis=1)[:, None])
        for i in range(3)
    ]
    records = []
    for epoch in EPOCHS:
        ds = (
            df
            if epoch >= settled
            else geometry(n, edges, stations, (commit == 0) | (commit > epoch))
        )
        bins = []
        for mask in masks:
            visible = mask & (ds >= 0)
            final = mask & (df >= 0)
            if np.any(ds[final] < d0[final]):
                raise ValueError("range shortened on a deletion-only arm")
            bins.append(
                dict(
                    enrolled=int(mask.sum()),
                    observed=int(visible.sum()),
                    sum_range=int(ds[visible].sum()),
                    initial_sum_for_observed=int(d0[visible].sum()),
                    cohort_count=int(final.sum()),
                    cohort_sum=int(ds[final].sum()),
                    cohort_initial=int(d0[final].sum()),
                )
            )
        records.append(
            dict(
                epoch=epoch,
                writes=int(np.sum((commit > 0) & (commit <= epoch))),
                bins=bins,
            )
        )
    times = np.r_[0, np.cumsum(clock["waits"])]

    def estimate(index):
        if not records[0]["bins"][index]["cohort_initial"]:
            return None
        scales = np.array(
            [
                r["bins"][index]["cohort_sum"] / r["bins"][index]["cohort_initial"]
                for r in records
            ]
        )
        rate = np.log(scales[2]) / times[4]
        tlo = np.r_[0, np.cumsum(clock["lower"])]
        thi = np.r_[0, np.cumsum(clock["upper"])]
        predictions = []
        for t in (8, 12, 16):
            val = scales[EPOCHS.index(t)]
            lo = scales[2] ** (tlo[t] / thi[4])
            hi = scales[2] ** (thi[t] / tlo[4])
            predictions.append(
                dict(
                    t=t,
                    scale=val,
                    predicted=np.exp(rate * times[t]),
                    rejected=bool(val < lo / 1.01 or val > 1.01 * hi),
                )
            )
        slopes = np.diff(scales[:3]) / np.diff(times[[0, 2, 4]])
        acc = 2 * (slopes[1] - slopes[0]) / times[4] / scales[1]
        return dict(
            first_h=rate,
            initial_acceleration=acc,
            scales=scales.tolist(),
            held=predictions,
        )

    return (
        dict(
            level=level,
            seed=seed,
            attempt_probability=probability,
            stations=stations,
            settled=settled,
            terminal_writes=int(np.count_nonzero(commit)),
            terminal_reachable=int((df >= 0).sum() - len(stations)),
            epochs=records,
            standard=[estimate(i) for i in range(3)],
        ),
        commit,
        sign,
    )


def run():
    clock = clocks()
    panels = []
    arrays = {}
    for level in range(3):
        for j in range(NSEEDS):
            seed = 280210 + 10 * level + j
            primary, commits, sign = trial(level, seed, 0.25, clock)
            quick, other, w = trial(level, seed, 1.0, clock)
            assert np.array_equal(sign, w)
            arrays[f"epoch_{seed}"] = commits
            arrays[f"sign_{seed}"] = sign
            arrays[f"quick_{seed}"] = other
            panels.extend((primary, quick))
    return dict(calibration=calibrate(), clock=clock, panels=panels), arrays


def main():
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("--write", action="store_true")
    report, arrays = run()
    text = json.dumps(report, indent=1) + "\n"
    if args.parse_args().write:
        (HERE / "results.json").write_text(text)
        np.savez_compressed(HERE / "records.npz", **arrays)
    else:
        if json.loads((HERE / "results.json").read_text()) != json.loads(text):
            raise ValueError("Range or clock report changed")
        with np.load(HERE / "records.npz") as saved:
            assert set(arrays) == set(saved)
            for k, v in arrays.items():
                np.testing.assert_array_equal(v, saved[k])
        print(
            "Verified 48 local repair histories, ranging surveys and independent dimer timing."
        )


if __name__ == "__main__":
    main()
