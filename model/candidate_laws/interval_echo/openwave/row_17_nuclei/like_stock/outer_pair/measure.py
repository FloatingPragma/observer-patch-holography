"""Clock-coincident graph distances and monitored first passages of the signed pair."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import shortest_path
from scipy.stats import beta

from ..measure import CONTRASTS, EDGES, SHOTS, shore_reads
from ..two import STEP
from ....row_16_binding.signed_pair.material import Assembly

HERE = Path(__file__).resolve().parent
SEED = 26112733
WAITS = (0, 1, 4, 16, 64)


def distances():
    adjacency = np.zeros((12, 12), dtype=int)
    for a, b in EDGES:
        adjacency[a, b] = adjacency[b, a] = 1
    dist = shortest_path(adjacency, unweighted=True, directed=False).astype(int)
    np.testing.assert_array_equal(np.bincount(dist.ravel()), [12, 60, 60, 12])
    return dist


def bright(assembly, ref, delays):
    vals, u = eigh(assembly.h.toarray())
    z = (vals - ref) * STEP / (2 * np.pi)
    weights = (np.sinc(delays * z) / np.sinc(z)) ** 2 / len(vals)
    rho = (u * weights) @ u.T
    return rho, vals, u


def rates(rho, shell):
    return np.bincount(shell, weights=rho.diagonal().real, minlength=4)[1:]


def propagation(values, u, t):
    if t == 0:
        return np.eye(len(values))
    return (u * np.exp(-1j * t * values)) @ u.T


def sample(rng, p):
    values = np.asarray(p, dtype=float)
    if values.min() < -2e-12 or abs(values.sum() - 1) > 1e-11:
        raise ValueError("unphysical cohort")
    values = np.maximum(values, 0)
    order = np.argsort(values, kind="stable")
    return rng.multinomial(SHOTS, values[order] / values.sum())[
        np.argsort(order)
    ].tolist()


def acquire():
    parent = json.loads((HERE.parent / "two_counts.json").read_text())
    assert parent["step"] == STEP and parent["delays"] == 131072
    assert shore_reads() == [-1, 0, 1]
    rng = np.random.default_rng(SEED)
    dist = distances()
    records = []
    for k, row in zip(CONTRASTS, parent["rows"], strict=True):
        assert k == row["k"]
        pair = Assembly(2, 0, k)
        single = Assembly(1, 0, k)
        supports = np.array([np.flatnonzero(w == 1) for w in pair.w])
        shell = dist[supports[:, 0], supports[:, 1]]
        rho, values, u = bright(pair, row["refs"][1], parent["delays"])
        prob = rho.trace()
        single_rho, _, _ = bright(single, row["refs"][0], parent["delays"])
        # Single-copy sharp support replaces an arbitrary simulator basis ordering.
        pd = np.zeros(12)
        pd[np.argmax(single.w, axis=1)] = single_rho.diagonal()
        copy_rates = np.bincount(
            dist.ravel(), weights=np.outer(pd, pd).ravel(), minlength=4
        )
        controls = dict(
            independent=sample(rng, [*copy_rates, 1 - copy_rates.sum()]),
            unfiltered=sample(rng, rates(np.eye(len(rho)) / len(rho), shell)),
        )
        terminal = []
        for t in WAITS:
            evo = propagation(values, u, t)
            groups = []
            for sharp in (False, True):
                joint = []
                for initial in (shell == 1, shell > 1):
                    state = rho * np.outer(initial, initial)
                    if sharp:
                        state = np.diag(np.diag(state))
                    joint.extend(rates(evo @ state @ evo.conj().T, shell))
                groups.append(sample(rng, [*joint, 1 - prob]))
            terminal.append(
                dict(
                    t=t,
                    unread=sample(
                        rng, [*rates(evo @ rho @ evo.conj().T, shell), 1 - prob]
                    ),
                    coarse=groups[0],
                    sharp=groups[1],
                )
            )
        passages = []
        for dt in (1, 4):
            evo = propagation(values, u, dt)
            near = shell == 1
            state = rho * np.outer(near, near)
            initial_near = state.trace()
            first = []
            for _ in range(16 // dt):
                evolved = evo @ state @ evo.conj().T
                first.append(evolved.diagonal().real[~near].sum())
                state = evolved * np.outer(near, near)
            counts = sample(
                rng, [*first, state.trace().real, prob - initial_near, 1 - prob]
            )
            passages.append(dict(step=dt, counts=counts))
        records.append(
            dict(
                k=k,
                ground_reference=row["refs"][1],
                one_reference=row["refs"][0],
                controls=controls,
                terminal=terminal,
                passages=passages,
            )
        )
    return dict(
        seed=SEED,
        starts=SHOTS,
        delays=parent["delays"],
        step=STEP,
        stock_shells=np.bincount(shell)[1:].tolist(),
        records=records,
    )


def interval(hits, n):
    if n == 0:
        return [0.0, 1.0]
    tail = 0.01 / (2 * 4096)
    return [
        float(beta.ppf(tail, hits, n - hits + 1)) if hits else 0.0,
        float(beta.isf(tail, hits + 1, n - hits)) if hits < n else 1.0,
    ]


def shells(c):
    n = sum(c)
    far = sum(c[1:])
    p = interval(far, n)
    outer = interval(c[2], n)
    return dict(
        n=n,
        shells=c,
        far=far / n if n else None,
        far_limits=p,
        outer=c[2] / n if n else None,
        mean_square=sum((i + 1) ** 2 * v for i, v in enumerate(c)) / n if n else None,
        mean_square_limits=[1 + 3 * a + 5 * b for a, b in zip(p, outer, strict=True)],
    )


def decode(data):
    if (data["seed"], data["starts"]) != (SEED, SHOTS) or [
        r["k"] for r in data["records"]
    ] != list(CONTRASTS):
        raise ValueError("wrong enrollment")
    out = []
    for record in data["records"]:
        c = record["controls"]
        inputs = [(c["independent"], 5), (c["unfiltered"], 3)]
        if [s["t"] for s in record["terminal"]] != list(WAITS) or [
            p["step"] for p in record["passages"]
        ] != [1, 4]:
            raise ValueError("missing appointment")
        terminal = []
        for t in record["terminal"]:
            for key, length in [("unread", 4), ("coarse", 7), ("sharp", 7)]:
                inputs.append((t[key], length))
            terminal.append(
                dict(
                    t=t["t"],
                    unread=shells(t["unread"][:3]),
                    coarse_from_near=shells(t["coarse"][:3]),
                    coarse_from_far=shells(t["coarse"][3:6]),
                    sharp_from_near=shells(t["sharp"][:3]),
                    sharp_from_far=shells(t["sharp"][3:6]),
                )
            )
        passages = []
        for p in record["passages"]:
            v = p["counts"]
            inputs.append((v, 16 // p["step"] + 3))
            n = sum(v[:-2])
            hits = sum(v[:-3])
            passages.append(
                dict(step=p["step"], first=hits, near_start=n, limits=interval(hits, n))
            )
        if any(
            len(v) != m
            or sum(v) != SHOTS
            or any(type(j) is not int or j < 0 for j in v)
            for v, m in inputs
        ):
            raise ValueError("missing starts")
        out.append(
            dict(
                k=record["k"],
                geometry=shells(c["unfiltered"]),
                independent=shells(c["independent"][1:4]),
                terminal=terminal,
                passages=passages,
            )
        )
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquire", type=Path)
    args = p.parse_args()
    if args.acquire:
        args.acquire.write_text(
            json.dumps(acquire(), separators=(",", ":"), allow_nan=False) + "\n"
        )
    else:
        print(
            json.dumps(decode(json.loads((HERE / "counts.json").read_text())), indent=2)
        )
