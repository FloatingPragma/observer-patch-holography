"""Coincidences between the distributing courier and two remote quantum records."""

import argparse
from hashlib import sha256
import json
from math import log, sqrt

import numpy as np

from . import paired_relay as pair, relay

HERE = relay.HERE
RECEIPT = HERE / "paired_arrivals.json"
SEED = 191205
SHOTS = 2**74
FLAGS = ("first", "second", "neither")
# Include the six cells, flag-only and flag-forgotten marginals simultaneously.
L = log(4 * 2 * 3 * 15 * (6 + 3 + 2) / 0.01)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def cohorted(board, route, paths, q, stocks):
    owners = np.arange(board.n) // 12
    first, second = (x[0] // 12 for x in paths)
    masks = (owners == first, owners == second, (owners != first) & (owners != second))
    states = np.zeros((len(stocks), 3, 4, 4), complex)
    for block, w in enumerate(pair.flights(board, route, paths, q)):
        for i, mask in enumerate(masks):
            v = w[mask].reshape(-1, 4, 4)
            c = np.einsum("pca,pdb->cdab", v, v.conj(), optimize=True)
            states[:, i] += np.einsum("cdab,sab->scd", c, stocks[:, block])
    np.testing.assert_allclose(
        states.trace(axis1=-2, axis2=-1).sum(axis=1), 1, atol=2e-12, rtol=0
    )
    if np.linalg.eigvalsh(states).min() < -3e-13:
        raise ValueError("Nonpositive reception block")
    return states


def banks(states, rng):
    rows = []
    for state in states:
        pauli = np.einsum("aij,fji->af", pair.OPS, state).real
        flags = state.trace(axis1=-2, axis2=-1).real
        outcomes = np.stack((flags + pauli, flags - pauli), axis=-1) / 2
        rows.append(
            [
                np.reshape(pair.draws(x.ravel(), SHOTS, rng), (3, 2)).tolist()
                for x in outcomes
            ]
        )
    # Numpy infers object from these counts, but keep JSON's integer type
    # explicit regardless of the sizes of individual bins.
    return [[[[int(x) for x in f] for f in axis] for axis in stock] for stock in rows]


def error(n):
    return (sqrt(2 * min(n, SHOTS - n) * L) + 2 * L) / SHOTS


def contrast(arms):
    first, second = arms
    center = (first - second) / SHOTS
    half = error(first) + error(second)
    return [center - half, center + half]


def decode(data):
    if (
        data["protocol_sha256"] != digest(HERE / "paired_arrivals_protocol.md")
        or data["source_sha256"] != digest(pair.RECEIPT)
        or data["seed"] != SEED
        or data["shots"] != SHOTS
    ):
        raise ValueError("Changed acquisition")
    out = {}
    for q in ("-1", "1"):
        a = np.array(data["banks"][q], dtype=object)
        if (
            a.shape != (3, 15, 3, 2)
            or any(type(x) is not int or x < 0 for x in a.flat)
            or np.any(a.sum(axis=(2, 3)) != SHOTS)
        ):
            raise ValueError("Incomplete joint count bank")
        comparisons = {}
        for idx, tag in ((1, "source"), (2, "control")):
            arms = a[[0, idx]]
            events = []
            for axis in range(15):
                for f in range(3):
                    for sign in (0, 1):
                        v = contrast(arms[:, axis, f, sign])
                        events.append((max(0, v[0], -v[1]), axis, f, sign, v))
            lower, axis, f, sign, v = max(events)
            # Flag probability is read on the first preassigned Pauli bank.
            flags = [contrast(arms[:, 0, i].sum(axis=1)) for i in range(3)]
            marginal = [contrast(arms[:, j, :, 0].sum(axis=1)) for j in range(15)]
            comparisons[tag] = dict(
                event_lower_bound=lower,
                axis=pair.LABELS[axis],
                flag=FLAGS[f],
                outcome=sign,
                interval=v,
                flags=flags,
                forgotten_largest_lower=max(max(0, a, -b) for a, b in marginal),
            )
        out[q] = comparisons
    return out


def acquire():
    prior = json.loads(pair.RECEIPT.read_text())
    basis = pair.meter.head.joint_basis()
    p = np.array([3 * np.eye(4) / 16, np.eye(4) / 16])
    r = p + pair.H * np.einsum("i,ibjk->bjk", prior["source"]["direction"], basis)
    board = relay.FORMATION.federation(1)
    route = relay.FORMATION.region(board)
    raw = {}
    distance = {}
    rng = np.random.default_rng(SEED)
    for q in (-1, 1):
        states = cohorted(board, route, prior["paths"], q, np.array([p, r, p]))
        diff = states[0] - states[1]
        distance[str(q)] = dict(
            complete=float(abs(np.linalg.eigvalsh(diff)).sum() / 2),
            forgotten=float(abs(np.linalg.eigvalsh(diff.sum(axis=0))).sum() / 2),
        )
        raw[str(q)] = banks(states, rng)
    data = dict(
        protocol_sha256=digest(HERE / "paired_arrivals_protocol.md"),
        source_sha256=digest(pair.RECEIPT),
        seed=SEED,
        shots=SHOTS,
        banks=raw,
        distances=distance,
    )
    data["measurements"] = decode(data)
    return data


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    m = p.add_mutually_exclusive_group()
    m.add_argument("--counts-only", action="store_true")
    m.add_argument("--write", action="store_true")
    a = p.parse_args()
    if a.counts_only:
        d = json.loads(RECEIPT.read_text())
        relay.tidal.verify(decode(d), d["measurements"])
    else:
        d = acquire()
        if a.write:
            with RECEIPT.open("x") as f:
                json.dump(d, f, separators=(",", ":"), allow_nan=False)
                f.write("\n")
        else:
            reference = json.loads(RECEIPT.read_text())
            relay.tidal.verify(decode(reference), reference["measurements"])
            np.testing.assert_allclose(
                [p["complete"] for p in d["distances"].values()],
                [p["complete"] for p in reference["distances"].values()],
                rtol=1e-8,
                atol=1e-13,
            )
    print(d["distances"])
    for q, obs in d["measurements"].items():
        print(q, obs["source"], obs["control"]["event_lower_bound"])
