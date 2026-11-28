"""Infer free-medium exterior responses from separately recorded source settings."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy import sparse
from scipy.linalg import eigh

from . import undisturbed_medium as medium

lab = medium.lab
HERE = medium.HERE
AGES = (1 / 32, 4, 64, 1024, 65536)
READS = 2**46


def ordinary_gain(tape):
    clocks = tape["clock"]
    return (
        lab.delta(
            lab.rate_interval(clocks["far_standard"]),
            lab.rate_interval(clocks["near_standard"]),
        )
        - 0.5
    )


def absolute_budget(board, marks, tape):
    """Simultaneous closed-neighborhood indicator bounds, one sample per arm."""
    gain = ordinary_gain(tape)
    fields = []
    degrees = []
    for name, mask in (("blank", []), ("written", marks)):
        counts = np.array(tape["source_counts"][name])
        assert len(counts) == board.n and min(counts) >= 0 and sum(counts) == 2**40
        p = counts / sum(counts)
        adjacency = board.adjacency().tolil()
        for u, v in mask:
            adjacency[u, v] = adjacency[v, u] = 0
        adjacency = adjacency.tocsr()
        mean = p + adjacency @ p
        degree = np.asarray(adjacency.sum(axis=0)).ravel()
        degrees.append(float(degree @ p))
        np.testing.assert_allclose(mean.sum(), 1 + degrees[-1])
        e = np.sqrt(lab.RISK / (2 * sum(counts)))
        fields.append(np.stack([np.maximum(0, mean - e), np.minimum(1, mean + e)]))
    # The same ordinary gain multiplies the difference on both boards.
    lower, upper = fields[1] - fields[0][::-1]
    positive = gain * np.array([lower.clip(min=0).sum(), upper.clip(min=0).sum()])
    negative = gain * np.array([(-upper).clip(min=0).sum(), (-lower).clip(min=0).sum()])
    return dict(
        positive=positive.tolist(),
        negative=negative.tolist(),
        mean_live_degrees=degrees,
        maximum_degree=max(map(len, board.neighbours)),
        common_gain=gain.tolist(),
        native_ceiling=4 * lab.K * 7,
    )


def instruments(board, shells, mask):
    hb = lab.cold.motion(board, mask, 1, lab.K)[0]
    h = hb.toarray()
    instruments = []
    for u, v in medium.cuts(board, shells):
        d = np.bincount(u, minlength=board.n) - np.bincount(v, minlength=board.n)
        w = board.adjacency() @ d + d
        b = np.zeros_like(h)
        b[u, v], b[v, u] = 0.5, -0.5
        correction = h @ b - b @ h
        observable = 4 * lab.K * np.diag(w) + correction
        instruments.append((w, correction, observable))
    return hb, instruments


def verify_release(board, mask, shells, meters):
    h, i, j = lab.ordinary.two_loads(board, lab.K, mask)
    h2 = (2 * h).astype(np.int64)
    assert (h2 != 2 * h).nnz == 0
    w = sparse.csr_matrix(
        (np.ones(len(j), dtype=np.int64), (np.arange(len(j)), j)),
        shape=(len(j), board.n),
    )
    a = h2 @ w
    b = h2 @ a
    for radius, (_, _, observable) in enumerate(meters, start=1):
        inside = shells[i] <= radius
        past = w[inside].T @ b[inside]
        numerator = (2 * (a[inside].T @ a[inside]) - past - past.T).toarray()
        np.testing.assert_array_equal(numerator, 4 * observable)


def averaged(matrix, eig, age, *, heisenberg=False):
    e, v = eig
    gap = (e[:, None] - e[None, :]) * age / 2
    phase = np.exp((1j if heisenberg else -1j) * gap) * np.sinc(gap / np.pi)
    return v @ ((v.T @ matrix @ v) * phase) @ v.T


def preparations(board, response):
    n = board.n
    pairs = []
    candidates = [
        (a, b)
        for a in range(n)
        for b in range(a + 1, n)
        if b not in board.neighbours[a]
    ]
    candidates.sort(key=lambda pair: (-round(abs(response[pair].real), 10), pair))
    available = set(range(n))
    for a, b in candidates:
        if a in available and b in available:
            pairs.append((a, b, 1 if response[a, b].real >= 0 else -1))
            available.difference_update((a, b))
    x = np.zeros((n, n))
    for a, b, sign in pairs:
        x[a, b] = x[b, a] = sign
    eye = np.eye(n)
    return dict(
        random_port=eye / n,
        rest=np.ones((n, n)) / n,
        paired_plus=(eye + x) / n,
        paired_minus=(eye - x) / n,
    ), pairs


def fringe_counts(c, rho, rng):
    upper = np.triu(abs(c), 1)
    norm = np.trace(abs(c)) + upper.sum()
    p = np.diag(rho).real
    active = (abs(np.diag(c)) @ p + np.sum(upper * (p[:, None] + p[None, :]))) / norm
    contrast = np.sum(c * rho.real) / norm
    probabilities = np.maximum(
        0, [(active + contrast) / 2, (active - contrast) / 2, 1 - active]
    )
    return rng.multinomial(READS, probabilities / probabilities.sum()).tolist()


def decode(data, tape):
    board = lab.cold.federation(0)
    marks, _ = lab.cold.read_map(board, tape["spectral_tape"])
    shells = lab.cold.distances(board, marks[0])
    gain = ordinary_gain(tape)
    meters = instruments(board, shells, [])[1]
    out = {}
    for name, rows in data["observations"].items():
        out[name] = {}
        for age, arms in rows.items():
            bands = []
            for index, (w, c, _) in enumerate(meters):
                norm = np.trace(abs(c)) + np.triu(abs(c), 1).sum()
                field, correction = [], []
                for tag in ("blank", "written"):
                    bank = arms[tag]
                    count = bank["sites"]
                    fringe = bank["fringes"][index]
                    assert (
                        len(count) == board.n
                        and min(count) >= 0
                        and sum(count) == READS
                    )
                    assert (
                        len(fringe) == 3 and min(fringe) >= 0 and sum(fringe) == READS
                    )
                    width = np.sqrt(lab.RISK / (2 * READS))
                    signs = np.array([-1, 1])
                    field.append(w @ count / READS + signs * np.ptp(w) * width)
                    correction.append(
                        norm * ((fringe[0] - fringe[1]) / READS + signs * 2 * width)
                    )
                f = lab.signed_product(gain, lab.delta(field[1], field[0]))
                z = lab.delta(correction[1], correction[0])
                bands.append(((f + z) / (board.n - 1)).tolist())
            out[name][age] = bands
    return out


def run(tape):
    board = lab.cold.federation(0)
    marks, _ = lab.cold.read_map(board, tape["spectral_tape"])
    shells = lab.cold.distances(board, marks[0])
    banks = {}
    eig = {}
    for name, mask in (("blank", []), ("written", marks)):
        h, meters = instruments(board, shells, mask)
        verify_release(board, mask, shells, meters)
        banks[name] = meters
        eig[name] = eigh(h.toarray())
    for a, b in zip(banks["blank"], banks["written"], strict=True):
        np.testing.assert_array_equal(a[2], b[2])

    def first(name):
        return averaged(banks[name][0][2], eig[name], 64, heisenberg=True)

    preps, pairs = preparations(board, first("written") - first("blank"))
    rng = np.random.default_rng(376489)
    observations, predicted = {}, {}
    for prep, rho in preps.items():
        rows, theory = {}, {}
        for age in AGES:
            arm, values = {}, []
            for tag in ("blank", "written"):
                evolved = averaged(rho, eig[tag], age)
                p = np.maximum(np.diag(evolved).real, 0)
                arm[tag] = dict(
                    sites=rng.multinomial(READS, p / p.sum()).tolist(),
                    fringes=[fringe_counts(c, evolved, rng) for _, c, _ in banks[tag]],
                )
                values.append(
                    np.array([np.sum(o * evolved.real) for _, _, o in banks[tag]])
                )
            rows[str(age)] = arm
            theory[str(age)] = ((values[1] - values[0]) / (board.n - 1)).tolist()
        observations[prep], predicted[prep] = rows, theory
    data = dict(
        protocol_sha256=sha256((HERE / "ISOLATED_MEDIUM.md").read_bytes()).hexdigest(),
        reference_sha256=sha256(
            (HERE / "mediated_release.json").read_bytes()
        ).hexdigest(),
        seed=376489,
        trials=READS,
        ports=board.n,
        ages=list(AGES),
        selected_pairs=pairs,
        stationary_budget=absolute_budget(board, marks, tape),
        observations=observations,
    )
    data["decoded"] = decode(data, tape)
    for prep, arm in predicted.items():
        for age, row in arm.items():
            lo, hi = np.asarray(data["decoded"][prep][age]).T
            assert np.all((row >= lo) & (row <= hi))
    for age in predicted["paired_plus"]:
        np.testing.assert_allclose(
            np.array(predicted["paired_plus"][age]) + predicted["paired_minus"][age],
            0,
            atol=1e-12,
        )
        np.testing.assert_allclose(predicted["random_port"][age], 0, atol=1e-12)
    return data


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    flags = cli.add_mutually_exclusive_group()
    flags.add_argument("--write", action="store_true")
    flags.add_argument("--counts-only", action="store_true")
    args = cli.parse_args()
    tape = json.loads((HERE / "mediated_release.json").read_bytes())
    path = HERE / "isolated_medium.json"
    if args.counts_only:
        data = json.loads(path.read_bytes())
        lab.cold.apparatus.check(data["decoded"], decode(data, tape))
        board = lab.cold.federation(0)
        marks, _ = lab.cold.read_map(board, tape["spectral_tape"])
        lab.cold.apparatus.check(
            data["stationary_budget"], absolute_budget(board, marks, tape)
        )
    else:
        data = run(tape)
        if args.write:
            path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
        else:
            lab.cold.apparatus.check(
                json.loads(path.read_bytes()), json.loads(json.dumps(data))
            )
    print(data["stationary_budget"])
    for prep, row in data["decoded"].items():
        print(prep, {age: band[0] for age, band in row.items()})
