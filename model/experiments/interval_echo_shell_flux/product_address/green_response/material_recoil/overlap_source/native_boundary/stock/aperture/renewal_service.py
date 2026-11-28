"""Unconditional native service, event sampling, and the clock's decorrelation time."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from . import renewal_clock as clock

HERE = clock.HERE
lab = clock.lab


def service(full, stock):
    W, V, *_ = stock.kernels()
    edges = stock.carrier()[3]
    ix = {tuple(w): i for i, w in enumerate(full)}
    rows, cols, values = [], [], []
    for i, word in enumerate(full):
        for j, (a, b) in enumerate(edges):
            code = 9 * word[4 + j] + 3 * word[a] + word[b]
            mix = W[:, code] + V[:, code]
            for k in np.flatnonzero(mix):
                out = word.copy()
                out[4 + j], ports = divmod(k, 9)
                out[a], out[b] = divmod(ports, 3)
                rows.append(i)
                cols.append(ix[tuple(out)])
                values.append(int(mix[k]))
    tickets = sparse.csr_matrix(
        (values, (rows, cols)), shape=(len(full), len(full)), dtype=int
    )
    assert (tickets != tickets.T).nnz == 0
    assert np.all(tickets.sum(1) == 384) and min(tickets.diagonal()) > 0
    assert connected_components(tickets, directed=False)[0] == 1
    return tickets / 384


def dynamics(kernel, flag, change, report):
    uniform = np.ones(len(flag)) / len(flag)
    advance = flag if report == "advance_at_tick" else np.ones(len(flag))
    move = sparse.eye(len(flag), format="csr") + sparse.diags(advance) @ (
        kernel - sparse.eye(len(flag), format="csr")
    )
    pi = uniform / advance
    pi /= pi.sum()
    np.testing.assert_allclose(pi @ move, pi, atol=3e-17)
    tick = pi * flag
    rate = tick.sum()
    tick /= rate
    np.testing.assert_allclose(tick if report == "advance_at_tick" else pi, uniform)
    # The source records seen by an integrating observer determine the event
    # records by rate weighting, including when word and residence are correlated.
    contrast = change.copy()
    diff = {}
    for i in range(4097):
        if i in (0, 16, 64, 256, 1024, 4096):
            diff[str(i)] = dict(
                rate_difference=float(contrast @ flag),
                total_variation=float(abs(contrast).sum() / 2),
            )
        contrast = move.T @ contrast
    centered = flag - rate
    var = float(pi @ (centered**2))
    weight = pi * centered / var
    # The reversible kernel gives the covariance at all lags from one
    # symmetric eigensystem. Evaluate the exact finite triangular sum of
    # each spectral mode at powers of two, including long exposures.
    root = np.sqrt(pi)
    symmetric = move.toarray() * root[:, None] / root
    np.testing.assert_allclose(symmetric, symmetric.T, atol=1e-15)
    eigen, vectors = np.linalg.eigh(symmetric)
    assert abs(eigen[-1] - 1) < 1e-12 and eigen[-2] < 1 - 1e-10
    strengths = (vectors.T @ (root * centered)) ** 2 / var
    assert strengths[-1] < 1e-19
    mode = eigen[:-1].astype(np.longdouble)
    weights = strengths[:-1].astype(np.longdouble)
    counts = {}
    checkpoints = {}
    for power in range(28):
        n = 2**power
        terms = n * (1 + mode) / (1 - mode) - 2 * mode * (1 - mode**n) / (1 - mode) ** 2
        eff = float(n * n / (weights @ terms))
        if n in (1, 16, 64, 256, 1024, 4096):
            checkpoints[str(n)] = eff
        for threshold in (11, 16, 21, 35, 54, 57, 108, 216):
            if eff >= threshold and str(threshold) not in counts:
                counts[str(threshold)] = n
    assert len(counts) == 8
    cov = 0.0
    subtotal = 0.0
    prop = centered.copy()
    for k in range(1, 257):
        subtotal += 1 + 2 * cov
        if str(k) in checkpoints:
            np.testing.assert_allclose(k * k / subtotal, checkpoints[str(k)], rtol=2e-9)
        prop = move @ prop
        cov += weight @ prop
    return (
        pi,
        tick,
        dict(
            stationary_tick_probability=float(rate),
            initial_contrast=diff,
            effective_samples_at_exposure=checkpoints,
            first_dyadic_slots_for_effective_samples=counts,
        ),
    )


def run():
    full, w, cells, *_ = lab.energy.pair.setup()
    stock = lab.energy.pair.stock
    P = service(full, stock)
    lookup = {tuple(word): i for i, word in enumerate(full)}
    ix = [lookup[tuple(word)] for word in w]
    delta = np.zeros(len(full))
    delta[ix] = np.asarray(lab.energy.least_change(cells) - 1, float) / len(w)
    mask = np.zeros(len(full), bool)
    mask[ix] = True
    bits = stock.MAG[full[:, 4:6]] == 2
    code = 2 * bits[:, 0] + bits[:, 1].astype(int)

    def table(p):
        counts = np.bincount(code[mask], weights=p[mask], minlength=4)
        return dict(
            accepted=float(counts.sum()),
            conditional=(counts / counts.sum()).tolist(),
            odds=float(counts[0] * counts[3] / (counts[1] * counts[2])),
        )

    out = {}
    for load in (1, 2):
        _, bank = clock.read_clock(stock.kernels()[2], [12, 17, 19, 20, 21], load)
        flag = bank[full[:, 1], 0]
        row = {}
        for name in ("advance_each_slot", "advance_at_tick"):
            pi, event, report = dynamics(P, flag, delta, name)
            report.update(random_slot_heads=table(pi), at_tick_heads=table(event))
            row[name] = report
        out[str(load)] = row
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "renewal_clock.json", HERE / "CLOCK_RELEASE.md")
        },
        words=len(full),
        service_tickets=384,
        connected_classes=1,
        original_head_odds=str(F(3, 11)),
        readers=out,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "renewal_service.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        lab.check(json.loads(dst.read_bytes()), result)
    print(
        result["words"],
        {
            t: {
                k: v["first_dyadic_slots_for_effective_samples"]["108"]
                for k, v in row.items()
            }
            for t, row in result["readers"].items()
        },
    )
