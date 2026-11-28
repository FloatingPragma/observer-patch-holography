"""Repeat the registered remote read with independent shot banks for both charges."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import log, sqrt

import numpy as np
from scipy import sparse

from . import relay

HERE = relay.HERE
RECEIPT = HERE / "relay_precision.json"
SHOTS, BATCH, SEED = 2**74, 2**62, 191131


def exact_reader(board, path, radii, charge):
    h = relay.generator(board, [(path[:2], 1)], charge)
    doubled = 2 * (h - h[0, 0] * sparse.eye(h.shape[0], format="csr"))
    integer_data = np.rint(doubled.data)
    np.testing.assert_allclose(doubled.data, integer_data, atol=2e-12, rtol=0)
    j = sparse.csr_matrix(
        (integer_data.astype(np.int64), doubled.indices, doubled.indptr), shape=h.shape
    )
    r = np.repeat(radii.astype(np.int64), 2)[:, None]
    pair = []
    for rail in (path[:2], path[1:]):
        s = np.zeros((2 * board.n, 2), np.int64)
        for port in rail:
            s[2 * port : 2 * port + 2] = np.eye(2, dtype=np.int64)
        v = j @ s
        corner = s.T @ (r * (j @ v))
        pair.append(2 * v.T @ (r * v) - corner - corner.T)
    # Two factors of J=2(H-cI) and the unnormalized two-port launch.
    return (pair[1] - pair[0]).tolist()  # numerator over 8


def draws(law, rng):
    law = np.array(law, dtype=float)
    if law.min() < -1e-13 or abs(law.sum() - 1) > 5e-12:
        raise ValueError("Invalid terminal law")
    law = np.maximum(law, 0)
    law /= law.sum()
    # Put the largest outcome last so tiny tails are drawn explicitly.
    order = np.argsort(law, kind="stable")
    batches = rng.multinomial(BATCH, law[order], size=SHOTS // BATCH)
    sums = batches.sum(axis=0, dtype=object)
    return [int(x) for x in sums[np.argsort(order)]]


def intervals(values, radii, steps, prediction):
    def contrast(lo, hi, scale):
        center = float((values[hi] - values[lo]) / scale)
        half = (radii[hi] + radii[lo]) / float(scale)
        return [center - half, center + half]

    ordinary, held = [
        contrast(i, i + 1, 2 * Fraction(steps[name]))
        for i, name in ((0, "ordinary"), (2, "held"))
    ]
    pred = [x * y for x in prediction for y in ordinary]
    ratio = None
    if ordinary[0] > 0 or ordinary[1] < 0:
        ratios = [v / u for u in ordinary for v in held]
        ratio = [min(ratios), max(ratios)]
    return dict(
        ordinary=ordinary,
        held=held,
        held_over_ordinary=ratio,
        compatible=None
        if ratio is None
        else (held[1] >= min(pred) and held[0] <= max(pred)),
        formed=[contrast(i, i + 1, 1) for i in (4, 6)],
    )


def radius(count, factor):
    # Binomial Bernstein tails for both outcomes. Inverting the lower tail
    # gives mu-k <= L+sqrt(L^2+2kL) <= 2L+sqrt(2kL).
    n = min(count, SHOTS - count)
    return (sqrt(2 * n * factor) + 2 * factor) / SHOTS


def shell_read(hist, factor):
    if any(type(n) is not int or n < 0 for n in hist) or sum(hist) != SHOTS:
        raise ValueError("Incomplete shell counts")
    suffix, numerator, half = SHOTS, 0, 0.0
    for n in hist[:-1]:
        suffix -= n
        numerator += suffix
        half += radius(suffix, factor)
    return Fraction(numerator, SHOTS), half


def coefficients(numerators):
    (a, b), (c, d) = numerators
    if b != c:
        raise ValueError("Expected a real symmetric release operator")
    return Fraction(a + d, 16), Fraction(b, 8), Fraction(a - d, 16)


def metadata():
    def digest(file):
        return sha256(file.read_bytes()).hexdigest()

    primary = json.loads(relay.RECEIPT.read_text())
    interval = primary["source"]["predicted_ratio_interval"]
    return dict(
        protocol_sha256=digest(HERE / "relay_precision.md"),
        primary_sha256=digest(relay.RECEIPT),
        seed=SEED,
        shots=SHOTS,
        stocks=list(relay.NAMES),
        path=primary["stations"][0]["path"],
        half_steps=primary["half_steps"],
        source_interval=[float(x) for x in interval.strip("[]").split(",")],
    )


def decode(data):
    expected = metadata()
    if any(data[k] != value for k, value in expected.items()):
        raise ValueError("Acquisition differs from registered source or protocol")
    primary = json.loads(relay.RECEIPT.read_text())
    for i, station in enumerate(primary["stations"]):
        for j, name in enumerate(("-1", "1")):
            old = station["appointments"]["8"]["probes"][name]
            np.testing.assert_allclose(
                data["operators_times_8"][i][j],
                8 * np.asarray(old["acceleration_operator"]),
                atol=1e-8,
                rtol=0,
            )
            if len(old["counts"][0][0][0]) != data["shells"]:
                raise ValueError("Missing terminal shells")
    shell = data["shell"]
    shape = np.asarray(shell, dtype=object).shape
    if len(shape) != 5 or shape[:4] != (2, 8, 2, 3) or shape[4] != data["shells"]:
        raise ValueError("Unexpected shell acquisition shape")
    if np.shape(data["pauli_plus"]) != (8, 3):
        raise ValueError("Unexpected tomography bank")
    nsettings = 24 + 2 * 8 * 2 * 3 * (shape[4] - 1)
    factor = log(4 * nsettings / 0.01)
    axes, axis_errors = [], []
    for plus in data["pauli_plus"]:
        if any(type(n) is not int or not 0 <= n <= SHOTS for n in plus):
            raise ValueError("Invalid Pauli counts")
        axes.append([Fraction(2 * n - SHOTS, SHOTS) for n in plus])
        axis_errors.append([2 * radius(n, factor) for n in plus])
    signs = [coefficients(m) for m in data["operators_times_8"][0]]
    by_charge = []
    for settings in shell:
        finite = []
        for stock in settings:
            value, error = Fraction(0), 0.0
            for side, rail in zip((-1, 1), stock, strict=True):
                for weight, hist in zip((1, -2, 1), rail, strict=True):
                    m, e = shell_read(hist, factor)
                    value += side * weight * m
                    error += abs(weight) * e
            finite.append((value / Fraction(relay.DT) ** 2, error / relay.DT**2))
        by_charge.append(finite)
    out = {}
    for name, weights in (
        ("vacancy", (1, 0)),
        ("excess", (0, 1)),
        ("even", (Fraction(1, 2), Fraction(1, 2))),
        ("odd", (Fraction(1, 2), Fraction(-1, 2))),
    ):
        c, x, z = (
            sum(w * sign[i] for w, sign in zip(weights, signs)) for i in range(3)
        )
        release = [
            (c + x * a[0] + z * a[2], abs(x) * e[0] + abs(z) * e[2])
            for a, e in zip(axes, axis_errors, strict=True)
        ]
        finite = [
            (
                sum(w * a[0] for w, a in zip(weights, arms)),
                sum(abs(w) * a[1] for w, a in zip(weights, arms)),
            )
            for arms in zip(*by_charge, strict=True)
        ]
        out[name] = {
            method: intervals(*zip(*pairs), data["half_steps"], data["source_interval"])
            for method, pairs in (("release", release), ("finite", finite))
        }
    return out


def allowable_weights(panels, source):
    """Necessary weights for a fixed mixture to carry the source at each time."""
    allowed = {}
    for method in ("release", "finite"):
        residuals = []
        for name in ("vacancy", "excess"):
            row = panels[name][method]
            products = [r * x for r in source for x in row["ordinary"]]
            residuals.append(
                (row["held"][0] - max(products), row["held"][1] - min(products))
            )
        # The residual interval at vacancy weight w is w*vac+(1-w)*exc.
        # Include every w where the lower envelope <=0 and the upper >=0.
        lo, hi = 0.0, 1.0
        for k, sign in ((0, 1), (1, -1)):
            end, start = (sign * pair[k] for pair in residuals)
            if min(start, end) > 0:
                hi = -1
            elif max(start, end) > 0:
                boundary = -start / (end - start)
                if end > start:
                    hi = min(hi, boundary)
                else:
                    lo = max(lo, boundary)
        allowed[method] = [lo, hi] if lo <= hi else None
    return allowed


def acquire():
    data = metadata()
    rng = np.random.default_rng(SEED)
    board = relay.FORMATION.federation(1)
    route = relay.FORMATION.region(board)
    radii, _, _ = relay.tidal.chart(board, route)
    paths = relay.receiver_paths(board, route, radii)
    np.testing.assert_equal(paths[0], data["path"])
    stocks = relay.sources(board, route)[0]
    states = relay.deliver(board, route, paths[0][:2], stocks)[-1]
    axes = np.einsum("kij,sji->sk", relay.FORMATION.PAULI, states).real
    data["pauli_plus"] = [
        [draws([(1 + a) / 2, (1 - a) / 2], rng)[0] for a in row] for row in axes
    ]
    data["operators_times_8"] = [
        [exact_reader(board, p, radii, q) for q in (-1, 1)] for p in paths
    ]
    data["shells"] = int(radii.max()) + 1
    data["shell"] = []
    for q in (-1, 1):
        accel, effects = relay.position_reader(board, paths[0], radii, q)
        np.testing.assert_allclose(
            8 * accel, data["operators_times_8"][0][int(q > 0)], atol=1e-9, rtol=0
        )
        laws = np.einsum("tmrij,sji->srtm", effects, states).real
        data["shell"].append([[[draws(t, rng) for t in r] for r in s] for s in laws])
    data["results"] = decode(data)
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--counts-only", action="store_true")
    args = parser.parse_args()
    if args.counts_only:
        data = json.loads(RECEIPT.read_text())
        relay.tidal.verify(decode(data), data["results"])
    else:
        data = acquire()
        if args.write:
            with RECEIPT.open("x") as f:
                json.dump(data, f, separators=(",", ":"), allow_nan=False)
                f.write("\n")
        else:
            relay.tidal.verify(data, json.loads(RECEIPT.read_text()))
    print(
        "Possible vacancy weights",
        allowable_weights(data["results"], data["source_interval"]),
    )
    for name, methods in data["results"].items():
        for method, panel in methods.items():
            print(name, method, panel["held_over_ordinary"], panel["formed"])
