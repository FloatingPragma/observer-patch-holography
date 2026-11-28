"""Two remote written heads as a joint receiver for a correlation source."""

import argparse
from hashlib import sha256
from itertools import product
import json

import numpy as np
from scipy.sparse.linalg import expm_multiply

from . import measure as meter, relay

HERE = relay.HERE
RECEIPT = HERE / "paired_relay.json"
H = 1 / 64
SHOTS = 2**66
SEED = 191203
AXES = (np.eye(2), *relay.FORMATION.PAULI)
LABELS = list(product(range(4), repeat=2))[1:]
OPS = np.array([np.kron(AXES[i], AXES[j]) for i, j in LABELS])


def paths(board, route, radii):
    home = set(np.array(route[1:3]) // 12)
    around = set()
    for a, b in board.external:
        i, j = a // 12, b // 12
        if i in home:
            around.add(j)
        if j in home:
            around.add(i)
    around -= home
    options = sorted(
        (radii[i], i, j, k)
        for i in range(board.n)
        if i // 12 in around
        for j in board.neighbours[i]
        if j // 12 == i // 12 and radii[j] == radii[i] + 1
        for k in board.neighbours[j]
        if k // 12 == i // 12 and radii[k] == radii[j] + 1
    )
    first = options[0][1:]
    second = next(row[1:] for row in options if row[1] // 12 != first[0] // 12)
    assert first == relay.receiver_paths(board, route, radii)[0]
    return first, second


def flights(board, route, receivers, q):
    # Four source basis columns, the remote coins initially 00. The tensor
    # product order is position, first source, second source, two receivers.
    w = np.zeros((16 * board.n, 4))
    w[16 * route[0] + 4 * np.arange(4)] = np.eye(4)
    for parity in (1, -1):
        contacts = [
            (route[:2], parity),
            (route[2:], 1),
            *((p[:2], 1) for p in receivers),
        ]
        h = relay.generator(board, contacts, q)
        out = w
        for _ in range(8):
            out = expm_multiply(-1j * h, out, traceA=0)
        np.testing.assert_allclose(out.conj().T @ out, np.eye(4), atol=2e-12, rtol=0)
        yield out.reshape(board.n, 4, 4, 4)


def channel(board, route, receivers, q):
    channels = []
    for wave in flights(board, route, receivers, q):
        v = wave.reshape(-1, 4, 4)
        channels.append(np.einsum("pca,pdb->cdab", v, v.conj(), optimize=True))
    channels = np.array(channels)
    np.testing.assert_allclose(
        channels.trace(axis1=1, axis2=2),
        np.broadcast_to(np.eye(4), (2, 4, 4)),
        atol=2e-12,
        rtol=0,
    )
    return channels


def deliver(c, stocks):
    return np.einsum("kcdab,skab->scd", c, stocks)


def response(states):
    return np.einsum("aij,sji->as", OPS, states).real


def spectrum(a):
    _, s, v = np.linalg.svd(a, full_matrices=True)
    rank = int(np.sum(s > max(1e-12, 1e-10 * s[0])))
    return s, v[rank:]


def witness(channels, basis, score):
    maps = np.array([response(deliver(c, basis)) for c in channels])
    marg = np.array([i * j == 0 for i, j in LABELS])
    scope = {}
    for key, m in (
        ("vacancy", maps[0]),
        ("excess", maps[1]),
        ("both", np.vstack(maps)),
        ("both_marginals", np.vstack(maps[:, marg])),
    ):
        s, kernel = spectrum(m)
        scope[key] = dict(
            singular_values=s.tolist(),
            numerical_rank=18 - len(kernel),
            unobserved_score_norm=float(np.linalg.norm(kernel @ score)),
        )
    _, null = spectrum(maps[1])
    v = null.T @ (null @ score)
    if np.linalg.norm(v) < 1e-14:
        v = null[0]
    v /= np.linalg.norm(v)
    direction = np.einsum("i,ibjk->bjk", v, basis)
    p = np.array([3 * np.eye(4) / 16, np.eye(4) / 16])
    r = p + H * direction
    np.testing.assert_allclose(meter.head.marginals(r - p), 0, atol=1e-15)
    lp = meter.head.station.log_blocks(p) - meter.head.station.log_blocks(r)
    coef = meter.head.station.pair(basis, lp)
    target = meter.head.station.pair(lp, p - r)
    assert target >= H**2
    return (p, r, coef), dict(
        direction=v.tolist(),
        scope=scope,
        expected_charge=float(target),
        basis_coefficients=coef.tolist(),
        source_min_eigenvalue=float(np.linalg.eigvalsh(r).min()),
        preparation=[meter.head.station.preparation(x) for x in (p, r)],
    )


def draws(p, shots, rng):
    a = np.array(p, dtype=float)
    if a.min() < -3e-12 or abs(a.sum() - 1) > 4e-12:
        raise ValueError("Invalid outcome law")
    a = np.maximum(a, 0)
    a /= a.sum()
    order = np.argsort(a, kind="stable")
    n = min(2**62, shots)
    return [
        int(x)
        for x in rng.multinomial(n, a[order], size=shots // n).sum(
            axis=0, dtype=object
        )[np.argsort(order)]
    ]


def source_counts(stocks, basis, rng):
    banks = []
    for state in stocks:
        means = meter.head.station.pair(basis, state)
        bank = []
        for i, m in enumerate(means):
            w = np.trace(state[i // 9]).real
            bank.append(draws([w / 2 + m, w / 2 - m, 1 - w], SHOTS, rng))
        banks.append(bank)
    return banks


def checked(counts, shape, shots):
    a = np.array(counts, dtype=object)
    if (
        a.shape != shape
        or any(type(x) is not int or x < 0 for x in a.flat)
        or np.any(a.sum(axis=-1) != shots)
    ):
        raise ValueError("Incomplete outcome bank")
    return a


def fixed_score_read(data):
    """Reuse the measured endpoint Pauli means for the fixed calibration score."""
    m = json.loads(meter.RECEIPT.read_text())["gates"]["0.0"]["pauli_coefficients"]
    c = np.asarray(m["original"])
    source = checked(data["source_counts"], (2, 18, 3), SHOTS)
    means = (source[:, :, 0] - source[:, :, 1]) / (2 * SHOTS)
    center = float(c @ (means[0] - means[1]))
    eps = np.sqrt(np.log(2 * (90 + 36) / 0.01) / (2 * SHOTS))
    half = float(2 * sum(abs(c)) * eps)
    return [center - half, center + half]


def decode(data):
    if (
        data["protocol_sha256"]
        != sha256((HERE / "paired_relay_protocol.md").read_bytes()).hexdigest()
        or data["shots"] != SHOTS
        or data["seed"] != SEED
        or data["stocks"] != list(relay.NAMES)
    ):
        raise ValueError("Unexpected source or registered protocol")
    plus = checked(data["readers"], (2, 3, 15, 2), SHOTS)
    source = checked(data["source_counts"], (2, 18, 3), SHOTS)
    # One simultaneous event for every mean (bounded), including controls.
    n = 2 * 3 * 15 + 2 * 18
    eps = np.sqrt(np.log(2 * n / 0.01) / (2 * SHOTS))
    c = np.array(data["source"]["basis_coefficients"])
    means = (source[:, :, 0] - source[:, :, 1]) / (2 * SHOTS)
    score = float(c @ (means[0] - means[1]))
    half = float(2 * sum(abs(c)) * eps)
    receivers = {}
    for i, name in enumerate(("vacancy", "excess")):
        b = (plus[i, :, :, 0] - plus[i, :, :, 1]) / SHOTS
        results = []
        for k in (1, 2):
            d = np.array(b[0] - b[k], dtype=float)
            # ||rho-sigma||_1/2 <= sqrt(sum_P delta_P^2)/2 for two qubits.
            upper = float(np.linalg.norm(np.abs(d) + 4 * eps) / 2)
            lower = max(0, float(np.max(np.abs(d)) - 4 * eps) / 2)
            results.append([lower, upper])
        receivers[name] = dict(source_distance=results[0], unchanged_control=results[1])
    return dict(
        modular_contrast=[score - half, score + half], receiver_distances=receivers
    )


def encode(z):
    return [z.real.tolist(), z.imag.tolist()]


def unpack(z):
    return np.array(z[0]) + 1j * np.array(z[1])


def run():
    board = relay.FORMATION.federation(1)
    route = relay.FORMATION.region(board)
    radii = relay.tidal.chart(board, route)[0]
    dest = paths(board, route, radii)
    basis = meter.head.joint_basis()
    ref, *_ = meter.head.jacobi.beam.aperture.source_and_changes()
    score = meter.head.station.log_blocks(
        meter.product_reference(ref)
    ) - meter.head.station.log_blocks(ref)
    channels = [channel(board, route, dest, q) for q in (-1, 1)]
    (p, r, _), source = witness(channels, basis, meter.head.station.pair(basis, score))
    stocks, steps, ratio, _ = relay.sources(board, route)
    rng = np.random.default_rng(SEED)
    banks = []
    distances = []
    installed = []
    for c in channels:
        recv = deliver(c, [p, r, p])
        b = response(recv).T
        banks.append(
            [[draws([(1 + x) / 2, (1 - x) / 2], SHOTS, rng) for x in v] for v in b]
        )
        distances.append(float(np.abs(np.linalg.eigvalsh(recv[0] - recv[1])).sum() / 2))
        installed.append(encode(deliver(c, stocks)))
    data = dict(
        protocol_sha256=sha256(
            (HERE / "paired_relay_protocol.md").read_bytes()
        ).hexdigest(),
        seed=SEED,
        shots=SHOTS,
        paths=[list(p) for p in dest],
        source=source,
        readers=banks,
        source_counts=source_counts((p, r), basis, rng),
        population_trace_distances=distances,
        stocks=list(relay.NAMES),
        installed=installed,
        half_steps=steps,
        source_interval=ratio,
    )
    data["measurement"] = decode(data)
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--counts-only", action="store_true")
    args = parser.parse_args()
    if args.counts_only:
        data = json.loads(RECEIPT.read_text())
        relay.tidal.verify(decode(data), data["measurement"])
    else:
        data = run()
        if args.write:
            with RECEIPT.open("x") as f:
                json.dump(data, f, separators=(",", ":"), allow_nan=False)
                f.write("\n")
        else:
            saved = json.loads(RECEIPT.read_text())
            # LAPACK is free to rotate singular directions, and extremely large
            # seeded counts can resolve last-bit differences in outcome laws.
            for k in ("paths", "stocks", "protocol_sha256"):
                assert data[k] == saved[k]
            np.testing.assert_allclose(
                data["source"]["direction"], saved["source"]["direction"], atol=1e-7
            )
            np.testing.assert_allclose(
                data["installed"], saved["installed"], atol=3e-12, rtol=0
            )
            for key, row in data["source"]["scope"].items():
                np.testing.assert_allclose(
                    row["singular_values"],
                    saved["source"]["scope"][key]["singular_values"],
                    atol=2e-12,
                    rtol=0,
                )
            relay.tidal.verify(decode(saved), saved["measurement"])
    print(
        data["paths"],
        {k: v["numerical_rank"] for k, v in data["source"]["scope"].items()},
    )
    print(data["measurement"], data["population_trace_distances"])
    print("Fixed source interval", fixed_score_read(data))
