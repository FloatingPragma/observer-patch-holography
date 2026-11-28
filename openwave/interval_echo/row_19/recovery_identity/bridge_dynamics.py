"""Independent propagation and population predictions for the archive and contact tests."""

import argparse
import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.special import gammaln, jv

from . import (
    bridge_radar as bridge,
    paired_relay as pair,
    paired_survey as radar,
    relay,
    stored_source as store,
)

RECEIPT = bridge.HERE / "bridge_dynamics.json"


def flight(board, route, paths, q, parity):
    """Compile from the two-sender generator and propagate by a Bessel expansion."""
    h2 = relay.RANGING.generator(board, route, q, 8, bool(parity))
    h = sparse.kron(h2, sparse.eye(4), format="csr")
    native = relay.FORMATION.local_parts(q, 8)
    for slot, path in enumerate(paths):
        edge = sorted(path[:2])
        for i, j in np.ndindex((2, 2)):
            change = native[i + 1, j + 1].copy()
            if i == j:
                change -= native[0, 0] + 8.5 * np.eye(2)
            else:
                change += 0.5 * np.eye(2)
            active = np.kron(
                np.eye(4 * 2**slot), np.kron(change, np.eye(2 ** (1 - slot)))
            )
            addr = sparse.csr_matrix(
                ([1], ([edge[i]], [edge[j]])), shape=(board.n, board.n)
            )
            h += sparse.kron(addr, active, format="csr")
    h -= h.diagonal().mean() * sparse.eye(h.shape[0], format="csr")
    radius = float(abs(h).sum(axis=1).max())
    x = h / radius
    z = 8 * radius
    order = int(np.ceil(1.6 * z + 55))
    log_bound = (
        (order + 1) * np.log(z / 2)
        - gammaln(order + 2)
        + np.log(2 / (1 - z / (2 * (order + 2))))
    )
    assert log_bound < -39
    coeff = 2 * (-1j) ** np.arange(order + 1) * jv(np.arange(order + 1), z)
    a = np.zeros((h.shape[0], 4))
    a[16 * route[0] + 4 * np.arange(4)] = np.eye(4)
    b = x @ a
    u = coeff[0] * a / 2 + coeff[1] * b
    for c in coeff[2:]:
        a, b = b, 2 * (x @ b) - a
        u += c * b
    return u.reshape(board.n, 4, 4, 4)


def return_operators(board, paths, ports, probe):
    """Project the full propagated wave off the detector after each missed arrival."""
    contacts = [p for path in paths for p in path[:2]]
    h = relay.RANGING.generator(board, contacts, probe, 8, False)
    n = h.shape[0]
    first = {}
    for target in range(3):
        senders = [s for s in range(3) if s != target]
        wave = np.zeros((n, 8))
        for i, s in enumerate(senders):
            p = ports[s]
            wave[4 * p : 4 * p + 4, 4 * i : 4 * i + 4] = np.eye(4)
            first[s, target] = [np.zeros((4, 4))]
        p = ports[target]
        for _ in range(10):
            wave = expm_multiply(-0.5j * h, wave, traceA=0)
            arrived = wave[4 * p : 4 * p + 4]
            for i, s in enumerate(senders):
                first[s, target].append(arrived[:, 4 * i : 4 * i + 4].copy())
            arrived[:] = 0
    means = []
    for i, j in radar.ROUTES:
        mean = 3 * np.eye(4, dtype=complex)
        for t in range(1, 11):
            for u in range(1, 12 - t):
                amplitude = first[j, i][u] @ first[i, j][t]
                mean += (t + u - 12) / 4 * (amplitude.conj().T @ amplitude)
        means.append(mean)
    return np.array(means)


def geometry(op, states):
    times = np.einsum("rij,sji->sr", op, states).real
    lane = {r: row for r, row in zip(radar.ROUTES, times.T, strict=True)}
    a, b, c = ((lane[i, j] + lane[j, i]) / 2 for i, j in ((0, 1), (1, 2), (0, 2)))
    area = np.sqrt((a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c)) / 4
    return times, area


def calculate():
    board = relay.FORMATION.federation(1)
    route = relay.FORMATION.region(board)
    arrival = json.loads(pair.RECEIPT.read_text())
    paths = arrival["paths"]
    stocks, steps, _, _ = relay.sources(board, route)
    accepted = np.array([store.measure.product_reference(s) for s in stocks[:4]])
    np.testing.assert_allclose(accepted[2], accepted[3], atol=1e-14, rtol=0)
    # These two products coincide algebraically at fixed one-head marginals.
    accepted[3] = accepted[2]
    full = np.concatenate((stocks, accepted))
    ref, *_ = store.measure.head.jacobi.beam.aperture.source_and_changes()
    k = -store.measure.head.station.log_blocks(ref)
    result, conversions, balanced = {}, {}, {}
    operators = {
        (i, p): return_operators(board, paths, path, p)
        for i, path in enumerate(paths)
        for p in (-1, 1)
    }
    for q in (-1, 1):
        received = np.zeros((len(full), 4, 4), complex)
        contact = np.zeros((len(full), 2))
        for b in (0, 1):
            u = flight(board, route, paths, q, b)
            v = u.reshape(-1, 4, 4)
            reduced = np.einsum("pca,pdb->cdab", v, v.conj(), optimize=True)
            received += np.einsum("cdab,sab->scd", reduced, full[:, b])
            endpoints = u[list(route[1:3])].reshape(2, -1, 4)
            # Real part of the endpoint coherence is the two-port interferometer.
            diag = sum(x.conj().T @ x for x in endpoints)
            off = endpoints[0].conj().T @ endpoints[1]
            coherence = off + off.conj().T
            for j, sgn in enumerate((1, -1)):
                contact[:, j] += (
                    0.5
                    * np.einsum("ij,sji->s", diag + sgn * coherence, full[:, b]).real
                )
        np.testing.assert_allclose(
            received[:8],
            pair.unpack(arrival["installed"][int(q == 1)]),
            atol=4e-12,
            rtol=0,
        )
        change = np.diff(contact[:8].reshape(4, 2, 2), axis=1)[:, 0, :]
        work = change @ np.array([8.0, 9.0])
        ratio = work[1] / work[0] * steps["ordinary"] / steps["held"]
        result[str(q)] = dict(
            probabilities=contact[:8].tolist(), work_ratio=float(ratio)
        )
        for (i, p), op in operators.items():
            vals, areas = geometry(op, received[:8])
            # Pure fixed-marginal endpoints exclusively. Balance using the
            # model probabilities; check the two resulting mixed quantum states.
            coeff = np.cross(change[1:, 0], change[1:, 1])
            coeff /= abs(coeff).sum()
            indices = np.array([[3, 2], [5, 4], [7, 6]])
            pairs = np.where((coeff < 0)[:, None], indices[:, ::-1], indices)
            weights = abs(coeff)
            mix = np.einsum("k,ktbij->tbij", weights, full[pairs])
            for s in range(2):
                np.testing.assert_allclose(
                    store.measure.head.marginals(mix[0])[s],
                    store.measure.head.marginals(mix[1])[s],
                    atol=1e-14,
                )
            matched = np.einsum("k,ktm->tm", weights, contact[pairs])
            np.testing.assert_allclose(matched[0], matched[1], atol=1e-14)
            mix_received = np.einsum("k,ktij->tij", weights, received[pairs])
            # Cancel the common identity term analytically for a normalized contrast.
            centered = op - 3 * np.eye(4)
            dr = np.einsum(
                "rij,ji->r", centered, mix_received[0] - mix_received[1]
            ).real
            balanced[f"{i},{q},{p}"] = dict(
                weights=coeff.tolist(),
                matched=matched.tolist(),
                range_differences=dr.tolist(),
                areas=areas.tolist(),
                ranges=vals.tolist(),
            )
            if i != 0:
                continue
            for name, f in store.fractions().items():
                state = f * full[:4] + (1 - f) * full[8:]
                output = f * received[:4] + (1 - f) * received[8:]
                _, area = geometry(op, output)
                modular = np.einsum("bij,sbji->s", k, state).real
                mu = store.two_slopes(modular, steps)
                g = store.two_slopes(area, steps)
                conversions.setdefault(name, {})[f"{q},{p}"] = dict(
                    source=list(mu), geometry=list(g)
                )
    return dict(contact=result, balance=balanced, archive=conversions)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    d = calculate()
    if args.write:
        with RECEIPT.open("x") as f:
            json.dump(d, f, separators=(",", ":"), allow_nan=False)
            f.write("\n")
    else:
        relay.tidal.verify(d, json.loads(RECEIPT.read_text()))
    print("Bridge work ratios", [v["work_ratio"] for v in d["contact"].values()])
    for name, v in d["archive"].items():
        src, geo = v["1,1"].values()
        print(name, src[1] / src[0], geo[1] / geo[0])
