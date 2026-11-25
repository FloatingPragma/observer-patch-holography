"""Ground-sector and charge-resolved endpoint readouts on physical open blank paths."""

import argparse
import json
from itertools import pairwise
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import eigh, eigh_tridiagonal, expm
from scipy.sparse import csr_matrix, eye, kron

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation

LAW = candidate.native_law
K = candidate.KAPPA
HERE = Path(__file__).resolve().parent
SHOTS = 2**22
SEED = 948473


def bonds():
    ex, con = LAW.local_term_parts()
    native = ex + K * con
    h = native[:9, :9]
    np.testing.assert_allclose(native[9:, :9], 0, atol=1e-13)
    np.testing.assert_allclose(h, LAW.seam_term(K), atol=1e-13)
    np.testing.assert_allclose(
        eigh(h, eigvals_only=True), [0, 0, 0, 8, 8, 9, 9, 32, 33], atol=1e-12
    )
    raising = np.diag(np.sqrt([2.0, 2.0]), -1)
    sx = (raising + raising.T) / 2
    sy = (raising - raising.T) / (2j)
    sz = np.diag([-1.0, 0, 1])
    dot = sum(np.kron(s, s) for s in [sx, sy, sz]).real
    p2 = (dot + 2 * np.eye(9)) @ (dot + np.eye(9)) / 6
    np.testing.assert_allclose(p2 @ p2, p2, atol=1e-13)
    assert np.linalg.matrix_rank(p2, tol=1e-10) == 5
    for axis in [sx, sy]:
        r = expm(-1j * np.pi * axis)
        rr = np.kron(r, r)
        for bond in [h, p2]:
            np.testing.assert_allclose(rr @ bond, bond @ rr, atol=1e-12)
    return h, p2, sx, sy


def material(length):
    car = federation(0)
    graph = nx.Graph(car.seams)
    dist = nx.single_source_shortest_path_length(graph, 0)
    end = min(p for p, d in dist.items() if d == max(dist.values()))
    path = min(nx.all_shortest_paths(graph, 0, end))[:length]
    assert len(path) == length
    active = {tuple(sorted(e)) for e in pairwise(path)}
    assert len(nx.Graph(graph.subgraph(path)).edges) == len(active)
    unused = [
        e for e in car.seams if set(e) & set(path) and tuple(sorted(e)) not in active
    ]
    # Dark cuts and erased private references are fixtures. A zero-mark even
    # collar shuts every incoming load column, including arbitrary entanglement.
    dark = np.zeros((99, 9))
    for bit in (0, 1):
        dark[
            (1 + LAW.LABELS.index((bit, 0))) * 9 : (2 + LAW.LABELS.index((bit, 0))) * 9
        ] = np.eye(9) / np.sqrt(2)
    ex, con = LAW.local_term_parts()
    np.testing.assert_allclose((ex + K * con) @ dark, 0, atol=1e-13)
    return dict(
        ports=path,
        shutters=len(unused),
        private_erasures=sum(int(car.degree[p] == 5) for p in path),
    )


def chain(bond, L):
    return sum(
        kron(kron(eye(3**i), csr_matrix(bond)), eye(3 ** (L - i - 2)), format="csr")
        for i in range(L - 1)
    )


def valence_frame(L):
    # Symmetric pair of spin-1/2 indices in the physical load order 0,1,2,
    # contracted with an oriented singlet. No endpoint turn matrices are inputs.
    triplet = np.array(
        [[[1, 0], [0, 0]], [[0, 1 / np.sqrt(2)], [1 / np.sqrt(2), 0]], [[0, 0], [0, 1]]]
    )
    tensor = np.sqrt(2 / 3) * triplet @ np.array([[0, 1], [-1, 0]])
    frame = np.eye(2)[None, :, :]
    for _ in range(L):
        frame = np.einsum("wab,ibc->wiac", frame, tensor).reshape(-1, 2, 2)
    return frame


def axes():
    return np.array([[1, 1j, 0], [0, 0, np.sqrt(2)], [-1, 1j, 0]]) / np.sqrt(2)


def endpoint(L, p2, sx, sy):
    # The two ordinary right-boundary reference labels close the finite wire.
    # The left logical input has unit physical Gram on an arbitrary superposition.
    frame = valence_frame(L)

    def wire(n):
        return valence_frame(n).transpose(0, 2, 1).reshape(-1, 2)

    w = wire(L)
    shorter = wire(L - 1)
    np.testing.assert_allclose(w.T @ w, np.eye(2), atol=1e-12)
    np.testing.assert_allclose(shorter.T @ shorter, np.eye(2), atol=1e-12)
    raw = frame.reshape(-1, 4)
    np.testing.assert_allclose(chain(p2, L) @ raw, 0, atol=1e-12)
    logical = []
    for a, S in enumerate([sx, sy]):
        physical = axes()[:, a].conj() @ w.reshape(3, -1)
        physical = physical.reshape(-1, 2)
        gate = np.sqrt(3) * shorter.T @ physical
        np.testing.assert_allclose(shorter @ gate, np.sqrt(3) * physical, atol=1e-12)
        np.testing.assert_allclose(gate.conj().T @ gate, np.eye(2), atol=1e-12)
        global_turn = expm(-1j * np.pi * S)
        transformed = raw.reshape((3,) * L + (4,))
        for j in range(L):
            transformed = np.moveaxis(
                np.tensordot(global_turn, transformed, axes=(1, j)), 0, j
            )
        np.testing.assert_allclose(
            transformed.reshape(-1, 4), raw @ np.kron(gate, gate.conj()), atol=1e-12
        )
        logical.append(gate)
    x, y = logical
    comm = x @ y @ x.conj().T @ y.conj().T
    np.testing.assert_allclose(comm, -np.eye(2), atol=1e-12)
    return float(np.trace(comm).real / 2)


def clock_weights(energies, center):
    # Coherent equally spaced delays followed by a uniform clock erasure.
    # The total energy spread lies in one unaliased period in these fixtures.
    dt = 1 / 32
    assert np.ptp(energies) < 2 * np.pi / dt
    theta = (energies - center) * dt / (2 * np.pi)
    return (np.sinc((2**28) * theta) / np.sinc(theta)) ** 2


def bank(columns, rng, n):
    # Each source column carries its admitted amplitude.
    f = columns.reshape(3, 3, -1, columns.shape[-1])
    rank = np.linalg.norm(columns) ** 2
    # Charge-twirled Cartesian xy+ test and coherent xy− test expressed
    # entirely in load records. Neither requires a coherent 00/22 superposition.
    bright = (np.linalg.norm(f[0, 0]) ** 2 + np.linalg.norm(f[2, 2]) ** 2) / (2 * rank)
    dark = np.linalg.norm(f[0, 2] - f[2, 0]) ** 2 / (2 * rank)
    probs = np.array([bright, dark, 1 - bright - dark])
    outcomes = rng.multinomial(n, np.maximum(probs, 0))
    return dict(
        bright=float(bright),
        dark=float(dark),
        other=float(probs[-1]),
        counts=outcomes.tolist(),
    )


def run():
    rng = np.random.default_rng(SEED)
    clock_rng = np.random.default_rng(SEED + 1)
    h, p2, sx, sy = bonds()
    # A native pi blank contact is SWAP; joint unordered load read after a
    # balanced time-branch parity analyzer picks the antisymmetric 02 event.
    np.testing.assert_allclose(expm(-1j * np.pi * h), LAW._load_ops()[0], atol=1e-12)
    results = []
    for L in (3, 4, 5, 6):
        row = dict(material=material(L), length=L)
        for key, bond in (("registered_blank", h), ("AKLT", p2)):
            H = chain(bond, L)
            # The full eigensolve is independent of the analytic ground frame.
            ev, vectors = eigh(H.toarray(), driver="evd")
            columns = vectors[:, abs(ev) < 1e-8]
            dim = columns.shape[1]
            assert dim == (3 if key == "registered_blank" else 4)
            p = dim / (3**L)
            admitted = int(rng.binomial(SHOTS, p))
            row[key] = dict(
                ground_dimension=dim,
                ground_admission=p,
                preparation_refused=SHOTS - admitted,
                endpoint_bank=bank(columns, rng, admitted),
            )
            weights = clock_weights(ev, 0.0)
            p_clock = weights.mean()
            accepted = int(clock_rng.binomial(SHOTS, p_clock))
            row[key]["finite_clock"] = dict(
                admission=float(p_clock),
                refused=SHOTS - accepted,
                excited_fraction=float(weights[abs(ev) >= 1e-8].sum() / weights.sum()),
                endpoint_bank=bank(vectors * np.sqrt(weights), clock_rng, accepted),
            )
            if key == "AKLT":
                row[key]["edge_commutator"] = endpoint(L, p2, sx, sy)
            else:
                for q in [-1, 1]:
                    words = [
                        int(
                            np.dot(
                                np.array([1 + q if i == j else 1 for i in range(L)]),
                                3 ** np.arange(L - 1, -1, -1),
                            )
                        )
                        for j in range(L)
                    ]
                    sub = H[words, :][:, words].toarray()
                    diag = (K + 0.5) * np.array([1] + [2] * (L - 2) + [1])
                    expected = (
                        np.diag(diag)
                        + np.diag([-0.5] * (L - 1), 1)
                        + np.diag([-0.5] * (L - 1), -1)
                    )
                    np.testing.assert_allclose(sub, expected, atol=1e-12)
                    vals, v = eigh_tridiagonal(diag, [-0.5] * (L - 1))
                    endpoint_rows = v[0:1, :2]
                    n = int(rng.binomial(SHOTS, v[0, 0] ** 2))
                    wt = clock_weights(vals, vals[0])
                    accepted = int(clock_rng.binomial(SHOTS, wt.mean()))
                    p_left = np.dot(wt, v[0] ** 2) / sum(wt)
                    hit = int(clock_rng.binomial(accepted, p_left))
                    row[key][str(q)] = dict(
                        clock_target=vals[0],
                        clock_admission=float(wt.mean()),
                        clock_excited_fraction=float(wt[1:].sum() / wt.sum()),
                        clock_left_probability=float(p_left),
                        clock_refusals=SHOTS - accepted,
                        clock_arrivals=hit,
                        clock_other=accepted - hit,
                        lowest_energy=vals[0],
                        second_energy=vals[1],
                        orbital_split=vals[1] - vals[0],
                        same_end_arrival_rank=int(np.linalg.matrix_rank(endpoint_rows)),
                        left_arrival_probability=v[0, 0] ** 2,
                        arrivals=n,
                        other_positions=SHOTS - n,
                    )
        results.append(row)
    return dict(
        starts_per_bank=SHOTS,
        seed=SEED,
        clock_seed=SEED + 1,
        clock_delays=2**28,
        clock_tick=1 / 32,
        chains=results,
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    replay = run()
    if ap.parse_args().write:
        (HERE / "readings.json").write_text(json.dumps(replay, indent=2) + "\n")
    else:
        print(json.dumps(replay, indent=2))
