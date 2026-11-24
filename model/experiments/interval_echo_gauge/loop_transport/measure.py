"""Native vacancy histories and sector-vacuum loop tests, with complete read records."""

import argparse
import hashlib
import json
from itertools import permutations
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, lobpcg, expm_multiply

from model.executable_model import cell, federation
from model.experiments.interval_echo_encounter.conjugation import defect_generator
from model.experiments.material_inertia.neutral_defect import generators, ordered_index

HERE = Path(__file__).resolve().parent
HOLDS = (0, 0.25, 0.5, 1, 2)


def faces(board):
    candidates = []
    for root in range(12):
        for a, b, c in permutations(board.neighbours[root], 3):
            if max(a, b, c) >= 12 or a >= c:
                continue
            if a in board.neighbours[b] and c in board.neighbours[b]:
                key = (
                    board.degree[root] != 5,
                    not any(board.degree[j] == 6 for j in (a, b, c)),
                    root,
                    a,
                    b,
                    c,
                )
                candidates.append((key, (root, a, b, c)))
    return min(candidates)[1]


def step_block(h, board, y, dest):
    # Read each native off-diagonal vacancy transition, including an excess at dest.
    other = [x for x in range(board.n) if x != y]
    after = [y if x == dest else x for x in other]
    src = ordered_index(np.array(other), y, board.n)
    dst = ordered_index(np.array(after), dest, board.n)
    block = h[dst][:, src].tocsr()
    assert (block + 0.5 * sparse.eye(len(src))).nnz == 0
    return dict(zip(other, after))


def paths(h, board, labels):
    root, a, b, c = labels
    f, g = (a, b, root), (b, c, root)
    records = {}
    maps = []
    for name, route in (("UV", f + g), ("VU", g + f)):
        pos = root
        order = [a, b, c]
        for dest in route:
            hop = step_block(h, board, pos, dest)
            order = [hop[j] for j in order]
            pos = dest
        maps.append(order)
        records[name] = dict(route=list(route), excess_from_abc=order)
    perm = np.eye(3)[:, [maps[1].index(z) for z in maps[0]]]
    assert np.trace(perm) == 0 and np.all(perm @ np.ones(3) == 1)
    return records


def flight(h, w, t):
    return expm_multiply((-1j * t) * h, w)


def measured_history(h, board, y, labels, path, t, rng):
    root, a, b, c = labels
    wave = np.zeros(h.shape[0], complex)
    wave[ordered_index(b, root, board.n)] = 1
    log = []
    survival = 1.0
    # Normalize each successful intermediate state to keep rare histories accurate.
    for vertex in path:
        wave = flight(h, wave, t)
        mask = y == vertex
        cond = float(np.sum(abs(wave[mask]) ** 2))
        log.append(survival * (1 - cond))
        survival *= cond
        wave[~mask] = 0
        wave /= np.sqrt(cond)
    prob = abs(wave) ** 2
    monitors = ordered_index(np.array([a, b, c]), root, board.n)
    other = np.ones(len(prob), bool)
    other[monitors] = False
    channels = np.r_[prob[monitors], prob[other].sum()]
    full = np.r_[log, survival * channels]
    assert abs(full.sum() - 1) < 2e-13
    return dict(
        success=survival,
        terminal_given_success=channels.tolist(),
        probabilities=full.tolist(),
        counts=rng.multinomial(2**42, full / full.sum()).tolist(),
    )


def loop_minus(w, board, root, a, b):
    result = np.zeros_like(w)
    u, v = ordered_index(a, root, board.n), ordered_index(b, root, board.n)
    result[u] = (w[u] - w[v]) / 2
    result[v] = -result[u]
    return result


def vacuum(h, board, x, y):
    if board.n == 12:
        energies, vectors = np.linalg.eigh(h.toarray())
        return float(energies[0]), vectors[:, 0]

    # The separable single-defect operator is a preconditioner only.
    # Its compression differs from the actual H by the adjacent-pair potential and swap.
    lap = board.laplacian().toarray() / 2
    ex, left = np.linalg.eigh(lap - np.diag(16 * board.degree))
    vy, right = np.linalg.eigh(lap + np.diag(16 * board.degree))
    inv = 1 / (ex[:, None] + vy[None, :] - (ex[0] + vy[0]) + 0.0001)

    def precondition(z):
        full = np.zeros((board.n, board.n))
        full[x, y] = z.ravel()
        return (left @ ((left.T @ full @ right) * inv) @ right.T)[x, y].reshape(z.shape)

    trial = left[x, 0] * right[y, 0]
    pre = LinearOperator(h.shape, matvec=precondition, dtype=float)
    e, v = lobpcg(h, trial[:, None], M=pre, largest=False, tol=2e-11, maxiter=650)
    w = v[:, 0] * np.sign(v[np.argmax(abs(v[:, 0])), 0])
    assert np.min(w) > 0
    # Preserve positivity while smoothing the tiny-amplitude tails. The extrema of
    # (Hw)/w then bracket the exact lowest eigenvalue by Perron comparison.
    for _ in range(1200):
        w += (e[0] * w - h @ w) / 80
    w /= np.linalg.norm(w)
    energy = float(w @ (h @ w))
    assert np.min(w) > 0 and np.linalg.norm(h @ w - energy * w) < 1e-9
    return energy, w


def vacuum_records(h, board, labels, w, energy, rng):
    root, a, b, c = labels
    p = loop_minus(w, board, root, a, b)
    expectation = float(np.vdot(p, p).real)
    v = p - expectation * w
    centered = h - energy * sparse.eye(h.shape[0], format="csr")
    correlations = [
        float(v @ expm_multiply(-t * centered, v)) if t else float(v @ v) for t in HOLDS
    ]
    assert abs(correlations[0] - expectation * (1 - expectation)) < 1e-13
    # The pulse is reflection in the selected antisymmetric excess readout.
    result = dict(
        mean=expectation,
        correlation=correlations,
        ground_residual=float(np.linalg.norm(h @ w - energy * w)),
        ground_bracket=[float(np.min((h @ w) / w)), float(np.max((h @ w) / w))],
        first_moment=float(v @ (centered @ v)),
    )
    if board.external:
        gate = next(j for j in (a, b, c) if board.degree[j] == 6)
        far = next(j for j in board.neighbours[gate] if j // 12 != gate // 12)
        rings = sorted(j for j in board.neighbours[far] if j // 12 == far // 12)
        f, g = next(
            (i, j) for i in rings for j in rings if i < j and j in board.neighbours[i]
        )
        reads = []
        sham = float(np.linalg.norm(loop_minus(w, board, far, f, g)) ** 2)
        for t in HOLDS:
            # Evolve only the local disturbance. Factor the vacuum phase out exactly.
            diff = expm_multiply((-1j * t) * centered, -2 * p) if t else -2 * p
            yes = float(np.linalg.norm(loop_minus(w + diff, board, far, f, g)) ** 2)
            reads.append(
                dict(
                    time=t,
                    sham=sham,
                    pulse=yes,
                    delta=yes - sham,
                    counts=[int(rng.binomial(2**24, x)) for x in (sham, yes)],
                )
            )
        result.update(receiver=[far, f, g], remote=reads)
    return result


def run():
    rng = np.random.default_rng(624701)
    boards = []
    for board in (cell(), federation(0)):
        x, y, _, kinetic, physical, _ = generators(board)
        labels = faces(board)
        for k, h in [(8, physical)] + ([(0, kinetic)] if board.n == 12 else []):
            routes = paths(h, board, labels)
            reads = []
            for t in (0.25, 0.5):
                for tag, r in routes.items():
                    reads.append(
                        dict(
                            order=tag,
                            hold=t,
                            **measured_history(h, board, y, labels, r["route"], t, rng),
                        )
                    )
            energy, w = vacuum(h, board, x, y)
            single = defect_generator(board, -1, k, 0)
            for record in routes.values():
                tour = [labels[0], *record["route"]]
                assert all(single[a, b] == -0.5 for a, b in zip(tour, tour[1:]))
            boards.append(
                dict(
                    ports=board.n,
                    sector=h.shape[0],
                    kappa=k,
                    labels=labels,
                    paths=routes,
                    history=reads,
                    ground_energy=energy,
                    single_vacancy_commutator=1,
                    marked_plane_half_commutator_trace=-0.5,
                    vacuum=vacuum_records(h, board, labels, w, energy, rng),
                )
            )
    return dict(
        protocol_sha256=hashlib.sha256((HERE / "DESIGN.md").read_bytes()).hexdigest(),
        seed=624701,
        history_trials=2**42,
        remote_trials=2**24,
        holds=HOLDS,
        boards=boards,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("x") as f:
        json.dump(run(), f, indent=2, allow_nan=False)
        f.write("\n")
