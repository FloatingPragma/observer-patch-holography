"""Native coincidence-triggered pair conversion between two polarizable wings."""

import argparse
from itertools import permutations, product
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import eigh, expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_20_electric.experiment import guide
from model.executable_model import federation

LAW = candidate.native_law
HERE = Path(__file__).resolve().parent
SHOTS = 2**25
TIMES = (0, 0.5, 2, 8, 32)
PAIRS = list(permutations(range(4), 2))
READY = [PAIRS.index((1, 2)), PAIRS.index((2, 1))]


def geometry():
    b = federation(0)
    _, path, _ = guide(b)
    pivot = next(
        i for i in range(5, len(path) - 5, 2) if path[i] // 12 != path[i + 1] // 12
    )
    src = path[pivot : pivot + 2]
    inc = [
        min(p for p in b.neighbours[x] if p not in path and p // 12 == x // 12)
        for x in src
    ]
    left = list(reversed(path[1:pivot]))
    right = path[pivot + 2 : -1]
    return b, [inc[0], *src, inc[1]], [left, right]


def word(values, extra, empty):
    a = np.ones(values, dtype=int)
    a[extra], a[empty] = 2, 0
    return a


def term(k):
    x, c = LAW.local_term_parts()
    return (x + k * c)[:9, :9]


def collision(k, board, ports):
    mat = np.zeros((12, 12))
    local = term(k)
    for n, (p, q) in enumerate(PAIRS):
        w = word(4, p, q)
        mat[n, n] = k * ((6 - board.degree[ports]) @ (w * w - 1))
        for i in range(3):
            a, b = w[i : i + 2]
            col = LAW.SRC.index((a, b))
            mat[n, n] += local[col, col]
            if a != b:
                new = w.copy()
                new[i : i + 2] = [b, a]
                pos = (np.flatnonzero(new == 2).item(), np.flatnonzero(new == 0).item())
                mat[PAIRS.index(pos), n] += local[LAW.SRC.index((b, a)), col]
    return mat


def wing(k, n, inner, cut=False):
    # All dimers have matched native degree. Their degree cost is a common constant.
    # Energies below subtract the all-unit central reference and the same dimer constants.
    rows = list(product((0, 1), repeat=n))
    local = term(k)
    h = np.zeros((len(rows), len(rows)))
    for j, r in enumerate(rows):
        w = [v for bit in r for v in (2 * bit, 2 - 2 * bit)]
        for i in range(0, 2 * n, 2):
            a, b = w[i : i + 2]
            col = LAW.SRC.index((a, b))
            other = rows.index(
                tuple(1 - v if z == i // 2 else v for z, v in enumerate(r))
            )
            h[other, j] = local[LAW.SRC.index((b, a)), col]
        for i, (a, b) in enumerate(zip([inner] + w, w)):
            if cut and i == 2:
                continue
            col = LAW.SRC.index((a, b))
            h[j, j] += local[col, col]
    return h


def basis(k, n):
    board, ports, _ = geometry()
    coll = collision(k, board, ports)
    hs = [wing(k, n, z) for z in (0, 1, 2)]
    e, v = eigh(hs[1])
    g = v[:, 0]
    eye = sparse.eye(len(g), format="csr")
    blocks = []
    for p, q in PAIRS:
        a, b = word(4, p, q)[1:3]
        blocks.append(sparse.kron(hs[a], eye) + sparse.kron(eye, hs[b]))
    incoming = sparse.kron(coll, sparse.eye(len(g) ** 2)) + sparse.block_diag(blocks)
    return incoming.tocsr(), hs, coll, e, g


def launch(n, sign, h, g):
    v = np.zeros((12, len(g) ** 2), complex)
    source = PAIRS.index((0, 3) if sign == 1 else (3, 0))
    v[source] = np.kron(g, g)
    shift = h.diagonal().mean()
    generator = -0.25j * (h - shift * sparse.eye(h.shape[0]))
    arrival = []
    v = v.ravel()
    for _ in range(16):
        v = expm_multiply(generator, v)
        parts = v.reshape(12, -1)
        arrival.append(parts[READY].copy())
        parts[READY] = 0
    return np.array(arrival), float(np.vdot(v, v).real)


def effects(n):
    dim = 2**n
    # Anti-swap on the outermost dimer equals its variable exchange energy.
    return np.kron(np.eye(dim // 2), np.array([[0.5, -0.5], [-0.5, 0.5]]))


def read_joint(v, left, right):
    """Rows of v are unnormalized conditional pure states of left and right wings."""
    w = v.reshape(-1, len(left), len(right))
    norm = float(np.vdot(w, w).real)
    lp = left @ w
    rp = w @ right.T
    pleft = float(np.vdot(w, lp).real)
    pright = float(np.vdot(w, rp).real)
    both = float(np.vdot(lp, rp).real)
    return np.array([norm - pleft - pright + both, pleft - both, pright - both, both])


def probabilities(arrival, ops, read, t):
    total = np.zeros(4)
    for i in (0, 1):
        eig = []
        for h in ops[i]:
            e, v = eigh(h)
            u = (v * np.exp(-1j * t * (e - e[0]))) @ v.T
            eig.append(u.conj().T @ read @ u)
        total += read_joint(arrival[:, i], *eig)
    return total


def sample(prob, rng):
    p = np.maximum(0, np.asarray(prob))
    # Stable rounding for the seeded multinomial, with the largest bin as remainder.
    order = np.argsort(p.round(10), kind="stable")
    q = p[order].round(11)
    q[-1] = 1 - q[:-1].sum()
    assert all(q >= 0) and abs(p.sum() - 1) < 1e-9
    counts = np.empty(len(q), np.int64)
    counts[order] = rng.multinomial(SHOTS, q)
    return dict(probabilities=p.round(11).tolist(), counts=counts.tolist())


def stats(k, n, sign, rng):
    h, hs, coll, e, g = basis(k, n)
    arrival, miss = launch(n, sign, h, g)
    hit = 1 - miss
    dims = len(g)
    neutral = hs[1]
    post = [(neutral, neutral)] * 2
    parked = [(hs[2], hs[0]), (hs[0], hs[2])]
    cut = [(wing(k, n, 1, True), wing(k, n, 1, True))] * 2
    ground = np.eye(dims) - np.outer(g, g)
    abs_energy = np.zeros(2)
    budget = []
    for i, idx in enumerate(READY):
        psi = arrival[:, i].reshape(-1, dims, dims)
        w = float(np.vdot(psi, psi).real)

        def mean(left, right):
            return float(np.vdot(psi, left @ psi + psi @ right.T).real)

        after = mean(neutral, neutral)
        before = mean(*parked[i]) + w * coll[idx, idx]
        budget.append((before - after) / hit)
        abs_energy += (
            np.array(
                [
                    mean(neutral - e[0] * np.eye(dims), 0 * neutral),
                    mean(0 * neutral, neutral - e[0] * np.eye(dims)),
                ]
            )
            / hit
        )
    opportunity = np.sum(abs(arrival) ** 2, axis=2)
    energy = read_joint(arrival.reshape(-1, dims * dims), ground, ground)
    record = dict(
        contrast=k,
        dimers_per_wing=n,
        launch_sign=sign,
        opportunity=sample(np.r_[opportunity.ravel(), miss], rng),
        wing_excitation=sample(np.r_[energy, miss], rng),
        wing_energy_above_ground_given_hit=abs_energy.tolist(),
        pair_energy_removed_per_hit=float(sum(budget)),
        detectors=[],
    )
    for arm, ops in (("converted", post), ("parked", parked), ("converted_cut", cut)):
        for t in TIMES:
            p = probabilities(arrival, ops, effects(n), t)
            record["detectors"].append(
                dict(arm=arm, t=t, **sample(np.r_[p, miss], rng))
            )
    cut_rng = np.random.default_rng(81127 + 12 * k + sign + 4 * n)
    parked_cut = [
        (wing(k, n, 2, True), wing(k, n, 0, True)),
        (wing(k, n, 0, True), wing(k, n, 2, True)),
    ]
    for t in TIMES:
        p = probabilities(arrival, parked_cut, effects(n), t)
        record["detectors"].append(
            dict(arm="parked_cut", t=t, **sample(np.r_[p, miss], cut_rng))
        )
    return record


def uniform(k, sign, rng):
    board, ports, _ = geometry()
    coll = collision(k, board, ports)
    src = PAIRS.index((0, 3) if sign == 1 else (3, 0))
    psi = np.eye(12, dtype=complex)[:, src]
    u = expm(-0.25j * coll)
    arriving = 0.0
    for _ in range(16):
        psi = u @ psi
        arriving += sum(abs(psi[READY]) ** 2)
        psi[READY] = 0
    return dict(
        contrast=k,
        launch_sign=sign,
        extinction_and_either_remote_excitation_zero=sample(
            [arriving, 0, 1 - arriving], rng
        ),
    )


def run():
    rng = np.random.default_rng(81124)
    board, ports, wings = geometry()
    return dict(
        seed=81124,
        shots=SHOTS,
        board=board.cells,
        ports=ports,
        wing_routes=[w[:6] for w in wings],
        uniform=[uniform(k, s, rng) for k in (8, 0) for s in (-1, 1)],
        cases=[stats(k, n, s, rng) for k in (8, 0) for n in (2, 3) for s in (-1, 1)],
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = json.loads(json.dumps(run()))
    path = HERE / "records.json"
    if args.write:
        path.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        from .verify import compare

        compare(data, json.loads(path.read_text()), "collision")
        print(
            "Verified all encounter, work, neutral-wing energy and remote coincidence records."
        )


if __name__ == "__main__":
    main()
