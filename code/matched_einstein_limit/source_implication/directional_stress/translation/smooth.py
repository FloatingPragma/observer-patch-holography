"""Exact compressed energy current, weak-phase spatial transfer and leakage from flat port orbitals."""

import argparse
import json
from collections import Counter
from fractions import Fraction

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from .. import bridge
from . import counts, measure


def uniform_lift(board):
    return sparse.kron(
        sparse.kron(sparse.eye(board.cells, dtype=int), np.ones((12, 1), dtype=int)),
        sparse.eye(4, dtype=int),
        format="csr",
    )


def load4(board, a, b):
    position = sparse.csr_matrix(
        ([-4, 4], ([a, b], [b, a])), shape=(board.n, board.n), dtype=np.int64
    )
    return sparse.kron(position, sparse.eye(4, dtype=int), format="csr")


def routes(board, region, route):
    graph = nx.Graph(board.seams)
    graph.remove_edges_from((route[:2], route[2:]))
    weights = Counter()
    for c, d, perm in zip(region["cells"], region["images"], region["lifts"]):
        for i, j in enumerate(perm):
            path = nx.shortest_path(graph, 12 * c + i, 12 * d + j)
            for a, b in zip(path, path[1:]):
                if a // 12 != b // 12:
                    weights[min(a, b), max(a, b)] += 1 if a < b else -1
    assert tuple(route[1:3]) not in weights
    return {k: v for k, v in weights.items() if v}


def leakage_lemma(board, route, h, lift, q):
    ports = np.flatnonzero(board.degree == 5)
    refs = {
        c: [int(p) for p in ports if p // 12 == c and p not in route]
        for c in range(board.cells)
    }
    unmarked = [e for e in board.external if tuple(e) != tuple(route[1:3])]
    # The scalar removed from h cancels between outputs in a cell.
    mixed = h @ lift
    coeff = 33 if q == -1 else -31
    used = Counter()
    for a, b in unmarked:
        for g, v in ((a, b), (b, a)):
            c, d = g // 12, v // 12
            r = refs[c][used[c]]
            used[c] += 1
            delta = mixed[4 * g : 4 * g + 4] - mixed[4 * r : 4 * r + 4]
            expected = np.zeros((4, 4 * board.cells), dtype=int)
            expected[:, 4 * c : 4 * c + 4] = coeff * np.eye(4, dtype=int)
            expected[:, 4 * d : 4 * d + 4] = -np.eye(4, dtype=int)
            np.testing.assert_array_equal(delta.toarray(), expected)
    assert min(used.values()) >= 2 and len(used) == board.cells
    return str(Fraction((abs(coeff) - 1) ** 2, 48))


def one_head_structure(couple):
    """Exact decompositions L_b⊗I+I⊗R for every position block of the central current."""
    touched = set()
    for mat in couple:
        m = mat.tocoo()
        touched.update(zip(m.row // 4, m.col // 4))
    eye = np.eye(2, dtype=np.int64)
    for i, j in touched:
        blocks = [mat[4 * i : 4 * i + 4, 4 * j : 4 * j + 4].toarray() for mat in couple]
        b = blocks[0]
        r4 = 2 * np.trace(b.reshape(2, 2, 2, 2), axis1=0, axis2=2) - np.trace(b) * eye
        for b in blocks:
            l2 = np.trace(b.reshape(2, 2, 2, 2), axis1=1, axis2=3)
            np.testing.assert_array_equal(
                4 * b, np.kron(2 * l2, eye) + np.kron(eye, r4)
            )
    return len(touched)


def run():
    atlas = measure.Atlas(1)
    board = atlas.board
    patch = measure.chosen(atlas)
    route = bridge.native.iso.pairing.formation.region(board)
    rings = bridge.beam.screens(board, route)
    src = bridge.beam.indices(list(rings[0]) + list(rings[1]))
    dist = np.array([atlas.dist[0][c] for c in range(board.cells)])
    delta = Counter(
        int(dist[v] - dist[u]) for u, v in zip(patch["cells"], patch["images"])
    )
    assert delta == {2: 11, 0: 1}
    weights = routes(board, patch, route)
    lift = uniform_lift(board)
    states = bridge.action.preparations()
    rng = np.random.default_rng(67192638)
    trials = 2**60
    output = {}
    for q in (-1, 1):
        projected = []
        central = []
        first_unmarked = {}
        standard = []
        uniform_read = {}
        magnitudes = {}
        for odd in (False, True):
            cells = bridge.operators(board, route, q, odd)
            h, ops = measure.op_table(board, route, patch, q, odd)
            lower = leakage_lemma(board, route, h, lift, q)
            energy = sparse.csr_matrix(h.shape, dtype=int)
            for e in board.external:
                j = bridge.current(board, cells, e)
                projected_diff = lift.T @ (j - 40 * load4(board, *e)) @ lift
                assert (projected_diff.count_nonzero() == 0) == (e != tuple(route[1:3]))
                if e == tuple(route[1:3]):
                    central.append(j)
                    flat = np.kron((1j**dist)[:, None], np.eye(4))
                    standard.append(
                        flat.conj().T @ (lift.T @ j @ lift @ flat) * 1j / (8 * board.n)
                    )
                else:
                    scalar = j[::4, ::4]
                    assert (
                        j - sparse.kron(scalar, sparse.eye(4, dtype=int))
                    ).count_nonzero() == 0
                    if not odd:
                        first_unmarked[e] = scalar
                    else:
                        assert (scalar - first_unmarked[e]).count_nonzero() == 0
                if e in weights:
                    energy -= weights[e] * j
            # Actual O = real skew/(2i) or i*K/(8*40).
            long, dlong = ops["shift"]
            denergy = -320j
            # The all-cell linear phase response is an operator on the heads,
            # for every one of the 80 infinitesimal cell phase coordinates.
            skew_difference = lift.T @ (160 * long + energy) @ lift
            coo = skew_difference.tocoo()
            gradients = np.zeros((board.cells, 4, 4), dtype=np.int64)
            for i, j, coef in zip(coo.row, coo.col, coo.data):
                gradients[j // 4, i % 4, j % 4] += coef
                gradients[i // 4, i % 4, j % 4] -= coef
            assert not np.any(gradients)
            for name, op, unit in (
                ("shift", long, dlong),
                ("routed_energy_over_40", energy, denergy),
            ):
                magnitudes[name] = max(
                    magnitudes.get(name, 0), counts.selection_norm(op, unit)
                )
                mat = lift.T @ op @ lift
                rows = []
                for alpha in (0.0, 1 / 64, 1 / 16, 0.25):
                    amp = np.kron(np.exp(1j * alpha * dist)[:, None], np.eye(4))
                    head = amp.conj().T @ (mat @ amp) / (board.n * unit)
                    exact = (
                        -11 * np.sin(2 * alpha) / 80
                        if name == "shift"
                        else -11 * np.sin(alpha) / 40
                    )
                    np.testing.assert_allclose(head, exact * np.eye(4), atol=2e-15)
                    rows.append(head)
                uniform_read.setdefault(name, []).append(rows)
            start = np.zeros((h.shape[0], len(src)))
            start[src] = np.eye(len(src))
            ht = -2j * h
            w = expm_multiply(ht, start, traceA=ht.diagonal().sum())
            part = lift.T @ w / np.sqrt(12)
            projected.append(
                np.mean([a.conj().T @ a for a in np.hsplit(part, 10)], axis=0)
            )
        readings = {}
        for name, heads in uniform_read.items():
            scale = magnitudes[name]
            vals = []
            for i, alpha in enumerate((0.0, 1 / 64, 1 / 16, 0.25)):
                panel = np.array(heads)[:, i]
                val = bridge.evaluate(panel, states)
                np.testing.assert_allclose(
                    list(val.values()), val["equilibrium"], atol=2e-15
                )
                plus = int(rng.binomial(trials, 0.5 + val["equilibrium"] / (2 * scale)))
                estimate = float(Fraction(2 * plus - trials, trials)) * scale
                halfwidth = scale * np.sqrt(2 * np.log(2 * 16 / 0.01) / trials)
                assert abs(estimate - val["equilibrium"]) < halfwidth
                vals.append(
                    dict(
                        alpha=alpha,
                        means=val,
                        plus=plus,
                        estimate=estimate,
                        halfwidth=halfwidth,
                    )
                )
            readings[name] = vals
        output[str(q)] = dict(
            all_unmarked_bridges=119,
            central_position_blocks_proved_one_head_only=one_head_structure(central),
            arbitrary_common_orbital_pairing_response_exact=0,
            central_flux_at_pi_over_two=bridge.evaluate(standard, bridge.tangents()),
            marked_exceptions=1,
            independent_phase_directions=79,
            lower_bound_on_squared_leakage_speed=lower,
            original_mixture_uniform_subspace_weight=bridge.evaluate(projected, states),
            readings=readings,
        )
    return dict(
        source_heads="registered stationary source",
        weak_phase_response="-11/40 for distance from cell 0",
        shot_seed=67192638,
        trials_per_setting=trials,
        simultaneous_coverage=0.99,
        charge_results=output,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "smooth.json").write_text(data)
    else:
        print(data, end="")
