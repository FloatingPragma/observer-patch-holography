"""Search for a spin-order response of a charge prepared with a coherent field circulation."""

import argparse
from itertools import combinations, product
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
# The archive replay launches this script by filename from a temporary directory.
# ruff: noqa: E402
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (
    measure as tree,
)


class Circulation:
    def __init__(self, fixture, board, spur=False):
        self.fixture = fixture
        self.bits = fixture.edge_bits.copy()
        self.board = board
        root, entry = fixture.route[0], fixture.route[-1]
        tips = sorted(set(board.neighbours[root]) - set(fixture.ports))
        children = sorted(fixture.tree.neighbors(entry))
        faces = {tuple(sorted(f)) for f in board.faces()}
        self.loops = []
        for anchor, neighbors in ((root, tips), (entry, children)):
            a, b = next(
                (a, b)
                for a, b in combinations(neighbors, 2)
                if tuple(sorted((anchor, a, b))) in faces
            )
            mask, cw = 0, 0
            for u, v in ((anchor, a), (a, b), (b, anchor)):
                edge = tuple(sorted((u, v)))
                bit = self.bits.setdefault(edge, 1 << len(self.bits))
                mask |= bit
                if v == edge[0]:
                    cw |= bit
            self.loops.append((mask, cw, (anchor, a, b)))
        seed_background = fixture.word[root]
        self.spur = None
        if spur:
            tip = self.loops[0][2][1]
            used = {p for e in self.bits for p in e}
            new = min(set(board.neighbours[tip]) - used)
            self.spur = tuple(sorted((tip, new)))
            bit = 1 << len(self.bits)
            self.bits[self.spur] = bit
            if (tip == self.spur[0]) == (fixture.sign == 1):
                seed_background |= bit
        seeds = []
        for s, r in product((0, 1), repeat=2):
            word = seed_background
            for parity, (mask, cw, _) in zip((s, r), self.loops):
                direction = cw if parity == 0 else cw ^ mask
                if fixture.sign == -1:
                    direction ^= mask
                word = (word & ~mask) | direction
            seeds.append((root, word))
        self.states = list(seeds)
        self.index = {x: i for i, x in enumerate(seeds)}
        hops, rings = set(), set()
        for i, state in enumerate(self.states):
            for moves, store in ((self.hops(state), hops), (self.flips(state), rings)):
                for other in moves:
                    if other not in self.index:
                        self.index[other] = len(self.states)
                        self.states.append(other)
                    j = self.index[other]
                    store.add(tuple(sorted((i, j))))
        self.hops_matrix = self.laplacian(hops, tree.guide.HOP)
        self.rings_matrix = self.laplacian(rings, 0.5)
        diag = [
            tree.guide.DEGREE_FIVE_COST * (6 - int(board.degree[p]))
            for p, _ in self.states
        ]
        self.hops_matrix += np.diag(diag)
        self.select = np.array(
            [
                fixture.boundary_head(w) - fixture.initial_inward == fixture.sign
                for _, w in self.states
            ]
        )
        np.testing.assert_array_equal(
            self.select, [p in fixture.receiver for p, _ in self.states]
        )
        self.indices = np.flatnonzero(self.select)
        # Formation words differ only by oriented loops, so their local Gauss reference agrees.
        address = sorted(set(p for edge in self.bits for p in edge))
        baseline = self.inward(seeds[0][1], address)
        for p, word in self.states:
            change = self.inward(word, address) - baseline
            expected = fixture.sign * (
                (np.array(address) == p).astype(int) - (np.array(address) == root)
            )
            np.testing.assert_array_equal(change, expected)

    def inward(self, word, addresses):
        heads = [a if word & bit else b for (a, b), bit in self.bits.items()]
        return np.array([heads.count(p) for p in addresses])

    def contact(self, state, edge):
        p, word = state
        a, b = edge
        bit = self.bits[edge]
        if p in edge and ((p == a) == bool(word & bit)) == (self.fixture.sign == 1):
            return (b if p == a else a), word ^ bit
        return state

    def hops(self, state):
        for e in self.bits:
            target = self.contact(state, e)
            if target != state:
                yield target

    def flips(self, state):
        p, w = state
        for mask, cw, _ in self.loops:
            if w & mask in (cw, cw ^ mask):
                yield p, w ^ mask

    def laplacian(self, edges, weight):
        h = np.zeros((len(self.states),) * 2)
        for i, j in edges:
            h[i, i] += weight
            h[j, j] += weight
            h[i, j] -= weight
            h[j, i] -= weight
        return h

    def rotation(self, perm):
        seq = self.fixture.pulse_schedule(perm)
        final = []
        for i in self.indices:
            state = self.states[i]
            for edge in seq:
                state = self.contact(state, edge)
            final.append(self.index[state])
        np.testing.assert_array_equal(sorted(final), self.indices)
        dest = np.searchsorted(self.indices, final)
        ports = [self.states[i][0] - self.fixture.base for i in self.indices]
        position_errors = sum(ports[dst] != perm[src] for dst, src in zip(dest, ports))
        square_errors = np.count_nonzero(dest[dest] != np.arange(len(dest)))
        return dest, int(position_errors), int(square_errors)


def receiver_channel(model, waves):
    """Trace the unobserved exterior after the flux click, retaining port and every receiver head."""
    local_mask = sum(
        bit
        for (a, b), bit in model.bits.items()
        if a in model.fixture.receiver and b in model.fixture.receiver
    )
    states = [model.states[i] for i in model.indices]
    outside = sorted({w & ~local_mask for _, w in states})
    inside = sorted({(p, w & local_mask) for p, w in states})
    tensor = np.zeros((len(outside), len(inside), 2), complex)
    for i, (p, w) in zip(model.indices, states):
        tensor[outside.index(w & ~local_mask), inside.index((p, w & local_mask))] = (
            waves[i]
        )
    amplitudes = tensor.reshape(len(outside), -1).T
    choi = amplitudes @ amplitudes.conj().T / 2
    pt = (
        choi.reshape(len(inside), 2, len(inside), 2)
        .transpose(0, 3, 2, 1)
        .reshape(choi.shape)
    )
    cross = (
        (tensor[:, :, 0] + tensor[:, :, 1]).T
        @ (tensor[:, :, 0] - tensor[:, :, 1]).conj()
        / 2
    )
    return {
        "herald_probability_mixed_input": float(np.trace(choi).real),
        "minimum_reference_partial_transpose_eigenvalue": float(
            np.linalg.eigvalsh(pt)[0]
        ),
        "receiver_even_odd_cross_operator_norm": float(np.linalg.norm(cross)),
        "exterior_words": len(outside),
        "receiver_records": len(inside),
    }


def source_mirror_certificate(model):
    if model.spur:
        return None
    _, a, b = model.loops[0][2]
    turn = {a: b, b: a}
    perm = []
    exterior_action = {}
    local_mask = sum(
        bit
        for (u, v), bit in model.bits.items()
        if u in model.fixture.receiver and v in model.fixture.receiver
    )
    for k, (p, word) in enumerate(model.states):
        transformed = 0
        for (u, v), bit in model.bits.items():
            x, y = turn.get(u, u), turn.get(v, v)
            if (x < y) == bool(word & bit):
                transformed |= model.bits[tuple(sorted((x, y)))]
        perm.append(model.index[turn.get(p, p), transformed])
        if model.select[k]:
            if (transformed & local_mask) != (word & local_mask) or turn.get(p, p) != p:
                raise ArithmeticError("Source mirror acts on the received packet")
            old, new = word & ~local_mask, transformed & ~local_mask
            if exterior_action.setdefault(old, new) != new:
                raise ArithmeticError("Exterior mirror depends on the receiver")
    for h in (model.hops_matrix, model.rings_matrix):
        np.testing.assert_array_equal(h[np.ix_(perm, perm)], h)
    np.testing.assert_array_equal(np.array(perm)[perm], np.arange(len(perm)))
    np.testing.assert_array_equal(perm[:4], [2, 3, 0, 1])
    return {
        "tips": [a, b],
        "exact_commutation": True,
        "exterior_word_permutation": sorted(exterior_action.items()),
    }


def closed_address_support(ports, first, second, a, b):
    """Largest support with correct position turns after any sequence of the two commands."""
    good = (ports[first] == np.asarray(a)[ports]) & (
        ports[second] == np.asarray(b)[ports]
    )
    while True:
        kept = good & good[first] & good[second]
        if np.array_equal(kept, good):
            return kept
        good = kept


def run_sign(board, route, sign, spur=False):
    fixture = tree.WrittenTree(board, route, sign)
    model = Circulation(fixture, board, spur)
    half = [
        p
        for p in sorted(tree.rotations())
        if any(p[i] != i for i in range(12)) and all(p[p[i]] == i for i in range(12))
    ]
    ops = [model.rotation(p) for p in half]
    pairs = [
        (i, j)
        for i, j in combinations(range(len(half)), 2)
        if all(half[i][half[j][k]] == half[j][half[i][k]] for k in range(12))
    ]
    if len(pairs) != 15:
        raise ValueError("Incomplete perpendicular-axis census")
    ports = np.array([model.states[k][0] - fixture.base for k in model.indices])
    ranks = np.bincount(ports, minlength=12)
    masks = [
        closed_address_support(ports, ops[i][0], ops[j][0], half[i], half[j])
        for i, j in pairs
    ]
    order_disagreement = [
        int(np.count_nonzero(ops[j][0][ops[i][0]][mask] != ops[i][0][ops[j][0]][mask]))
        for (i, j), mask in zip(pairs, masks)
    ]
    panels = []
    for rate, parity in product((1, 0), (1, -1)):
        h = model.hops_matrix + rate * model.rings_matrix
        values, vec = np.linalg.eigh(h)
        sharp_heads = np.zeros((len(values), 2))
        sharp_heads[:2, 0] = 1 / np.sqrt(2)
        sharp_heads[2:4, 1] = 1 / np.sqrt(2)
        columns = vec.T @ sharp_heads
        prep = np.zeros(len(values))
        prep[:4] = [0.5, 0.5, parity / 2, parity / 2]
        coeffs = vec.T @ prep
        readings = []
        for time in tree.TIMES:
            wave = vec @ (np.exp(-1j * time * values) * coeffs)
            np.testing.assert_allclose(np.vdot(wave, wave), 1.0, atol=1e-13)
            local = wave[model.indices]
            arrival = float(np.vdot(local, local).real)
            outcomes = []
            for i, j in pairs:
                xy, yx = ops[j][0][ops[i][0]], ops[i][0][ops[j][0]]
                first, second = local[np.argsort(xy)], local[np.argsort(yx)]
                constructive = float(np.linalg.norm((first + second) / 2) ** 2)
                destructive = float(np.linalg.norm((first - second) / 2) ** 2)
                np.testing.assert_allclose(
                    constructive + destructive, arrival, atol=2e-13
                )
                outcomes.append([constructive, destructive, max(0.0, 1 - arrival)])
            channel = receiver_channel(
                model, vec @ (np.exp(-1j * time * values)[:, None] * columns)
            )
            readings.append(
                {
                    "time": time,
                    "arrival": arrival,
                    "per_pair": outcomes,
                    "source_qubit_to_receiver": channel,
                    "geometrically_valid_weight": [
                        float(np.vdot(local[mask], local[mask]).real) for mask in masks
                    ],
                }
            )
        panels.append(
            {
                "ring_rate": 0.5 * rate,
                "source_circulation_phase": parity,
                "readings": readings,
            }
        )
    return {
        "sign": sign,
        "source_spur": model.spur,
        "source_mirror": source_mirror_certificate(model),
        "dimension": len(model.states),
        "heralded_sharp_states": len(model.indices),
        "triangles": [list(face) for _, _, face in model.loops],
        "axis_pairs": [list(p) for p in pairs],
        "half_turn_calibration": [
            {
                "permutation": list(p),
                "wrong_address": o[1],
                "wrong_word_after_square": o[2],
            }
            for p, o in zip(half, ops)
        ],
        "rank_by_local_address": ranks.tolist(),
        "rank_preserving_half_turns": sum(
            np.array_equal(ranks, ranks[list(p)]) for p in half
        ),
        "closed_address_subspace_sizes": [int(mask.sum()) for mask in masks],
        "order_disagreements_on_closed_address_subspace": order_disagreement,
        "calibrated_axis_pairs": sum(
            ops[i][1:] == (0, 0) and ops[j][1:] == (0, 0) for i, j in pairs
        ),
        "outcomes": ["constructive", "destructive", "no signed flux"],
        "panels": panels,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, required=True)
    args = cli.parse_args()
    board, route = tree.guide.longest_unpinned_route()
    results = [
        run_sign(board, route, s, spur) for spur in (False, True) for s in (-1, 1)
    ]
    with args.output.open("x") as out:
        json.dump({"sectors": results}, out, indent=2, allow_nan=False)
        out.write("\n")
