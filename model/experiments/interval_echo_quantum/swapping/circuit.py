"""Four occupied physical ports, with addressed gates on the 240-port federation."""

from collections import defaultdict
from functools import lru_cache
from itertools import pairwise, product

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import carrier
from model.experiments.interval_echo_consensus.empty_full import measure as source
from model.experiments.interval_echo_quantum import experiment as old

LAW = candidate.native_law
MENU = {
    "chsh": ((0, np.pi / 2), (np.pi / 4, -np.pi / 4)),
    "bilocal": ((np.pi / 4, -np.pi / 4), (np.pi / 4, -np.pi / 4)),
}
MODES = (
    "native",
    "right_total",
    "both_total",
    "right_which",
    "cut_alice",
    "separable_bob",
)


def apparatus():
    board = carrier.federation(0)
    owners = nx.Graph()
    owners.add_edges_from(sorted((p // 12, q // 12) for p, q in board.external))
    distances = nx.single_source_shortest_path_length(owners, 0)
    end = min(owner for owner, distance in distances.items() if distance == 4)
    line = min(nx.all_shortest_paths(owners, 0, end))
    alice, left, bob, right, charlie = line
    written = {tuple(12 * s + i for i in (0, 1)) for s in (left, right)}
    goals = [
        (left, (0, 5), alice, (0, 5)),
        (left, (1, 7), bob, (0, 5)),
        (right, (0, 5), bob, (1, 7)),
        (right, (1, 7), charlie, (0, 5)),
    ]
    routes = []
    for s, src, t, dst in goals:
        corridor = nx.Graph()
        corridor.add_edges_from(
            e
            for e in board.seams
            if e not in written and all(p // 12 in (s, t) for p in e)
        )
        for i, j in zip(src, dst, strict=True):
            path = min(nx.all_shortest_paths(corridor, 12 * s + i, 12 * t + j))
            route = list(pairwise(path))
            routes.append(route + route[-2::-1])
    return line, routes


@lru_cache(maxsize=None)
def pulse(time):
    ex, co = LAW.local_term_parts()
    h = ex + candidate.KAPPA * co
    np.testing.assert_allclose(h[LAW.N_SRC :, : LAW.N_SRC], 0, atol=1e-14)
    return expm(-1j * time * h[: LAW.N_SRC, : LAW.N_SRC])


def contact(state, edge, matrix):
    """Use the actual two-qutrit gate on each sparse occupation word."""
    out = defaultdict(complex)
    for word, amp in state.items():
        local = tuple(word.count(p) for p in edge)
        other = tuple(p for p in word if p not in edge)
        column = matrix[:, LAW.SRC.index(local)]
        for n in np.flatnonzero(abs(column) > 1e-13):
            a, b = LAW.SRC[n]
            key = tuple(sorted(other + (edge[0],) * a + (edge[1],) * b))
            out[key] += amp * column[n]
    return {w: a for w, a in out.items() if abs(a) > 1e-13}


def phase(state, p, angle):
    return {w: a * np.exp(1j * angle * w.count(p)) for w, a in state.items()}


def rails(owner):
    return [(12 * owner + a, 12 * owner + b) for a, b in ((0, 5), (1, 7))]


def prepare(owner, mode, phase_gate=phase):
    r, s = rails(owner)
    edge = (r[0], s[0])
    state = {tuple(sorted((r[1], s[1]))): 1}
    for rail in (r, s):
        state = phase_gate(contact(state, rail, pulse(np.pi / 2)), rail[0], -np.pi / 2)
    for _ in range(2):
        state = phase_gate(contact(state, edge, pulse(refocus_time())), r[0], np.pi)
    state = contact(state, s, pulse(np.pi / 2))
    arms = [state]
    if mode == "which":
        arms = [
            {w: a for w, a in state.items() if sum(w.count(p) for p in edge) == n}
            for n in (0, 1, 2)
        ]
    write = source.total_sensitive() if mode == "total" else LAW.echo()
    blocks = write[:, : LAW.N_SRC].reshape(-1, LAW.N_SRC, LAW.N_SRC)
    return [
        v
        for arm in arms
        for block in blocks
        if np.any(block)
        if (v := contact(arm, edge, block))
    ]


def refocus_time():
    return np.pi / (4 * (candidate.KAPPA + 0.5))


def distribute(line, routes, mode, phase_gate=phase):
    _, left, _, right, _ = line
    lmode = "total" if mode == "both_total" else "native"
    rmode = {"right_total": "total", "both_total": "total", "right_which": "which"}.get(
        mode, "native"
    )
    # All nine load words implement SWAP with the same phase.
    exact_swap = np.eye(LAW.N_SRC)[[LAW.SRC.index((b, a)) for a, b in LAW.SRC]]
    np.testing.assert_allclose(pulse(np.pi), exact_swap, atol=3e-14)
    perm = list(range(240))
    for i, steps in enumerate(routes):
        if mode == "cut_alice" and i < 2:
            continue
        for a, b in steps:
            perm[a], perm[b] = perm[b], perm[a]
    destination = np.argsort(perm).tolist()
    branches = []
    for left_state, right_state in product(
        prepare(left, lmode, phase_gate), prepare(right, rmode, phase_gate)
    ):
        branches.append(
            {
                tuple(sorted(destination[p] for p in a + b)): u * v
                for a, u in left_state.items()
                for b, v in right_state.items()
            }
        )
    return branches


def hadamard(state, rail, phase_gate):
    state = phase_gate(state, rail[0], -np.pi / 2)
    state = contact(state, rail, pulse(np.pi / 2))
    return phase_gate(state, rail[0], -np.pi / 2)


def bob_analyzer(state, bob, phase_gate=phase):
    left, right = rails(bob)
    state = hadamard(state, right, phase_gate)
    for _ in range(2):
        state = contact(state, (left[0], right[0]), pulse(refocus_time()))
        state = phase_gate(state, left[0], np.pi)
    for p in (left[0], right[0]):
        state = phase_gate(state, p, np.pi / 2)
    return hadamard(hadamard(state, right, phase_gate), left, phase_gate)


def read(word, rail):
    # The fresh local echo reports spectator minus active load. It acts on a
    # distinct memory. Its effect is checked on both signed and quiet inputs.
    delta = word.count(rail[1]) - word.count(rail[0])
    return {1: 0, -1: 1, 0: 2}[delta]


def central_record(branches, bob, mode, phase_gate=phase):
    output = defaultdict(list)
    for branch in branches:
        state = (
            branch if mode == "separable_bob" else bob_analyzer(branch, bob, phase_gate)
        )
        groups = defaultdict(dict)
        for w, a in state.items():
            b = tuple(read(w, rail) for rail in rails(bob))
            groups[b][w] = a
        for b, arm in groups.items():
            # A quiet central pointer has the fixed bit zero. Keep the arms
            # incoherent after recording each distinct raw pointer.
            bits = tuple(1 if pointer == 1 else 0 for pointer in b)
            output[bits].append(arm)
    return output


def outer(state, rail, angle, phase_gate):
    return contact(
        phase_gate(state, rail[0], -np.pi / 2), rail, pulse(angle % (2 * np.pi))
    )


def distributions(line, routes, mode, menu, phase_gate=phase):
    """Joint measured mass P(two central bits, terminal pointers | x,z)."""
    alice, _, bob, _, charlie = line
    central = central_record(
        distribute(line, routes, mode, phase_gate), bob, mode, phase_gate
    )
    ar, cr = rails(alice)[0], rails(charlie)[0]
    ax, cz = MENU[menu]
    prob = np.zeros((2, 2, 2, 2, 3, 3))
    for (b, d), arms in central.items():
        for x, z in product(range(2), repeat=2):
            for arm in arms:
                measured = outer(
                    outer(arm, ar, ax[x], phase_gate), cr, cz[z], phase_gate
                )
                for w, amp in measured.items():
                    prob[b, d, x, z, read(w, ar), read(w, cr)] += abs(amp) ** 2
    np.testing.assert_allclose(prob.sum(axis=(0, 1, 4, 5)), 1, atol=1e-13)
    return prob


def settings_probability(owner):
    # Distinct unused local seam chosen geometrically after the central record.
    rail = tuple(12 * owner + p for p in (2, 3))
    state = contact({(rail[1],): 1}, rail, pulse(np.pi / 2))
    return [
        sum(abs(a) ** 2 for w, a in state.items() if read(w, rail) == b) for b in (0, 1)
    ]


def check_readers():
    old.verify_read()
