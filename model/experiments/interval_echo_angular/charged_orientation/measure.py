"""Measure the signed-defect orbital channel, controlled rotations and two-packet exchange."""

import argparse
import hashlib
import json
import math
from collections import deque
from itertools import combinations, product
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.electron_field.dressed_exchange import (
    search as written,
)
from model.executable_model import cell, federation
from model.executable_model.symmetry import rotations

HERE = Path(__file__).resolve().parent
TIMES = (0, 2, 8, 32, 128)
SHOTS = 2**42
SEED = 841117
EPS = math.sqrt(math.log(2 * 100000 / 0.01) / (2 * SHOTS))
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
AXES = [np.linalg.eigh(p)[1][:, ::-1] for p in PAULI]


def path(adjacency, start, end, excluded=()):
    parents = {start: None}
    queue = deque([start])
    forbidden = set(excluded)
    while queue:
        at = queue.popleft()
        if at == end:
            route = [end]
            while parents[route[-1]] is not None:
                route.append(parents[route[-1]])
            return route[::-1]
        for other in adjacency[at]:
            if other not in forbidden and other not in parents:
                parents[other] = at
                queue.append(other)
    return None


def design(board):
    eligible = {}
    for owner in range(board.cells):
        gates = set(
            np.flatnonzero(board.degree[12 * owner : 12 * owner + 12] == 6).tolist()
        )
        start = min(gates)
        turns = [
            tuple(r)
            for r in rotations()
            if {r[p] for p in gates} == gates
            and r[start] != start
            and r[r[r[start]]] == start
        ]
        if turns:
            turn = min(turns)
            eligible[owner] = (
                turn,
                [12 * owner + j for j in (start, turn[start], turn[turn[start]])],
            )
    bridges = [e for e in board.external if all(p // 12 in eligible for p in e)]
    a, b = min(bridges)
    source, first = a // 12, b // 12
    neighbours = [set() for _ in range(board.cells)]
    for a, b in board.external:
        neighbours[a // 12].add(b // 12)
        neighbours[b // 12].add(a // 12)
    second = sorted(
        set.union(*(neighbours[i] for i in neighbours[source]))
        .difference(neighbours[source], {source})
        .intersection(eligible)
    )
    return source, [first] + second[:1], eligible


def code_columns(gateways, n):
    result = np.zeros((n, 2), complex)
    wave = np.exp(2j * np.pi * np.arange(3) / 3) / np.sqrt(3)
    result[gateways] = np.column_stack([wave, wave.conj()])
    return result


def generator(board, sign, strength):
    """Subtract the uniform material energy and a 6κ shift from the exact single-defect block."""
    idx = [law.SRC.index(w) for w in ((1 + sign, 1), (1, 1 + sign))]
    blank = law.seam_term(strength)[np.ix_(idx, idx)]
    np.testing.assert_allclose(blank, [[strength + 0.5, -0.5], [-0.5, strength + 0.5]])
    return board.laplacian() / 2 + sparse.diags(
        2 * sign * strength * (6 - board.degree).astype(float)
    )


def step(h, states, t):
    if t == 0:
        return states.copy()
    return expm_multiply(-1j * t * h, states, traceA=-1j * t * h.diagonal().sum())


def counts(probabilities, rng):
    p = np.asarray(probabilities, dtype=float)
    if min(p) < -1e-12 or abs(p.sum() - 1) > 2e-12:
        raise ArithmeticError("Detector lost an original-start outcome")
    p = np.maximum(p, 0)
    return rng.multinomial(SHOTS, p / p.sum()).tolist()


def tomography(local, outgoing, rng):
    transfer = outgoing.conj().T @ local
    groups, states = [], []
    for prep in AXES:
        for v in prep.T:
            received = transfer @ v
            accepted = np.linalg.norm(received) ** 2
            present = np.linalg.norm(local @ v) ** 2
            panel = [
                counts(
                    [
                        *(abs(analyzer.conj().T @ received) ** 2),
                        present - accepted,
                        1 - present,
                    ],
                    rng,
                )
                for analyzer in AXES
            ]
            groups.append(panel)
            f = np.array(panel) / SHOTS
            spin = f[:, 0] - f[:, 1]
            states.append(
                (f[2, :2].sum() * np.eye(2) + np.einsum("i,ijk->jk", spin, PAULI)) / 2
            )
    x, y, z = np.asarray(states).reshape(3, 2, 2, 2)
    cross = ((x[0] - x[1]) + 1j * (y[0] - y[1])) / 2
    choi = np.block([[z[0], cross], [cross.conj().T, z[1]]]) / 2
    pt = choi.reshape(2, 2, 2, 2).swapaxes(1, 3).reshape(4, 4)
    measured = float(-np.linalg.eigvalsh(pt)[0])
    return dict(
        counts=groups,
        measured_negative_pt=measured,
        resolved=bool(measured > 5 * EPS),
        exact_negative_pt=float(abs(np.linalg.det(transfer)) / 2),
        accepted=float(np.trace(transfer.conj().T @ transfer).real / 2),
        present=float(np.trace(local.conj().T @ local).real / 2),
    )


def compile_turn(perm):
    """Factor a proper cell permutation into adjacent port swaps on the protected cell."""
    edges = []
    done = set()
    for first in range(len(perm)):
        if first in done:
            continue
        cycle = []
        at = first
        while at not in done:
            cycle.append(at)
            done.add(at)
            at = perm[at]
        for last in cycle[1:]:
            route = path(cell().neighbours, first, last)
            forward = list(zip(route, route[1:]))
            edges += forward + forward[-2::-1]
    moved = list(range(12))
    for a, b in edges:
        moved[a], moved[b] = moved[b], moved[a]
    if moved != np.argsort(perm).tolist():
        raise ArithmeticError("Template rotation was not implemented")
    return edges


def rotated(psi, edges):
    out = psi.copy()
    for a, b in edges:
        out[[a, b]] = out[[b, a]]
    return out


def half_turns():
    ident = tuple(range(12))
    order_two = [
        r
        for r in sorted(rotations())
        if r != ident and tuple(r[r[k]] for k in ident) == ident
    ]
    for a, b in combinations(order_two, 2):
        if all(a[b[i]] == b[a[i]] for i in ident):
            return compile_turn(a), compile_turn(b)
    raise ValueError("No perpendicular half-turns in the template")


def rotation_read(local, turn, halves, rng):
    ab, ba = halves[0] + halves[1], halves[1] + halves[0]
    words = {
        "idle": ([], []),
        "120": ([], turn),
        "360": ([], turn * 3),
        "720": ([], turn * 6),
        "xy_vs_yx": (ab, ba),
    }
    c3 = np.diag(np.exp(1j * np.array([-1.0, 1.0]) * np.pi / 3))
    sx, sy = -1j * PAULI[0], -1j * PAULI[1]
    one = np.eye(2)
    spin_arms = {
        "idle": (one, one),
        "120": (one, c3),
        "360": (one, c3 @ c3 @ c3),
        "720": (one, np.linalg.matrix_power(c3, 6)),
        "xy_vs_yx": (sy @ sx, sx @ sy),
    }
    results = {}
    for name, (first, second) in words.items():
        A, B = rotated(local, first), rotated(local, second)
        panel = []
        for side in range(2):
            bright = np.linalg.norm((A[:, side] + B[:, side]) / 2) ** 2
            dark = np.linalg.norm((A[:, side] - B[:, side]) / 2) ** 2
            panel.append(counts([bright, dark, 1 - bright - dark], rng))
        results[name] = {"orbital": panel}
        spin_first, spin_second = spin_arms[name]
        spin_panel = []
        for side in range(2):
            A_spin = np.outer(A[:, side], spin_first[:, 0])
            B_spin = np.outer(B[:, side], spin_second[:, 0])
            bright = np.linalg.norm((A_spin + B_spin) / 2) ** 2
            dark = np.linalg.norm((A_spin - B_spin) / 2) ** 2
            spin_panel.append(counts([bright, dark, 1 - bright - dark], rng))
        results[name]["imposed_spin_half"] = spin_panel
    return results


def exchange_word(board, src, dst):
    """Route the six possible defect locations via a separate unoccupied cell."""
    parked = []
    opening = []
    for ports in (src, dst):
        owner = ports[0] // 12
        available = set(range(12 * owner, 12 * (owner + 1))) - set(ports)
        support = set(ports)
        for port in ports:
            nearby = sorted(set(board.neighbours[port]) & available)
            if not nearby:
                raise ValueError("Missing empty home side port for the exchange")
            side = nearby[0]
            available.remove(side)
            if side in support:
                raise ArithmeticError("The parking operation needs a vacant port")
            opening.append((port, side))
            parked.append(side)
            support.remove(port)
            support.add(side)
    goals = parked[3:] + parked[:3]
    for owner in range(board.cells):
        if owner in (src[0] // 12, dst[0] // 12):
            continue
        holding = [
            j for j in range(12 * owner, 12 * (owner + 1)) if board.degree[j] == 5
        ][:6]
        positions = parked.copy()
        steps = []
        complete = True
        for dests in (holding, goals):
            for i, dest in enumerate(dests):
                route = path(
                    board.neighbours,
                    positions[i],
                    dest,
                    set(positions) - {positions[i]},
                )
                if route is None:
                    complete = False
                    break
                steps.extend(zip(route, route[1:]))
                positions[i] = dest
            if not complete:
                break
        if complete:
            return opening + steps + opening[::-1]
    raise ValueError("Six supported orbitals cannot be separated during exchange")


def apply_pair_word(src, dst, word):
    supported = set(src + dst)
    for x, y in word:
        if x in supported and y in supported:
            raise ArithmeticError(
                "A programmed swap touches two possible defect addresses"
            )
        if x in supported or y in supported:
            supported.symmetric_difference_update((x, y))
    if supported != set(src + dst):
        raise ArithmeticError("The exchange did not restore the supported addresses")
    labels = list(product(src, dst))
    native, odd = np.zeros((9, 9)), np.zeros((9, 9))
    for col, (a, b) in enumerate(labels):
        fermion = 1
        for x, y in word:
            before = (a, b)
            if {a, b} == {x, y}:
                raise ValueError("Two occupied trajectories share a contact")
            transposition = {x: y, y: x}
            a, b = transposition.get(a, a), transposition.get(b, b)
            if (before[0] < before[1]) != (a < b):
                fermion *= -1
        if a != dst[src.index(labels[col][0])] or b != src[dst.index(labels[col][1])]:
            raise ArithmeticError("Tagged origins failed the intended exchange")
        row = labels.index((b, a))
        # Use the same cell-ordered basis at both ends, independent of port numbering.
        native[row, col], odd[row, col] = 1, fermion
    return native, odd


def reference_fidelities(channel, rng):
    """Each input mode is entangled with its own retained reference."""
    ident = np.eye(4)
    swap = ident[[0, 2, 1, 3]]
    local = abs(np.vdot(ident, channel) / 4) ** 2
    crossed = abs(np.vdot(swap, channel) / 4) ** 2
    # Two independently prepared ensembles realize the alternative Bell pairings.
    return dict(
        original_pairing=float(local),
        crossed_pairing=float(crossed),
        counts=[counts([x, 1 - x], rng) for x in (local, crossed)],
    )


def exchange_reads(native, odd, rng, label="registered"):
    states = {
        "same": np.array([1.0, 0, 0, 0]),
        "unlike": np.array([0.0, 1, 0, 0]),
        "symmetric": np.array([0.0, 1, 1, 0]) / np.sqrt(2),
        "antisymmetric": np.array([0.0, 1, -1, 0]) / np.sqrt(2),
    }
    out = {}
    for label, op in ((label, native), ("odd_worldlines", odd)):
        panel = {}
        for key, v in states.items():
            w = op @ v
            panel[key] = counts(
                [np.linalg.norm(v + w) ** 2 / 4, np.linalg.norm(v - w) ** 2 / 4], rng
            )
        out[label] = dict(parity=panel, references=reference_fidelities(op, rng))
    return out


def run():
    from scipy.linalg import expm

    swap = np.eye(9)[[law.SRC.index((b, a)) for a, b in law.SRC]]
    for strength in (candidate.KAPPA, 0):
        np.testing.assert_allclose(
            expm(-1j * np.pi * law.seam_term(strength)), swap, atol=3e-14
        )
    rng = np.random.default_rng(SEED)
    rx, ry = half_turns()
    native = []
    for level in range(3):
        board = federation(level)
        source, dests, eligible = design(board)
        src = eligible[source][1]
        initial = code_columns(src, board.n)
        dest = eligible[dests[0]][1]
        word = exchange_word(board, src, dest)
        perm, odd = apply_pair_word(src, dest, word)
        localcode = code_columns([0, 1, 2], 3)
        productcode = np.kron(localcode, localcode)
        pair = productcode.conj().T @ perm @ productcode
        signed = productcode.conj().T @ odd @ productcode
        np.testing.assert_allclose(pair, np.eye(4)[[0, 2, 1, 3]], atol=2e-14)
        np.testing.assert_allclose(signed, -pair, atol=2e-14)
        arms = []
        for sign in (-1, 1):
            for kappa in (candidate.KAPPA, 0):
                h = generator(board, sign, kappa)
                wave = initial.copy()
                previous = 0
                appointments = []
                for t in TIMES:
                    wave = step(h, wave, t - previous)
                    previous = t
                    np.testing.assert_allclose(
                        wave.conj().T @ wave, np.eye(2), atol=2e-11
                    )
                    readers = []
                    for owner in dests:
                        p, gates = eligible[owner]
                        aperture = wave[12 * owner : 12 * (owner + 1)]
                        code = code_columns([g % 12 for g in gates], 12)
                        data = tomography(aperture, code, rng)
                        data["turns"] = rotation_read(
                            aperture, compile_turn(p), (rx, ry), rng
                        )
                        readers.append(data)
                    appointments.append(dict(exposure=t, receivers=readers))
                arms.append(dict(sign=sign, kappa=kappa, appointments=appointments))
        native.append(
            dict(
                cells=board.cells,
                source=source,
                source_gateways=src,
                receivers=dests,
                receiver_gateways=[eligible[i][1] for i in dests],
                exchange_word=word,
                exchange={str(s): exchange_reads(pair, signed, rng) for s in (-1, 1)},
                arms=arms,
            )
        )
    extension = []
    for sign, head in product((-1, 1), (False, True)):
        fixture = written.Fixture(sign, head)
        result = written.six_pulse_exchange(fixture)["paired_buffer_eight_pulse_search"]
        field = np.array(result["native_doublet_return_matrix"])
        fermion = np.array(result["odd_permutation_doublet_return_matrix"])
        extension.append(
            dict(
                sign=sign,
                head_at_left=head,
                pulses=result["best_schedule"],
                readout=exchange_reads(field, fermion, rng, "written_head"),
            )
        )
    return dict(
        protocol_sha256=hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        shots=SHOTS,
        seed=SEED,
        frequency_radius=EPS,
        pt_radius=5 * EPS,
        native=native,
        written_extension=extension,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    output = parser.parse_args().output
    result = run()
    with output.open("x") as out:
        json.dump(result, out, separators=(",", ":"), allow_nan=False)
        out.write("\n")
