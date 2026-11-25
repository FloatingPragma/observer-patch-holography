"""Native preparation, free channel and local process tomography on two owner loops."""

from collections import defaultdict
from itertools import product
from pathlib import Path
import json

import networkx as nx
import numpy as np
from scipy.linalg import expm
from sympy import Matrix, Poly

from ..live_trimer.measure import LAW, blocks, geometry
from ..winding_memory import measure as winding

HERE = Path(__file__).resolve().parent
EDGES = ((0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3))
SHOTS = 2**24
K = int(winding.flight.candidate.KAPPA)
TIMES = (0, 1, 2, 4, 8, 12, 16)
PAULI = [
    np.array(x, complex)
    for x in ([[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]])
]


def stock(odd):
    """Write either contrast-sign parity from a balanced physical vacancy."""
    native = sum(k * t for k, t in zip((1, K), LAW.local_term_parts(), strict=True))
    rails = [LAW.SRC.index(p) for p in ((0, 1), (1, 0))]
    base = expm(-0.5j * np.pi * native[np.ix_(rails, rails)])[:, 0]
    # The vacancy's material private reference gives +8 relative phase rate.
    dwell = (np.pi / 2 if odd else 3 * np.pi / 2) / K
    base[1] *= np.exp(1j * K * dwell)
    state = LAW.echo()[:, rails] @ base
    columns = np.zeros((99, 2))
    for c, site in enumerate(rails):
        for d in (-1, 1):
            columns[9 + 9 * LAW.LABELS.index((c, d)) + site, c] = (
                d if odd else 1
            ) / np.sqrt(2)
    amplitudes = columns.T @ state
    np.testing.assert_allclose(columns @ amplitudes, state, atol=1e-12)
    np.testing.assert_allclose(abs(amplitudes) ** 2, 0.5, atol=1e-12)
    full = expm(-1j * np.pi / (K + odd) * native)
    np.testing.assert_allclose(abs(columns[:, 0] @ full @ columns[:, 1]), 1, atol=1e-12)


def fixture():
    board, ports = geometry()
    edges = [(ports[a], ports[b]) for a, b in EDGES]
    active = {tuple(sorted(e)) for e in edges}
    side = [
        min(
            p
            for p in board.neighbours[ports[i]]
            if p not in ports and p // 12 == ports[i] // 12 and board.degree[p] == 5
        )
        for i in (1, 4)
    ]
    cuts = {e for e in board.seams if set(ports).intersection(e)} - active
    staging = next(p for p in board.neighbours[ports[0]] if p not in ports + side)
    shuttle = tuple(sorted((staging, ports[0])))
    assert shuttle in cuts
    forbidden = set(ports + side + [staging])
    donor, mate = next(
        e
        for e in board.seams
        if all(
            p not in forbidden and not forbidden.intersection(board.neighbours[p])
            for p in e
        )
    )
    store_cuts = {
        e
        for e in board.seams
        if set(e).intersection((donor, mate)) and e != (donor, mate)
    }
    bypass = nx.Graph(board.seams)
    bypass.remove_edges_from(list(store_cuts) + edges)
    quiet = bypass.copy()
    quiet.remove_edges_from(cuts - {shuttle})
    storage = next(
        p
        for p in sorted(quiet)
        if p not in forbidden | {donor, mate}
        and not any(p in e for e in cuts | store_cuts)
        and nx.has_path(quiet, p, staging)
    )
    loads = np.ones(board.n, int)
    loads[storage], loads[donor] = 0, 2

    def route(graph, a, b):
        path = nx.shortest_path(graph, a, b)
        for x, y in winding.flight.closed_walk(path):
            loads[x], loads[y] = loads[y], loads[x]

    # Isolate the compensator. The blank dimer parks the excess on either end
    # while a shutter is written at unit load on the other end.
    for e in store_cuts:
        if donor in e:
            loads[mate], loads[donor] = loads[donor], loads[mate]
        assert all(loads[v] == 1 for v in e)
        if donor in e:
            loads[mate], loads[donor] = loads[donor], loads[mate]
    for port in ports:
        if board.degree[port] == 5:
            route(bypass, storage, port)
            assert loads[port] == 0
            route(bypass, storage, port)
    for a, b in edges[:6]:
        route(bypass, storage, a)
        assert loads[a] == 0 and loads[b] == 1
        # Phase the second coordinate with the first one isolated. Neither
        # herald consumes a vacancy, and the corrected winding coin is zero.
        isolated = bypass.copy()
        isolated.remove_node(a)
        assert any(nx.has_path(isolated, b, ref) for ref in side)
        route(bypass, storage, a)
    source_lane = tuple(sorted((ports[1], side[0])))
    receiver_lane = tuple(sorted((ports[4], side[1])))
    assert source_lane in cuts and receiver_lane in cuts
    for e in cuts - {shuttle, source_lane}:
        assert all(loads[p] == 1 for p in e)
    route(quiet, storage, staging)
    loads[staging], loads[ports[0]] = loads[ports[0]], loads[staging]
    # Wind the first marked bond while writing the final shutter, then undo it.
    loads[ports[0]], loads[ports[1]] = loads[ports[1]], loads[ports[0]]
    assert all(loads[p] == 1 for p in shuttle)
    loads[ports[0]], loads[ports[1]] = loads[ports[1]], loads[ports[0]]
    assert loads[ports[0]] == 0 and loads[donor] == 2
    # The emitter keeps its reference lane blank through preparation and
    # closes it only after both logical words put the vacancy on the gateway.
    assert all(loads[p] == 1 for p in source_lane)
    # Once loaded, all external ports outside the isolated compensating
    # dimer stay filled. The receiver herald puts both lane ends at one.
    sealed = nx.Graph(board.seams)
    sealed.remove_edges_from(cuts | store_cuts)
    assert set(nx.node_connected_component(sealed, donor)) == {donor, mate}
    assert set(nx.node_connected_component(sealed, ports[0])) == set(ports)
    assert loads[side[1]] == 1 and loads[side[0]] == 1
    winding.flight.seam_data(K)
    for odd in (False, True):
        stock(odd)
    return dict(
        ports=ports,
        outer_phase_references=side,
        phase_lanes=[source_lane, receiver_lane],
        outer_shutters=len(cuts),
        final_shutter=shuttle,
        stock_vacancy=storage,
        compensating_dimer=[donor, mate],
        compensator_shutters=len(store_cuts),
        erased_private_ports=[p for p in ports if board.degree[p] == 5],
    )


def build(odd):
    """Follow every nonzero integer local column, retaining the six direction coins."""
    words = [(0, 0)]
    indices = {words[0]: 0}
    columns = []
    for vacancy, heads in words:
        outgoing = [defaultdict(int) for _ in EDGES]
        for e, (left, right) in enumerate(EDGES):
            if vacancy not in (left, right):
                continue
            marked = e < 6
            coin = (heads >> e) & 1 if marked else 0
            h = (
                blocks()[1, -1 if odd and e == 0 else 1]
                if marked
                else blocks()["blank"]
            )
            site = LAW.SRC.index(tuple(int(p != vacancy) for p in (left, right)))
            source = 9 * coin + site
            for target in np.flatnonzero(h[:, source]):
                loads = LAW.SRC[target % 9]
                position = (left, right)[loads.index(0)]
                mask = heads ^ ((coin ^ (target // 9)) << e) if marked else heads
                outgoing[e][position, mask] += int(h[target, source])
            filled = 9 * coin + LAW.SRC.index((1, 1))
            outgoing[e][vacancy, heads] -= int(h[filled, filled])
        columns.append(outgoing)
        for col in outgoing:
            for word, value in col.items():
                if value and word not in indices:
                    indices[word] = len(words)
                    words.append(word)
    terms = np.zeros((7, len(words), len(words)), int)
    for j, columns_at_j in enumerate(columns):
        for e, col in enumerate(columns_at_j):
            for word, value in col.items():
                if value:
                    terms[e, indices[word], j] = value
    np.testing.assert_array_equal(terms, terms.transpose(0, 2, 1))
    assert len(words) == 16
    return words, indices, terms


def source_reflection(words, ix):
    def flip(word):
        site, pattern = word
        left = pattern & 7
        backwards = ((left & 1) << 2) | (left & 2) | ((left & 4) >> 2)
        return ({1: 2, 2: 1}.get(site, site), (pattern & 56) | (7 ^ backwards))

    return [ix[flip(w)] for w in words]


def unitary_gate(terms, edge, duration):
    return expm(-1j * duration * terms[edge] / 8)


def phase_lane(relative):
    """Native shutter, relay to the private reference, dwell and return."""
    echo = LAW.echo()
    blank_filled = np.eye(99)[:, LAW.SRC.index((1, 1))]
    silent_filled = echo @ blank_filled
    np.testing.assert_allclose(echo @ silent_filled, blank_filled, atol=1e-13)
    for term in LAW.local_term_parts():
        np.testing.assert_allclose(term @ silent_filled, 0, atol=1e-13)
    blank, private = winding.flight.seam_data(K)
    relay = expm(-1j * np.pi * blank)
    rate = K * (private[3] - private[0])
    dwell = (relative % (2 * np.pi)) / rate
    # Site order is (tip, remote reference). Its ordinary private term
    # addresses the second site's load against its private zero rail.
    phase = np.diag([np.exp(-1j * K * private[3 * b] * dwell) for a, b in LAW.SRC])
    contact = relay.conj().T @ phase @ relay
    indices = [LAW.SRC.index((n, 1)) for n in (0, 1)]
    filled_phase = contact[indices[1], indices[1]]
    expected = np.eye(9, dtype=complex)[:, indices] @ np.diag(
        [np.exp(1j * relative), 1]
    )
    np.testing.assert_allclose(contact[:, indices] / filled_phase, expected, atol=1e-12)
    return contact[indices[0], indices[0]] / filled_phase


def face(words, ix, terms, odd, owner, axis, sign):
    """Unwind, route one logical coordinate to a private phase and rotate the local bond."""
    base = 0 if owner == 0 else 3
    first = base
    angle = {0: np.pi / 4, 1: np.pi / 4, 2: 0 if sign == 1 else np.pi / 2}[axis]
    # An addressed active contact prepares [cos(angle), i sin(angle)] up to phase.
    relative = {0: -sign * np.pi / 2, 1: 0 if sign == 1 else np.pi, 2: 0}[axis]
    partial = unitary_gate(terms, first, angle / (4.5 if first == 0 and odd else 4))
    phase = phase_lane(relative)
    routed_phase = np.diag([phase if site == base + 1 else 1 for site, _ in words])
    # A full native swap lends the charged tip to a blank degree-five reference.
    full = [unitary_gate(terms, e, np.pi / 8) for e in (first + 1, first + 2)]
    return full[1] @ full[0] @ routed_phase @ partial


def kraus(u, ix):
    source = [ix[0, 0], ix[0, 7]]
    return [u[np.ix_([ix[3, a + 8 * b] for b in (0, 7)], source)] for a in (0, 7)]


def choi(operators):
    return sum(np.outer(k.T.ravel(), k.T.ravel().conj()) for k in operators) / 2


def pt(matrix):
    return matrix.reshape(2, 2, 2, 2).transpose(2, 1, 0, 3).reshape(4, 4)


def tomography(words, ix, terms, odd, u, rng):
    accepted = [ix[3, a + 8 * b] for a, b in product((0, 7), repeat=2)]
    incoming = [ix[0, 0], ix[0, 7]]
    outputs = []
    records = []
    for axis in range(3):
        for sign in (1, -1):
            prep = face(words, ix, terms, odd, 0, axis, sign)[:, incoming[0]]
            code = prep[incoming]
            np.testing.assert_allclose(PAULI[axis] @ code, sign * code, atol=1e-12)
            assert abs(np.linalg.norm(code) - 1) < 1e-12
            # Close the source tip's blank reference contact by its unit-pair
            # echo after the completed circuit has returned BOTH arms to port 0.
            assert all(words[i][0] == 0 for i in np.flatnonzero(abs(prep) > 1e-12))
            y = u @ prep
            gated = np.zeros_like(y)
            # Nondestructive material herald only. Each accepted word has a
            # filled receiving tip, and the sealed dimer keeps the outer
            # reference filled. The inverse unit-pair echo opens this lane
            # before any decoding pulse is applied.
            assert accepted == [ix[3, a + 8 * b] for a, b in product((0, 7), repeat=2)]
            assert set(accepted) == {i for i, w in enumerate(words) if w[0] == 3}
            gated[accepted] = y[accepted]
            total = float(np.vdot(gated, gated).real)
            means = []
            bank = []
            for read in range(3):
                # Undo the receiving preparation. Its inverse maps the positive
                # state to the gateway and the negative one to the first local tip.
                decoder = face(words, ix, terms, odd, 1, read, 1).conj().T
                z = decoder @ gated
                plus = float(
                    sum(abs(z[j]) ** 2 for j, w in enumerate(words) if w[0] == 3)
                )
                # The other logical output is at the first tip, not another head oracle.
                minus = float(
                    sum(abs(z[j]) ** 2 for j, w in enumerate(words) if w[0] == 4)
                )
                np.testing.assert_allclose(plus + minus, total, atol=1e-12)
                p = np.maximum([plus, minus, 1 - total], 0)
                p /= sum(p)
                counts = rng.multinomial(SHOTS, p)
                bank.append(counts.tolist())
                means.append((counts[0] - counts[1]) / SHOTS)
            trace = sum(c[0] + c[1] for c in bank) / (3 * SHOTS)
            outputs.append(
                (trace * np.eye(2) + sum(y * s for y, s in zip(means, PAULI))) / 2
            )
            records.append(bank)
    plus_z, minus_z = outputs[4:]
    cross = (outputs[0] - outputs[1] + 1j * (outputs[2] - outputs[3])) / 2
    j = np.block([[plus_z, cross], [cross.conj().T, minus_z]]) / 2
    # Bonferroni-Hoeffding for 18 independent settings. Each parity is in
    # [-1,1], each success indicator in [0,1]. Bound operator error by
    # the Frobenius norm of the four reconstructed 2x2 blocks after transpose.
    eps = np.sqrt(2 * np.log(2 * 36 / 0.01) / SHOTS)
    local = np.sqrt(((eps / 2) ** 2 + 3 * eps**2) / 2)
    error = np.sqrt(2 * local**2 + 2 * (2 * local) ** 2) / 2
    return dict(
        starts_per_basis=SHOTS,
        counts=records,
        measured_min=float(np.linalg.eigvalsh(pt(j))[0]),
        simultaneous99_radius=float(error),
    )


def evaluate(odd):
    words, ix, terms = build(odd)
    h8 = sum(terms)
    flip = source_reflection(words, ix)
    invariant = np.array_equal(h8[np.ix_(flip, flip)], h8)
    assert invariant == (not odd)
    energies = Poly(Matrix(h8.tolist()).charpoly().as_expr())
    code = np.eye(len(words), dtype=int)[
        :, [ix[3, a + 8 * b] for a, b in product((0, 7), repeat=2)]
    ]
    sub = code.T @ h8 @ code
    np.testing.assert_array_equal(sub, 36 * np.eye(4))
    escape = h8 @ code - code @ sub
    np.testing.assert_array_equal(escape.T @ escape, 1040 * np.eye(4))
    dwell = np.linalg.svd(code.T @ expm(-0.25j * h8 / 8) @ code, compute_uv=False) ** 2
    rows = []
    for t in TIMES:
        u = expm(-1j * t * h8 / 8)
        j = choi(kraus(u, ix))
        rows.append(
            dict(
                time=t,
                success=float(j.trace().real),
                min_pt=float(np.linalg.eigvalsh(pt(j))[0]),
            )
        )
    u = expm(-2j * h8)
    operators = kraus(u, ix)
    success = np.linalg.eigvalsh(sum(a.conj().T @ a for a in operators))
    return dict(
        input_success_range=success.tolist(),
        reflection=invariant,
        dimension=len(words),
        simple_spectrum=energies.gcd(energies.diff()).degree() == 0,
        conditional_mean=4.5,
        conditional_variance=16.25,
        dwell_quarter=[float(dwell[-1]), float(dwell[0])],
        records=rows,
        tomography=tomography(
            words, ix, terms, odd, u, np.random.default_rng(13079 + odd)
        ),
    )


def run():
    from .certificate import witness

    return dict(
        apparatus=fixture(),
        even=evaluate(False),
        odd=evaluate(True),
        exact_witness=witness(),
    )


if __name__ == "__main__":
    (HERE / "results.json").write_text(json.dumps(run(), indent=2) + "\n")
