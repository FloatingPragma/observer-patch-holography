"""Ten routed rails, sharp echo readers, and addressed vacancy-dimer gates."""

from itertools import combinations, pairwise

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.face_frame import (
    shared_report,
)
from model.executable_model import federation
from model.experiments.interval_echo_consensus.signed_effects import record_transfer


def dimer_gates():
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + candidate.KAPPA * contrast
    blank = np.eye(law.DIM)[:, :9]
    np.testing.assert_allclose(h @ blank, blank @ h[:9, :9], atol=2e-13)
    swap = np.eye(9)[:, [law.SRC.index((b, a)) for a, b in law.SRC]]
    np.testing.assert_allclose(expm(-1j * np.pi * h[:9, :9]), swap, atol=2e-13)
    # At a bridge with exchange suppressed, every integer contrast gives phase one.
    frozen = candidate.KAPPA * contrast
    np.testing.assert_allclose(frozen @ blank, blank @ frozen[:9, :9], atol=2e-13)
    np.testing.assert_allclose(
        expm(-1j * np.pi * frozen[:9, :9]), np.eye(9), atol=2e-13
    )
    indices = [law.SRC.index(pair) for pair in ((0, 1), (1, 0))]
    code = blank[:, indices]
    contact = code.T @ h @ code
    np.testing.assert_allclose(h @ code, code @ contact, atol=2e-13)
    # A unit of private contrast to its empty reference at a degree-five marked port.
    reference = np.diag([candidate.KAPPA, 0])
    flip = expm(-1j * np.pi * contact)
    phase = expm(-1j * np.pi / 8 * reference)
    change_basis = np.eye(2)
    for kind, _, _, dwell in record_transfer.basis_rotation(0, 1):
        term = contact if kind == "contact" else reference
        change_basis = expm(-1j * np.pi * float(dwell) * term) @ change_basis
    np.testing.assert_allclose(flip, [[0, 1], [1, 0]], atol=2e-13)
    np.testing.assert_allclose(-phase, np.diag([1, -1]), atol=2e-13)
    ideal_h = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    np.testing.assert_allclose(
        change_basis / (np.sqrt(2) * change_basis[0, 0]), ideal_h, atol=2e-13
    )
    return flip, phase, change_basis


def dimers(board, owner):
    edges = []
    for a, b in board.seams:
        if a // 12 == b // 12 == owner:
            orders = [pair for pair in ((a, b), (b, a)) if board.degree[pair[1]] == 5]
            if orders:
                edges.append(min(orders))
    return next(
        group
        for group in combinations(sorted(edges), 5)
        if len({p for edge in group for p in edge}) == 10
    )


def layout():
    graph, shared = shared_report.layout()
    board = federation(0)
    pairs = {key: dimers(board, shared[key]) for key in ("alice", "bob")}
    reader = {}
    for key, bank in pairs.items():
        owner = shared[key]
        local = graph.subgraph(range(12 * owner, 12 * (owner + 1)))
        reader[key] = [
            [(shared_report.path(local, port, g), r) for g, r in shared["readers"][key]]
            for dimer in bank
            for port in dimer
        ]
    protected = {r for standards in shared["readers"].values() for _, r in standards}
    transit = graph.copy()
    transit.remove_nodes_from(protected)
    paths = [
        shared_report.path(transit, src, dst)
        for send, receive in zip(pairs["bob"], pairs["alice"], strict=True)
        for src, dst in zip(send, receive, strict=True)
    ]
    contacts = []
    for trail in paths:
        forward = list(pairwise(trail))
        contacts.extend(forward + forward[:-1][::-1])
    permutation = list(range(board.n))
    crossed_off = permutation.copy()
    for a, b in contacts:
        assert graph.has_edge(a, b)
        permutation[a], permutation[b] = permutation[b], permutation[a]
        if a // 12 == b // 12:
            crossed_off[a], crossed_off[b] = crossed_off[b], crossed_off[a]
    expected = list(range(board.n))
    for path in paths:
        a, b = path[0], path[-1]
        expected[a], expected[b] = expected[b], expected[a]
    assert permutation == expected
    eligible = [
        p
        for p in range(board.n)
        if p // 12 not in (shared["bob"], shared["alice"]) and p not in protected
    ]
    stock = np.ones(board.n, dtype=int)
    for standards in shared["readers"].values():
        stock[standards[0][1]], stock[standards[1][1]] = 0, 2
    stock[eligible[:5]] = 2
    for a, _ in pairs["bob"]:
        stock[a] = 0
    assert stock.sum() == board.n
    geometry = dict(
        owners={k: shared[k] for k in pairs},
        dimers=pairs,
        transfer_paths=paths,
        read_paths=reader,
        stock_doubles=eligible[:5],
        contact_count=len(contacts),
        bridge_contacts=sum(a // 12 != b // 12 for a, b in contacts),
    )
    return geometry, stock, permutation, crossed_off


def read(word, routes, lookup):
    """Twenty restored echo queries, in dimer, endpoint, then standard order."""
    before = word.copy()
    clicks = []
    for queries in routes:
        for path, ref in queries:
            shared_report.transpose(word, path)
            gateway = path[-1]
            key = int(word[gateway]), int(word[ref]), gateway > ref
            clicks.append(int(lookup[key]))
            shared_report.transpose(word, path)
    np.testing.assert_array_equal(word, before)
    return "".join(map(str, clicks))


def decode(clicks):
    if len(clicks) != 20:
        return None
    bits = []
    for start in range(0, 20, 4):
        # At the zero and double standards a vacancy reads 01, a filler 00.
        pair = clicks[start : start + 4]
        if pair not in ("0100", "0001"):
            return None
        bits.append(int(pair == "0001"))
    return sum(bit << j for j, bit in enumerate(bits))


def flights():
    geometry, stock, actual, suppressed = layout()
    lookup = shared_report.certify()
    rows = []
    for launch in range(32):
        source = stock.copy()
        for j, (a, b) in enumerate(geometry["dimers"]["bob"]):
            if launch & (1 << j):
                source[a], source[b] = source[b], source[a]
        sent = read(source, geometry["read_paths"]["bob"], lookup)
        delivered = read(source[actual], geometry["read_paths"]["alice"], lookup)
        blocked = read(source[suppressed], geometry["read_paths"]["alice"], lookup)
        rows.append(
            dict(launch=launch, sender=sent, receiver=delivered, blocked=blocked)
        )
    return geometry, rows
