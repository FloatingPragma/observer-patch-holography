"""Compile the gauge paper's two flux bands and binary admission detector."""

import json
from pathlib import Path
from functools import cache
from math import pi

import numpy as np
from scipy.linalg import block_diag, expm

from model.executable_model import cell, federation
from model.experiments.interval_echo_angular.charged_orientation.flux_doublet.measure import (
    face_edges,
    faces,
    induced_hamiltonian,
    links,
)
from model.experiments.interval_echo_consensus.signed_effects.full_address import (
    measure as address,
)
from .apparatus import ANGLES, K, LAW, SERVICE, GUARD, contact, path, reference

FLUXES = (7, 2)  # vacancy, excess above the same unit background
CHARGES = (-3, 2)
ACCEPT = (1, 3, 5)
RAILS = (*range(12), 238)
N = len(RAILS)


@cache
def phase_pair(vacancy, excess):
    """True local ternary gate, stripping just the unit-background scalar."""
    A = ((vacancy - excess) / 4) % (2 * pi)
    B = (A - vacancy) % (2 * pi)
    half = contact(K, B / (2 * (K + 0.5)))
    flip = reference(K, pi / K)
    gate = flip @ half @ flip @ half @ reference(K, A / K)
    cols = [LAW.SRC.index((n, 1)) for n in range(3)]
    diag = gate[cols, cols]
    np.testing.assert_allclose(
        gate[:, cols], np.eye(9)[:, cols] * diag, atol=1e-13, rtol=0
    )
    out = diag[[0, 2]] / diag[1]
    np.testing.assert_allclose(
        out, np.exp(1j * np.array([vacancy, excess])), atol=1e-13, rtol=0
    )
    return out


@cache
def service_route(port):
    return path(federation(0), RAILS[port], SERVICE, address.RESERVED)


@cache
def contact_route(a, b):
    """Move b next to a, then restore; an occupied transit port is permitted."""
    edges = path(federation(0), RAILS[b], RAILS[a], address.RESERVED)
    return edges[:-1], edges[-1]


def phase(port, vacancy, excess):
    return ("phase", port, vacancy, excess)


def seam(a, b, t):
    return ("seam", a, b, t)


def selected_phase(port, sign, angle):
    return phase(port, angle if sign == 0 else 0, angle if sign == 1 else 0)


def pair_rotation(a, b, sign, theta, arg):
    """Realization of [[c,-conj(z)s],[z s,c]], z=exp(i arg), on one sign."""
    pre = selected_phase(b, sign, pi / 2 - arg)
    post = selected_phase(b, sign, arg - pi / 2)
    toggle = selected_phase(b, 1 - sign, pi)
    return [
        pre,
        seam(a, b, theta),
        toggle,
        seam(a, b, theta),
        toggle,
        post,
        phase(a, 17 * theta, 17 * theta),
        phase(b, 17 * theta, 17 * theta),
    ]


def inverse(word):
    out = []
    for kind, a, b, t in reversed(word):
        if kind == "seam":
            out.append(seam(a, b, (-t) % (2 * pi)))
        else:
            out.append(phase(a, -b, -t))
    return out


def evaluate(word):
    """Evaluate both sectors of the *same* sequence, including non-adjacent routes."""
    ops = np.array([np.eye(N, dtype=complex) for _ in FLUXES])
    for kind, a, b, t in word:
        if kind == "phase":
            service_route(a)
            ops[:, a, :] *= phase_pair(b, t)[:, None]
        else:
            contact_route(a, b)
            gate = contact(K, t % (2 * pi))
            for s, n in enumerate((0, 2)):
                xy = [LAW.SRC.index(x) for x in ((n, 1), (1, n))]
                ops[s, [a, b]] = gate[np.ix_(xy, xy)] @ ops[s, [a, b]]
    return ops


def unitary_word(unitary, sign):
    """Adjacent QR elimination with conditional common-sign Givens pulses."""
    u = unitary.copy()
    reduction = []
    for col in range(N):
        for j in range(N - 1, col, -1):
            a, b = u[j - 1 : j + 1, col]
            if abs(b) < 1e-14:
                continue
            theta = np.arctan2(abs(b), abs(a))
            arg = np.angle(-b) - (np.angle(a) if abs(a) > 1e-15 else 0)
            z = np.exp(1j * arg)
            c, s = np.cos(theta), np.sin(theta)
            R = np.array([[c, -z.conjugate() * s], [z * s, c]])
            u[j - 1 : j + 1] = R @ u[j - 1 : j + 1]
            reduction.extend(pair_rotation(j - 1, j, sign, theta, arg))
    diagonal = np.diag(u)
    np.testing.assert_allclose(u, np.diag(diagonal), atol=1e-12, rtol=0)
    return [
        selected_phase(j, sign, np.angle(d)) for j, d in enumerate(diagonal)
    ] + inverse(reduction)


def logical_frame():
    return json.loads((Path(__file__).parent / "frame.json").read_text())


def flat_word(word):
    """Expand every routed phase/reference and contact to a positive native wait."""
    for kind, a, b, t in word:
        if kind == "phase":
            route = service_route(a)
            A = ((b - t) / 4) % (2 * pi)
            B = (A - b) % (2 * pi)
            local = [
                ("reference", SERVICE, SERVICE, A / K),
                ("seam", SERVICE, GUARD, B / (2 * (K + 0.5))),
                ("reference", SERVICE, SERVICE, pi / K),
                ("seam", SERVICE, GUARD, B / (2 * (K + 0.5))),
                ("reference", SERVICE, SERVICE, pi / K),
            ]
        else:
            route, edge = contact_route(a, b)
            local = [("seam", *edge, t % (2 * pi))]
        for i, j in route:
            yield "swap", i, j, pi
        yield from local
        for i, j in reversed(route):
            yield "swap", i, j, pi


def flat_frame(word):
    """Independently expand and play the frame on 240 physical position rows."""
    board = federation(0)
    initial = np.eye(board.n, dtype=complex)[:, RAILS]
    ops = np.array([initial.copy() for _ in FLUXES])
    count = 0
    swaps = np.eye(9)[[LAW.SRC.index((b, a)) for a, b in LAW.SRC]]
    np.testing.assert_allclose(contact(K, pi), swaps, rtol=0, atol=1e-13)
    for kind, a, b, t in flat_word(word):
        count += 1
        if kind == "reference":
            assert board.degree[a] == 5 and a not in address.RESERVED
            for s, n in enumerate((0, 2)):
                ops[s, a, :] *= np.exp(-1j * K * t * (n * n - 1))
        else:
            assert b in board.neighbours[a]
            if kind == "swap":
                assert a not in address.RESERVED and b not in address.RESERVED
                ops[:, [a, b]] = ops[:, [b, a]]
            else:
                assert not ({a, b} & (address.RESERVED - {GUARD}))
                gate = contact(K, t)
                for s, n in enumerate((0, 2)):
                    xy = [LAW.SRC.index(x) for x in ((n, 1), (1, n))]
                    ops[s, [a, b]] = gate[np.ix_(xy, xy)] @ ops[s, [a, b]]
    outside = [j for j in range(board.n) if j not in RAILS]
    assert np.max(abs(ops[:, outside])) < 2e-12
    return ops[:, RAILS], count


def target(k):
    levels, modes = np.linalg.eigh(induced_hamiltonian(links(k)))
    rank = np.count_nonzero(levels < levels[0] + 1e-8)
    basis = modes[:, :rank]
    P = basis @ basis.conj().T
    return levels, rank, modes, P


def frame(modes, rank):
    """Both codes share one binary pass set; embed the unused spare transparently."""
    pass_slots = list(ACCEPT[:rank])
    reject_slots = [p for p in range(N) if p not in ACCEPT]
    used = pass_slots + reject_slots[: 12 - rank]
    leftover = next(p for p in range(N) if p not in used)
    u = np.zeros((N, N), complex)
    u[used, :12] = modes.conj().T
    u[leftover, 12] = 1
    return u


def membership():
    address.local_check()
    calls = [
        call
        for call in address.schedule()
        if call["bit"] == 0 and call["declared_address"] in ACCEPT
    ]
    E = []
    for sgn in (-1, 1):
        flags = []
        for p in RAILS:
            loads, scalar = address.execute_basis(calls, sgn, p)
            assert abs(scalar - address.COMMON_PHASE ** len(calls)) < 1e-13
            raw = address.read_word(loads)
            assert raw in (0, 1) and np.array_equal(
                loads, address.initial_word(sgn, p, raw)
            )
            flags.append(raw)
        E.append(np.diag(flags))
    return E, calls


def edge_check(resolution):
    board = cell()
    phase_fields = [links(k) for k in FLUXES]
    single_step = [np.eye(12, dtype=complex) for _ in FLUXES]
    physical = [np.zeros((12, 12), complex) for _ in FLUXES]
    dt = 1 / (2 * resolution)
    for idx, (a, b) in enumerate([*board.seams, *reversed(board.seams)]):
        actual = phase_pair(*(2 * pi * field[a, b] / 20 for field in phase_fields))
        service_route(b)
        contact_route(a, b)
        for s, n in enumerate((0, 2)):
            xy = [LAW.SRC.index(x) for x in ((n, 1), (1, n))]
            D = np.diag([1, actual[s]])
            G = D @ contact(K, dt)[np.ix_(xy, xy)] @ D.conj().T
            single_step[s][[a, b]] = G @ single_step[s][[a, b]]
            if idx < len(board.seams):
                local = LAW.seam_term(K)[np.ix_(xy, xy)]
                physical[s][np.ix_([a, b], [a, b])] += D @ local @ D.conj().T
    residual = []
    for i, k in enumerate(FLUXES):
        h = induced_hamiltonian(links(k))
        np.testing.assert_allclose(physical[i], h, atol=2e-12, rtol=0)
        residual.append(
            float(
                np.linalg.norm(
                    np.linalg.matrix_power(single_step[i], resolution) - expm(-1j * h),
                    2,
                )
            )
        )
    return residual


def report():
    targets = [target(k) for k in FLUXES]
    word = logical_frame()
    actual = evaluate(word)
    undo = evaluate(inverse(word))
    flat, waits = flat_frame(word)
    np.testing.assert_allclose(flat, actual, atol=3e-12, rtol=0)
    mask, calls = membership()
    acts = []
    zeno = []
    generators = [[] for _ in range(12)]

    def opnorm(m):
        return float(np.linalg.norm(m, 2))

    for s, (_, r, modes, P) in enumerate(targets):
        np.testing.assert_allclose(undo[s], actual[s].conj().T, atol=8e-12, rtol=0)
        yes = (undo[s] @ mask[s] @ actual[s])[:, :12]
        no = (undo[s] @ (np.eye(N) - mask[s]) @ actual[s])[:, :12]
        err = max(
            opnorm(yes - np.vstack([P, np.zeros(12)])),
            opnorm(no - np.vstack([np.eye(12) - P, np.zeros(12)])),
        )
        assert err < 1e-10
        Pactual = yes[:12]
        V = modes[:, :r]
        q = CHARGES[s]
        h = r / 12
        np.testing.assert_allclose(P.diagonal(), h, atol=1e-12, rtol=0)
        for i in range(12):
            generators[i].append(q * np.outer(V[i].conj(), V[i]))
        bright = P[:, 0] / np.sqrt(h)
        leakage = 1 - np.linalg.norm(P @ (np.diag([-1] + [1] * 11) @ bright)) ** 2
        assert abs(leakage - 4 * h * (1 - h)) < 1e-12
        acts.append(
            dict(
                rank=int(r),
                native_load=(-1, 1)[s],
                encoded_charge=q,
                face_flux=FLUXES[s],
                leverage=h,
                sign_kick_leakage=float(leakage),
                branch_error=err,
            )
        )
        theta = pi / 3
        ideal = expm(-1j * theta * generators[0][s])
        runs = []
        for n in (1, 4, 16, 64):
            D = np.eye(12, dtype=complex)
            D[0, 0] = phase_pair(*(-c * theta / n for c in CHARGES))[s]
            Kpass = V.conj().T @ np.linalg.matrix_power(Pactual @ D, n) @ V
            runs.append(
                dict(
                    n=n,
                    error=opnorm(Kpass - ideal),
                    bright_survival=float(
                        np.linalg.norm(Kpass @ (V.conj().T @ bright)) ** 2
                    ),
                )
            )
        assert all(y["error"] < x["error"] for x, y in zip(runs, runs[1:]))
        zeno.append(runs)
    matrices = [block_diag(*blocks) for blocks in generators]
    gram = np.array([[np.trace(x @ y).real for y in matrices] for x in matrices])
    assert (
        np.linalg.matrix_rank(gram, tol=1e-8) == 12
        and max(abs(np.trace(m)) for m in matrices) < 1e-12
    )
    for triangle in faces():
        assert (
            abs(
                np.prod(
                    [
                        np.exp(1j * (ANGLES[a] - ANGLES[b]))
                        for a, b in face_edges(triangle)
                    ]
                )
                - 1
            )
            < 1e-12
        )
    assert target(0)[1] == 1
    errors = {str(n): edge_check(n) for n in (8, 16, 32)}
    return dict(
        sectors=acts,
        compiled_frame_logical_pulses=len(word),
        membership_native_waits=len(address.pulse_word(calls)),
        minimum_current_gram_eigenvalue=float(np.linalg.eigvalsh(gram)[0]),
        gram_rank=12,
        symmetric_product_unit_time_errors=errors,
        zeno=zeno,
        expanded_frame_native_waits=waits,
        zero_flux_leading_rank=1,
    )
