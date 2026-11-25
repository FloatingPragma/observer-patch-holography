"""A geometric flux, its charged orbital and the native written-collar comparator."""

import argparse
from collections import deque
from functools import cache
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation, rotations
from model.executable_model.carrier import icosahedron_from_coordinates, isomorphism
from model.experiments.interval_echo_lepton_register.measure import contact, phase
from ..measure import compile_turn, path
from .passive import native_flux_test
from .bands import band_counts
from model.experiments.interval_echo_lepton_register import measure as koide
from sympy import Poly, Symbol, cyclotomic_poly

HERE = Path(__file__).resolve().parent
K = candidate.KAPPA
LAW = candidate.native_law
CELL = cell()
MOD = 20
N = CELL.n
TAU = 2 * np.pi
SHOTS = 2**18


def face_edges(face):
    return list(zip(face, (*face[1:], face[0])))


@cache
def faces():
    points, edges = icosahedron_from_coordinates()
    points = points[isomorphism(CELL.seams, edges)]
    result = []
    for a, b, c in CELL.faces():
        if np.linalg.det(points[[a, b, c]]) < 0:
            b, c = c, b
        result.append((a, b, c))
    return tuple(result)


def links(q, root=0):
    """Set integer edge phases by routing each unit of face flux to a root face."""
    triangles = faces()
    parent = {root: None}
    walk = [root]
    for f in walk:
        for edge in face_edges(triangles[f]):
            g = next(
                i for i, face in enumerate(triangles) if edge[::-1] in face_edges(face)
            )
            if g not in parent:
                parent[g] = (f, edge[::-1])
                walk.append(g)
    size = np.ones(MOD, dtype=int)
    link = np.zeros((N, N), dtype=int)
    for f in reversed(walk[1:]):
        g, (a, b) = parent[f]
        link[a, b] = q * size[f]
        link[b, a] = -link[a, b]
        size[g] += size[f]
    for face in triangles:
        assert sum(link[a, b] for a, b in face_edges(face)) % MOD == q % MOD
    return link


def compensated_turn(permutation, link):
    """The permuted field agrees after a local phase u_j-u_i measured on each seam."""
    labels = {0: 0}
    todo = deque([0])
    while todo:
        i = todo.popleft()
        for j in CELL.neighbours[i]:
            phase_j = (
                labels[i] + link[permutation[i], permutation[j]] - link[i, j]
            ) % MOD
            if j in labels:
                assert labels[j] == phase_j
            else:
                labels[j] = int(phase_j)
                todo.append(j)
    return np.array([labels[i] for i in range(N)])


@cache
def half_pairs():
    ids = np.arange(N)
    halves = [
        p
        for p in rotations()
        if tuple(ids) != p and np.array_equal(np.array(p)[list(p)], ids)
    ]
    return tuple(
        (r, s)
        for r, s in combinations(halves, 2)
        if all(r[s[i]] == s[r[i]] for i in range(N))
    )


def turn_word(p, phases, owner=0):
    # The compiler supplies positive native π interchanges on blank collars.
    out = np.eye(N, dtype=complex)
    for i, angle in enumerate(phases):
        phase(out, i, TAU * int(angle) / MOD, 1, owner=owner)
    for i, j in compile_turn(p):
        contact(out, i, j, np.pi, 1)
    return out


def spatial_order(link):
    frames = {r: compensated_turn(r, link) for r in rotations()}
    phases = []
    for r, s in half_pairs():
        relative = (
            frames[s] + frames[r][list(s)] - frames[r] - frames[s][list(r)]
        ) % MOD
        assert len(set(relative)) == 1
        phases.append(int(relative[0]))
    return phases


def induced_hamiltonian(link, speed=1):
    h = np.zeros((N, N), complex)
    exchange, contrast = LAW.local_term_parts()
    native = (speed * exchange + K * contrast)[:9, :9]
    i, j = map(LAW.SRC.index, [(0, 1), (1, 0)])
    for a, b in CELL.seams:
        h[a, a] += native[i, i]
        h[b, b] += native[j, j]
        h[b, a] = native[j, i] * np.exp(1j * TAU * link[a, b] / MOD)
        h[a, b] = h[b, a].conjugate()
    return h


def split_steps(link, steps):
    """Strang sequence. Contact sums use total simulated exposure one."""
    dt = 1 / steps
    half_step = np.eye(N, dtype=complex)
    for a, b in [*CELL.seams, *reversed(CELL.seams)]:
        angle = TAU * link[a, b] / MOD
        phase(half_step, b, -angle, 1)
        contact(half_step, a, b, dt / 2, 1)
        phase(half_step, b, angle, 1)
    return np.linalg.matrix_power(half_step, steps)


def sample(probability, rng):
    p = np.maximum(0, np.round(probability.reshape(-1), 12))
    return rng.multinomial(SHOTS, p / p.sum()).reshape(probability.shape).tolist()


def arm_read(left, right, rng):
    out = abs(np.stack((left + right, left - right))) ** 2 / 4
    return dict(probability=out.round(12).tolist(), counts=sample(out, rng))


def band_experiment(q, rng):
    link = links(q)
    h = induced_hamiltonian(link)
    e, v = np.linalg.eigh(h)
    ground = v[:, e < e[0] + 1e-8]
    projector = ground @ ground.conj().T
    # Real delta-station sources projected by the exact spectral response.
    sharp = (
        [
            0,
            next(
                i
                for i in range(1, N)
                if abs(projector[0, i]) < projector[0, 0].real - 1e-6
            ),
        ]
        if ground.shape[1] > 1
        else [0]
    )
    ports = sharp
    r, s = half_pairs()[0]
    up, vp = (compensated_turn(g, link) for g in (r, s))
    R, S = turn_word(r, up), turn_word(s, vp)
    plain_r, plain_s = turn_word(r, np.zeros(N)), turn_word(s, np.zeros(N))
    np.testing.assert_allclose(R @ h, h @ R, atol=1e-11)
    np.testing.assert_allclose(S @ h, h @ S, atol=1e-11)
    shots = []
    for port in ports:
        source = projector[:, port] / np.sqrt(projector[port, port].real)
        shots.append(
            dict(
                port=port,
                ideal_spectral_admission=float(projector[port, port].real),
                compensated=arm_read(R @ S @ source, S @ R @ source, rng),
                bare=arm_read(
                    plain_r @ plain_s @ source, plain_s @ plain_r @ source, rng
                ),
            )
        )
    phases = spatial_order(link)
    assert phases == spatial_order(links(q, 13))
    # A local rephasing of every link changes each lift, but cancels in its commutator.
    rechart = np.arange(N) ** 2
    assert phases == spatial_order(link + rechart[None, :] - rechart[:, None])
    return dict(
        quanta=q,
        face_phase_units=int(q),
        spectrum=e.round(10).tolist(),
        ground_dimension=ground.shape[1],
        gap=float(e[ground.shape[1]] - e[0]),
        commutator_units=phases,
        gram_determinant=float(np.linalg.det(projector[np.ix_(ports, ports)]).real),
        reads=shots,
    )


def release(link):
    h = induced_hamiltonian(link)
    e, v = np.linalg.eigh(h)
    basis = v[:, e < e[0] + 1e-8]
    board = federation(0)
    initial = np.zeros((board.n, len(basis.T)), complex)
    initial[:N] = basis
    native = board.laplacian() / 2 + diags(K * (2 * board.degree - 6))
    evolved = expm_multiply(-1j * native, initial)
    retention = np.linalg.svd(basis.conj().T @ evolved[:N], compute_uv=False) ** 2
    return dict(
        native_time=1,
        minimum_retention=float(min(retention)),
        maximum_retention=float(max(retention)),
    )


def character_census(link):
    """Character norm of the magnetic permutation action, reduced in Q[zeta_20]."""
    power = np.zeros(MOD, dtype=int)
    for r in rotations():
        u = compensated_turn(r, link)
        fixed = [v for v in range(N) if r[v] == v]
        for i in fixed:
            for j in fixed:
                power[(u[i] - u[j]) % MOD] += 1
    z = Symbol("z")
    poly = Poly(sum(int(c) * z**i for i, c in enumerate(power)), z)
    residue = poly.rem(Poly(cyclotomic_poly(MOD, z), z))
    assert residue.degree() == 0 and int(residue.TC()) % 60 == 0
    return dict(
        trace_product_coefficients=power.tolist(),
        commutant_dimension=int(residue.TC()) // 60,
    )


def filtered_band(link, rng):
    h = induced_hamiltonian(link)
    e, v = np.linalg.eigh(h)
    band = v[:, e < e[0] + 1e-8]
    p = band @ band.conj().T
    u = split_steps(link, 128)
    good = (np.eye(N) + np.exp(1j * e[0]) * u) / 2
    bad = (np.eye(N) - np.exp(1j * e[0]) * u) / 2
    r, s = half_pairs()[0]
    R, S = [turn_word(g, compensated_turn(g, link)) for g in (r, s)]
    panels = []
    for port in (0, 1):
        source = np.eye(N, dtype=complex)[:, port]
        refuse = []
        for _ in range(12):
            refuse.append(np.linalg.norm(bad @ source) ** 2)
            source = good @ source
        admitted = float(np.vdot(source, source).real)
        phase_reads = (
            abs(
                np.stack(
                    (R @ S @ source + S @ R @ source, R @ S @ source - S @ R @ source)
                )
            )
            ** 2
            / 4
        )
        complete = np.r_[refuse, phase_reads.ravel()]
        np.testing.assert_allclose(complete.sum(), 1, atol=3e-10)
        panels.append(
            dict(
                port=port,
                spectral_appointment=float(e[0]),
                admission=admitted,
                band_fidelity=float(np.linalg.norm(p @ source) ** 2 / admitted),
                all_probabilities=complete.round(12).tolist(),
                counts=sample(complete, rng),
            )
        )
    edge = koide.addresses(CELL.seams[0])
    analyzer = koide.reader(edge, 2 / 9, 1)
    a, _, plus, minus = edge
    modes = analyzer.conj().T[:, [a, minus, plus]]
    overlap = band.conj().T @ modes
    weights = (abs(overlap) ** 2).sum(axis=0)
    common = np.linalg.eigvalsh(overlap.conj().T @ overlap)
    return dict(
        filters=panels,
        family_character_order=["heavy", "light", "middle"],
        fraction_in_ground_band=weights.tolist(),
        principal_overlap_squared=common.tolist(),
    )


def charged_source(owner=0):
    """Return a reproducible port frame for the owner-zero magnetic doublet.

    Phases use the root-face-zero chart. Gram-Schmidt on projections of ports
    0 and 1 fixes the two vectors without choosing eigensolver phases.
    """
    link = links(1)
    energies, eigenvectors = np.linalg.eigh(induced_hamiltonian(link))
    eigenbasis = eigenvectors[:, energies < energies[0] + 1e-8]
    projector = eigenbasis @ eigenbasis.conj().T
    first = projector[:, 0] / np.sqrt(projector[0, 0].real)
    second = projector[:, 1] - first * np.vdot(first, projector[:, 1])
    columns = np.column_stack((first, second / np.linalg.norm(second)))
    turns = tuple(
        turn_word(g, compensated_turn(g, link), owner=owner) for g in half_pairs()[0]
    )
    return columns, projector, turns


def exchange_in_band(rng):
    board = federation(0)
    _, other = next(e for e in board.external if e[0] < 12)
    src = list(range(N))
    dst = [N * (other // N) + p for p in src]
    mapping = np.arange(board.n)
    gates = []
    for a, b in zip(src, dst, strict=True):
        walk = path(board.neighbours, a, b)
        edges = list(zip(walk, walk[1:]))
        for x, y in edges + edges[:-1][::-1]:
            gates.append((x, y))
            mapping[[x, y]] = mapping[[y, x]]
    np.testing.assert_equal(mapping[src], dst)
    np.testing.assert_equal(mapping[dst], src)
    others = sorted(set(range(board.n)) - set(src + dst))
    np.testing.assert_equal(mapping[others], others)
    # Every π contact, including two occupied or two empty ends, is an
    # ordinary ternary interchange. Apply the computed complete port map to
    # two canonically ordered vacancy addresses, retaining both endpoint reads.
    swap = np.eye(9)[:, [LAW.SRC.index((b, a)) for a, b in LAW.SRC]]
    np.testing.assert_allclose(expm(-1j * np.pi * LAW.seam_term(K)), swap, atol=5e-13)
    basis, _, _ = charged_source()
    product_basis = np.kron(basis, basis)
    after = np.zeros_like(product_basis)
    for i, a in enumerate(src):
        for j, b in enumerate(dst):
            ia, ib = src.index(int(mapping[b])), dst.index(int(mapping[a]))
            after[N * ia + ib] = product_basis[N * i + j]
    effective = product_basis.conj().T @ after
    np.testing.assert_allclose(effective, np.eye(4)[[0, 2, 1, 3]], atol=1e-12)
    same = np.array([1, 0, 0, 0])
    opposite = np.array([0, 1, 0, 0])
    anti = np.array([0, 1, -1, 0]) / np.sqrt(2)
    reads = []
    for name, state in [
        ("aligned", same),
        ("unlike", opposite),
        ("antisymmetric", anti),
    ]:
        initial = product_basis @ state
        final = after @ state
        prob = (
            np.array(
                [
                    np.linalg.norm(initial + final) ** 2,
                    np.linalg.norm(initial - final) ** 2,
                ]
            )
            / 4
        )
        reads.append(
            dict(
                source=name,
                exit_probabilities=prob.round(12).tolist(),
                counts=sample(prob, rng),
            )
        )
    return dict(
        owners=[0, other // N],
        native_pi_contacts=len(gates),
        exchange_gates=gates,
        reads=reads,
    )


def run():
    rng = np.random.default_rng(47015)
    panels = [band_experiment(q, rng) for q in (0, 1, -1, 2)]
    link = links(1)
    h = induced_hamiltonian(link)
    target = expm(-1j * h)
    compiled = [
        dict(
            slices=k,
            operator_error=float(np.linalg.norm(split_steps(link, k) - target, 2)),
        )
        for k in (8, 32, 128)
    ]
    return dict(
        shots_per_bank=SHOTS,
        integer_phase_denominator=MOD,
        panels=panels,
        compiled_contact_error=compiled,
        guard_release=release(link),
        rotation_character=character_census(link),
        certified_bands=band_counts(link, CELL.seams),
        coherent_preparations=filtered_band(link, rng),
        aligned_doublet_exchange=exchange_in_band(rng),
        passive_collar=native_flux_test(rng, SHOTS),
        exchange_times_three_spectrum=np.linalg.eigvalsh(induced_hamiltonian(link, 3))
        .round(10)
        .tolist(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = run()
    if parser.parse_args().write:
        (HERE / "results.json").write_text(json.dumps(data, indent=2) + "\n")
    else:
        print(json.dumps(data, indent=2))
