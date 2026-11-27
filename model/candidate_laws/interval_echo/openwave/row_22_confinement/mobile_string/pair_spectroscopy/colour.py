"""A signed stored pair and the calibrated two-owner conjugate colour receiver."""

import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls import (
    orbital_singlet as orbital,
)
from model.executable_model import carrier

SHOTS = 2**22


def colour_bases():
    j, k = np.indices((3, 3))
    omega = np.exp(2j * np.pi / 3)
    return [np.eye(3)] + [
        omega ** ((m * j * j + j * k) % 3) / np.sqrt(3) for m in range(3)
    ]


def native_analyzer(frame, path):
    """Adjacent native pulses and single-port phases send three admitted columns to rails."""
    cols = frame.copy().astype(complex)
    transform = np.eye(len(cols), dtype=complex)
    dwells = []
    for col in range(3):
        for pos in range(len(path) - 2, col - 1, -1):
            sites = np.array(path[pos : pos + 2])
            a, b = cols[sites, col]
            if abs(b) < 1e-14:
                continue
            theta = np.arctan2(abs(b), abs(a))
            t = 2 * theta
            # Full local loaded 01/10 or 12/21 generator, (8.5 I-X/2).
            co, si = np.cos(theta), 1j * np.sin(theta)
            pulse = np.exp(-8.5j * t) * np.array([[co, si], [si, co]])
            prep = np.diag([np.exp(-1j * np.angle(a)), -1j * np.exp(-1j * np.angle(b))])
            post = np.exp(8.5j * t) * np.diag([1, 1j])
            gate = post @ pulse @ prep
            cols[sites] = gate @ cols[sites]
            transform[sites] = gate @ transform[sites]
            dwells.append((sites.tolist(), float(t)))
    response = transform[path[:3]] @ frame
    np.testing.assert_allclose(abs(response), np.eye(3), atol=3e-13)
    return transform[path[:3]], len(dwells)


def path_on_owner(board, c):
    graph = nx.Graph(board.seams).subgraph(range(12 * c, 12 * c + 12))

    def extend(path):
        if len(path) == 12:
            return path
        for v in sorted(set(graph[path[-1]]) - set(path)):
            found = extend([*path, v])
            if found:
                return found
        return None

    return [v - 12 * c for v in extend([12 * c])]


def layout():
    board = carrier.federation(1)
    transforms, locations = orbital.placements(board)
    ca, cb = locations[0][1][:2]
    _, frame, _ = orbital.compact_levels((0, 1, 2))
    a, b = (frame[list(transforms[c])] for c in (ca, cb))
    # Use the same canonical row at both receiving sites. Exact equal-norm rows
    # use the smallest label. The tolerance separates rational projector diagonals.
    norm = (frame**2).sum(axis=1)
    row = int(np.flatnonzero(norm >= norm.max() - 1e-12)[0])
    u = 12 * ca + list(transforms[ca]).index(row)
    v = 12 * cb + list(transforms[cb]).index(row)
    partner = min(q for q in board.neighbours[u] if q // 12 == ca)
    graph = nx.Graph(board.seams)
    graph.remove_node(u)
    route = nx.shortest_path(graph, partner, v)
    gate = min(q for q in range(12 * cb, 12 * cb + 12) if board.degree[q] == 6)
    return board, (ca, cb), (a, b), (u, partner, v, gate), route


def route_product(board, vacancy, path):
    native = expm(-1j * np.pi * law.seam_term(8))
    loads = np.ones(board.n, int)
    loads[vacancy], loads[path[0]] = 0, 2
    for a, b in zip(path, path[1:]):
        if b not in board.neighbours[a] or vacancy in (a, b):
            raise ValueError("Invalid physical extraction path")
        before = (int(loads[a]), int(loads[b]))
        after = before[::-1]
        column = native[:, law.SRC.index(before)]
        np.testing.assert_allclose(
            column, np.eye(9)[:, law.SRC.index(after)], atol=4e-14
        )
        loads[a], loads[b] = after
    return int(np.flatnonzero(loads == 0)[0]), int(np.flatnonzero(loads == 2)[0])


def graph_check(board, owners, frames):
    full = []
    residual = []
    for q, owner, z in zip((-1, 1), owners, frames, strict=True):
        v = np.zeros((board.n, 3))
        v[12 * owner : 12 * owner + 12] = z
        hdiag = 8 * (6 * q * q + 2 * q * (6 - board.degree))
        energy = 3 + 8 * (6 * q * q + 2 * q)
        residual.append(
            float(abs(board.laplacian() @ v / 2 + (hdiag - energy)[:, None] * v).max())
        )
        full.append(v)
    row, col = np.array(board.seams).T
    # Amplitudes for each of the nine two-owner colour products on each oriented seam.
    u = np.r_[row, col]
    v = np.r_[col, row]
    local = full[0][u, :, None] * full[1][v, None, :]
    # Noninteracting product action is exact: no two nonzero components occupy
    # the same port or touch. This also holds for the whole 9-dimensional span.
    overlap = abs(full[0][:, :, None] * full[1][:, None, :]).max()
    assert overlap == 0 and abs(local).max() == 0
    return dict(
        oriented_seams=len(u),
        pairs_checked=9 * len(u),
        max_local_overlap=float(abs(local).max()),
        single_residuals=residual,
        energies=[35, 67],
        joint_energy=102,
    )


def memory_check():
    source = np.eye(law.DIM)[:, law.SRC.index((0, 2))]
    uniform = np.eye(law.DIM)[:, law.SRC.index((1, 1))]
    echo = law.echo()
    stored = echo @ source
    hx, hc = law.local_term_parts()
    mark = np.repeat([0] + [int(abs(d) == 2) for _, d in law.LABELS], 9)
    blank = np.repeat([1] + [0] * 10, 9)

    def comm(h):
        return float(abs((blank[:, None] - blank[None, :]) * h).max())

    def reads(psi):
        probs = abs(psi.reshape(11, 9)) ** 2
        return dict(
            port_11=float(probs[:, law.SRC.index((1, 1))].sum()),
            marked=float(mark @ abs(psi) ** 2),
            sign_bins=[
                float(
                    sum(
                        probs[i + 1].sum()
                        for i, (_, tag) in enumerate(law.LABELS)
                        if tag == d
                    )
                )
                for d in (-2, 0, 2)
            ],
            released_02=float(probs[0, law.SRC.index((0, 2))]),
        )

    return dict(
        tag=reads(stored),
        opposite_tag=reads(echo @ np.eye(law.DIM)[:, law.SRC.index((2, 0))]),
        sham=reads(uniform),
        inverse_on_tag=reads(echo @ stored),
        inverse_on_sham=reads(echo @ uniform),
        tag_free_residual=float(np.linalg.norm((hx + 8 * hc) @ stored)),
        blank_block_commutators=[comm(hx), comm(hc)],
    )


def tomography(rho, readers, rng):
    counts = []
    for a, b in readers:
        u = np.kron(a, b)
        weights = (u @ rho @ u.conj().T).diagonal().real
        probs = np.r_[weights, 1 - weights.sum()].clip(0.0)
        counts.append(rng.multinomial(SHOTS, probs / probs.sum()).tolist())
    # Each setting contributes equal-colour minus one quarter of admission.
    # Bounded sample range is 1 in each of four independent settings.
    estimate = sum(
        sum(c[i] for i in (0, 4, 8)) - (SHOTS - c[-1]) / 4 for c in counts
    ) / (3 * SHOTS)
    return dict(counts=counts, inferred_singlet=estimate)


def acquisition():
    board, owners, frames, sites, route = layout()
    a, b = frames
    ca, cb = owners
    u, partner, v, gate = sites
    pulses = []
    readers = []
    for basis in colour_bases():
        ea, na = native_analyzer(a @ basis, path_on_owner(board, ca))
        eb, nb = native_analyzer(b @ basis.conj(), path_on_owner(board, cb))
        readers.append((ea @ a, eb @ b))
        pulses.append([na, nb])
    endpoint, positive = route_product(board, u, route)
    graph = nx.Graph(board.seams)
    graph.remove_node(u)
    gateway_path = nx.shortest_path(graph, partner, gate)
    _, gateway = route_product(board, u, gateway_path)
    direction = np.kron(a[endpoint % 12], b[positive % 12])
    gateway_vector = np.kron(a[endpoint % 12], b[gateway % 12])
    direct = np.outer(direction, direction)
    singlet = np.eye(3).reshape(-1) / np.sqrt(3)
    coherent = np.outer(singlet, singlet)
    mixture = np.diag(singlet**2)
    absent = np.zeros((9, 9))
    turned = np.kron(np.eye(3), np.roll(np.eye(3), 1, axis=0))

    def phase_state(x):
        assert abs(np.linalg.norm(x) - 1) < 1e-12
        return np.outer(x, x.conj())

    rng = np.random.default_rng(181832)
    backgrounds = {
        "singlet": coherent,
        "incoherent_same_colour": mixture,
        "one_end_turn": phase_state(turned @ singlet),
        "blank": absent,
        "written_11": absent,
        "local_inverse_02": absent,
        "inverse_blank_11": absent,
        "routed_02": direct,
        "gateway_02": np.outer(gateway_vector, gateway_vector),
    }
    result = {}
    for key, rho in backgrounds.items():
        result[key] = {}
        times = (
            (0, 1, 16)
            if key in ("routed_02", "local_inverse_02", "written_11")
            else (0,)
        )
        for t in times:
            # Exact compressed all-contact phase e^(-102 i t) times identity;
            # hence the admitted density and the complete pair refusal are stationary.
            result[key][str(t)] = tomography(rho, readers, rng)
    count = sum(len(z) for z in result.values())
    radius = float(np.sqrt(2 * np.log(2 * count / 0.01) / SHOTS) / 3)
    # Begin a separate attachment to the already prepared colour singlet. All
    # nonzero orbitals avoid the distant tag owner and its two interior sites.
    edge = min(
        (a, b) for a, b in board.seams if a // 12 == b // 12 and a // 12 not in owners
    )
    tagged = {}
    tag_rng = np.random.default_rng(181833)
    memory = memory_check()
    for key, read in [
        ("blank", memory["sham"]),
        ("written", memory["tag"]),
        ("inverse_written", memory["inverse_on_tag"]),
    ]:
        obs = tomography(coherent, readers, tag_rng)
        obs.update(
            tag_bin="written -2" if read["marked"] > 0.5 else "no magnitude-two mark",
            remote_load_word="02" if read["released_02"] > 0.5 else "11",
        )
        tagged[key] = obs
    return dict(
        seed=181832,
        shots_per_setting=SHOTS,
        owners=list(owners),
        sites=list(map(int, sites)),
        routed_positive_ports=list(map(int, route)),
        active_colour_turn_counts=pulses,
        gateway_route=list(map(int, gateway_path)),
        colour_checks=graph_check(board, owners, frames),
        memory=memory,
        independent_singlet_tag=dict(
            seed=181833,
            ports=list(map(int, edge)),
            motion="pre-existing compact pair; reads immediately after the stated intervention",
            common_untouched_ordinary_pair_energy=102,
            records=tagged,
        ),
        joint_99_radius=radius,
        detector="nine resolved colour pairs then joint complement",
        routed_exact=float(singlet @ direct @ singlet),
        records=result,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.parse_args().output.write_text(
        json.dumps(acquisition(), indent=1, allow_nan=False) + "\n"
    )
