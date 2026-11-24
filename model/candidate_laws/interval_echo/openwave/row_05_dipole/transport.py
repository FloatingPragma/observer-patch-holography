"""Charge and current transported in the orthogonal mobile relative of a compact orbital."""

import argparse
from functools import lru_cache
import json
import networkx as nx
import numpy as np
from scipy import sparse
from scipy.linalg import null_space
from scipy.sparse.linalg import expm_multiply
from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from .experiment import HERE, CHARGE, defect_generator, gates_of, portable, source

SHOTS = 2**27
TIMES = (0, 2, 8, 32, 64, 128)


def draw(probs, rng):
    p = np.clip(np.round(probs, 12), 0, 1)
    if abs(sum(p) - 1) > 1e-9:
        raise ArithmeticError("Incomplete mobile record")
    return rng.multinomial(SHOTS, p / p.sum()).tolist()


@lru_cache(maxsize=None)
def modes(gates):
    compact, areas, _, moment = source(gates)
    constraints = cell().adjacency().toarray() + np.eye(12)
    dark = null_space(np.vstack((constraints, np.eye(12)[list(gates)])))
    bright = null_space(np.vstack((constraints, dark.T)))
    values, u = np.linalg.eigh(bright.T @ moment @ bright)
    if values[-1] < 1e-10:
        return dark, None, areas
    p = bright @ u[:, -1]
    pivot = np.flatnonzero(abs(p) > max(abs(p)) - 1e-10)[0]
    p *= np.exp(-1j * np.angle(p[pivot]))
    return dark, p, areas


def menu(v, board, owner, rng):
    gates = gates_of(board, owner)
    dark, p, areas = modes(gates)
    weight = float(np.vdot(v, v).real)
    compact = float(np.sum(abs(dark.T @ v) ** 2))
    plus = abs(np.vdot(p, v)) ** 2 if p is not None else 0
    minus = abs(p @ v) ** 2 if p is not None else 0
    pc = draw([plus, minus, compact, weight - plus - minus - compact, 1 - weight], rng)
    currents = []
    for a, b in cell().seams:
        pos = abs(v[a] - 1j * v[b]) ** 2 / 2
        neg = abs(v[a] + 1j * v[b]) ** 2 / 2
        currents.append(draw([pos, neg, 1 - pos - neg], rng))
    cr = np.array(currents)
    mu = CHARGE * (cr[:, 0] - cr[:, 1]) @ areas / (2 * SHOTS)
    return dict(mode_counts=pc, current_counts=currents, moment=mu)


def links(board):
    graph = nx.Graph((a // 12, b // 12) for a, b in board.external)
    paths = [
        min(
            tuple(p)
            for v, d in nx.single_source_shortest_path_length(
                graph, 0, cutoff=2
            ).items()
            if d == k
            for p in nx.all_shortest_paths(graph, 0, v)
        )
        for k in (1, 2)
    ]
    return paths


def kick(board, target, h, rng):
    a, b = next((a, b) for a, b in board.external if a // 12 == 0 and b // 12 == target)
    g = gates_of(board, 0)
    z = gates_of(board, target)
    records = []
    for label, u, v in [
        ("bright", modes(g)[1], modes(z)[1]),
        ("compact", source(g)[0], source(z)[0]),
    ]:
        if u is None or v is None:
            continue
        phase = (
            np.angle(u[a]) - np.angle(v[b % 12])
            if min(abs(u[a]), abs(v[b % 12])) > 1e-10
            else 0
        )
        v = v * np.exp(1j * phase)
        mat = -0.5 * u[a].conj() * v[b % 12]
        vec = np.zeros((board.n, 3), complex)
        for j, q in enumerate((0, np.pi / 2, np.pi)):
            vec[:12, j] = u / np.sqrt(2)
            vec[12 * target : 12 * target + 12, j] = np.exp(1j * q) * v / np.sqrt(2)
        acted = h @ vec
        en = np.sum(vec.conj() * acted, axis=0).real
        width = np.linalg.norm(acted - vec * en, axis=0)
        rows = []
        for t in (1 / 128, 1 / 64):
            out = expm_multiply(-1j * t * h, vec)
            for j in range(3):
                vout = out[12 * target : 12 * target + 12, j]
                inside = np.vdot(vout, vout).real
                ov = np.vdot(vec[:, j], out[:, j])
                counts = [
                    draw([(1 + x) / 2, (1 - x) / 2], rng) for x in (ov.real, ov.imag)
                ]
                decoded = np.array(counts) @ np.array([1, -1]) / SHOTS
                rows.append(
                    dict(
                        arm=j,
                        time=t,
                        ramsey_xy=counts,
                        energy_read=-np.arctan2(decoded[1], decoded[0]) / t,
                        receiver_counts=draw([inside, 1 - inside], rng),
                    )
                )
        records.append(
            dict(
                orbital=label,
                gateway=[a, b],
                matrix=[mat.real, mat.imag],
                frame_angle=phase,
                mean_energies=en.tolist(),
                widths=width.tolist(),
                reads=rows,
            )
        )
    return records


def study(board, k, rng):
    gates = gates_of(board, 0)
    c = source(gates)[0]
    _, p, _ = modes(gates)
    if p is None:
        raise ValueError("Canonical bright orbital is not chiral")
    h = defect_generator(board, -1, k, 0).astype(complex) - (3 + 4 * k) * sparse.eye(
        board.n
    )
    v = np.zeros((board.n, 5), complex)
    v[:12] = np.column_stack([p, p.conj(), c, c.conj(), (p + c) / np.sqrt(2)])
    dst = links(board)
    # Only the two bright packets have to be evolved. The compact states have
    # zero shifted energy, which the full material operator verifies directly.
    assert np.linalg.norm(h @ v[:, 2:4]) < 1e-10
    wave = expm_multiply(-1j * h, v[:, :2], start=0, stop=128, num=65)
    reads = []
    for t in TIMES:
        z = wave[t // 2]
        wave_t = np.column_stack([z, v[:, 2:4], (v[:, 2] + z[:, 0]) / np.sqrt(2)])
        for path in dst:
            target = path[-1]
            local = wave_t[12 * target : 12 * target + 12]
            reads.append(
                dict(
                    time=t,
                    path=path,
                    sources=[menu(col, board, target, rng) for col in local.T],
                )
            )
    return dict(
        owners=board.cells, kappa=k, reads=reads, kicks=kick(board, dst[0][-1], h, rng)
    )


def run():
    rng = np.random.default_rng(510245)
    data = [
        study(federation(level), k, rng)
        for level in (0, 2)
        for k in (candidate.KAPPA, 0)
    ]
    # Simultaneous error for individual signed scores; use twice this
    # radius for a difference between directions of the source.
    n = 4 * len(TIMES) * 2 * 5 * 2
    log = np.log(2 * n / 0.01)
    return dict(
        shots=SHOTS,
        seed=510245,
        sources=[
            "bright_plus",
            "bright_minus",
            "compact_plus",
            "compact_minus",
            "mixed_plus",
        ],
        moment_radius=float(np.sqrt(log / SHOTS)),
        mode_score_radius=float(np.sqrt(2 * log / SHOTS)),
        boards=data,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    data = portable(run())
    path = HERE / "transport.json"
    if parser.parse_args().write:
        path.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    elif json.loads(path.read_text()) != data:
        raise ValueError("Transport read differs")
    print(
        "Verified charge-resolved bright circulation, compact controls and inter-owner phase kicks."
    )


if __name__ == "__main__":
    main()
