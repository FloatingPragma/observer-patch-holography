"""Can a phased incident load or a coherent local receiver distinguish the hidden phase?"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .measure import (
    DEPTH,
    HERE,
    PHASE_STATES,
    SEED,
    SHOTS,
    TIMES,
    angular_frame,
    compile_wave,
    federation,
    inputs,
    operators,
    portals,
    product_ball,
    seed_packet,
)


def receiver_blocks(packet, first, second, receiver, board):
    """Traces the two pure source branches over all ports outside the receiver."""
    inside_a = first // 12 == receiver
    inside_b = second // 12 == receiver
    both = inside_a & inside_b
    neither = ~(inside_a | inside_b)
    singles = inside_a ^ inside_b
    local = np.where(inside_a[singles], first[singles], second[singles]) % 12
    exterior = np.where(inside_a[singles], second[singles], first[singles])
    waves = np.zeros((2, 12, board.n), complex)
    waves[:, local, exterior] = packet[singles].T
    pair = np.zeros((2, 66, 1), complex)
    pair[:, portals(first[both], second[both], 12 * receiver) - 13, 0] = packet[both].T
    empty = packet[neither].T[:, None, :]

    def component_matrices(branches):
        a, b = branches
        ab = a @ b.conj().T
        mean = (a @ a.conj().T + b @ b.conj().T) / 2
        y = (ab - ab.conj().T) / (2j)
        return mean, y

    return [component_matrices(vectors) for vectors in (empty, waves, pair)]


def local_readers(board, gate, cell):
    start = gate // 12 * 12
    cohort = range(start, start + 12)
    seams = [
        (a - start, b - start)
        for a in cohort
        for b in board.neighbours[a]
        if a < b < start + 12
    ]
    standards = np.eye(12)
    labels, waves = [], []
    for a, b in seams:
        for name, c in (("real", 1.0), ("circular", 1j)):
            labels.append(dict(local_ports=[int(a), int(b)], analyzer=name))
            waves.append(
                compile_wave(cell, (standards[a] + c * standards[b]) / np.sqrt(2))
            )
    return labels, np.array(waves)


def count_binary(laws, rng):
    assert min(laws.min(), 1 - laws.max()) > -1e-12
    yes = rng.binomial(SHOTS, np.clip(laws, 0, 1))
    return np.stack([yes, SHOTS - yes], axis=-1)


def evaluate(records, traces, level, index):
    rng = np.random.default_rng([SEED, level, 20, index])
    effects = np.array([r["effect"] for r in records])
    p, s = effects.T
    calibration = count_binary(np.array([p - s, p + s]), rng)
    yes = calibration[:, :, 0] / SHOTS
    uncertainty = np.sqrt((yes * (1 - yes)).sum(axis=0) / SHOTS)
    best = int(np.argmax(np.abs(yes[1] - yes[0]) / np.maximum(uncertainty, 1 / SHOTS)))
    held_laws = p[best] + s[best] * PHASE_STATES[1:, 1]
    held = count_binary(held_laws, rng)
    held_p = held[:, 0] / SHOTS
    error = math.sqrt((held_p * (1 - held_p)).sum() / SHOTS)
    return dict(
        settings=traces,
        calibration_records=records,
        calibration_counts=calibration.tolist(),
        selected_setting=best,
        held_counts=held.tolist(),
        predicted_held=held_laws.tolist(),
        held_difference=float(np.diff(held_p)[0]),
        held_standard_error=error,
    )


def measure_board(level):
    board = federation(level)
    cell = angular_frame()[0]
    dark, gates = inputs(board)
    partners = [
        min(p for p in board.neighbours[g] if p // 12 == g // 12) for g in gates
    ]
    first, second = product_ball(board, np.r_[gates, partners])
    kinetic, contrast, keys = operators(board, first, second)
    gateway = seed_packet(keys, board.n, dark, gates)
    neighbor = seed_packet(keys, board.n, dark, partners)
    phased = (gateway + 1j * neighbor) / np.sqrt(2)
    for g, h in zip(gates, partners, strict=True):
        wave = (np.eye(12)[g % 12] + 1j * np.eye(12)[h % 12]) / np.sqrt(2)
        np.testing.assert_allclose(compile_wave(cell, wave), wave, atol=2e-12)
    packet = np.column_stack([gateway, phased])
    readers = [local_readers(board, g, cell) for g in gates]
    step = -0.25j * (kinetic + contrast)
    records, traces = [[], []], [[], []]
    for t in TIMES:
        packet = expm_multiply(step, packet, traceA=step.diagonal().sum())
        for k in range(2):
            for j, gate in enumerate(gates):
                branches = packet[:, 6 * k + 2 * j : 6 * k + 2 * j + 2]
                blocks = receiver_blocks(branches, first, second, gate // 12, board)
                np.testing.assert_allclose(
                    [
                        sum(matrix.trace().real for matrix in members)
                        for members in zip(*blocks)
                    ],
                    [1.0, 0.0],
                    atol=1e-13,
                )
                difference = PHASE_STATES[2, 1] - PHASE_STATES[1, 1]
                quantum = (
                    difference
                    * sum(np.abs(np.linalg.eigvalsh(y)).sum() for _, y in blocks)
                    / 2
                )
                counting = (
                    difference * sum(np.abs(np.diag(y)).sum() for _, y in blocks) / 2
                )
                assert counting <= quantum + 1e-15
                setting = dict(
                    gateway=int(gate), probe_partner=int(partners[j]), exposure=t
                )
                traces[k].append(
                    dict(
                        **setting,
                        trace_distance=float(quantum),
                        occupation_tv=float(counting),
                    )
                )
                labels, an = readers[j]
                row = np.array(
                    [
                        np.einsum("si,ij,sj->s", an.conj(), matrix, an).real
                        for matrix in blocks[1]
                    ]
                ).T
                for label, effect in zip(labels, row, strict=True):
                    records[k].append(dict(**setting, **label, effect=effect.tolist()))
    return dict(
        cells=board.cells,
        configurations=len(keys),
        state_norm_bound=2 * 3.0**DEPTH / math.factorial(DEPTH),
        studies={
            name: evaluate(records[k], traces[k], level, k)
            for k, name in enumerate(("gateway", "phased_gate_and_neighbor"))
        },
    )


def run():
    return dict(
        protocol_sha256=hashlib.sha256(
            (HERE / "PHASE_RECEIVER.md").read_bytes()
        ).hexdigest(),
        seed=SEED,
        shots_per_ensemble=SHOTS,
        boards=[measure_board(level) for level in range(3)],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = run()
    if args.output:
        with args.output.open("x") as stream:
            json.dump(data, stream, separators=(",", ":"), allow_nan=False)
            stream.write("\n")
    else:
        from model.candidate_laws.interval_echo.openwave.graviton_identity.repair_origin.coarse_response import (
            check_replay,
        )

        check_replay(json.loads((HERE / "phase_receiver.json").read_text()), data)
    for board in data["boards"]:
        for name, s in board["studies"].items():
            print(
                board["cells"],
                name,
                max(r["trace_distance"] for r in s["settings"]),
                max(r["occupation_tv"] for r in s["settings"]),
                s["held_difference"],
                s["held_standard_error"],
            )
