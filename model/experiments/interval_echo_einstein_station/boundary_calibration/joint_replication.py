"""Fresh held occupation and circular counts on a fixed joint aperture."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .measure import (
    HERE,
    PHASE_STATES,
    SEED,
    compile_wave,
    federation,
    inputs,
    operators,
    product_ball,
    seed_packet,
)
from .two_gateways import joint_blocks, local_carrier, outside_route

REPEATS = 2**43


def contrast(rows, selected):
    counts = np.asarray(rows)[:, selected].sum(axis=-1)
    probabilities = counts / REPEATS
    return dict(
        yes_counts=counts.tolist(),
        difference=float(np.diff(probabilities)[0]),
        standard_error=float(
            np.sqrt(np.sum(probabilities * (1 - probabilities)) / REPEATS)
        ),
    )


def counts(mean, variation, parameters, rng):
    laws = mean + np.asarray(parameters)[:, None] * variation
    np.testing.assert_allclose(laws.sum(axis=-1), 1.0, rtol=0, atol=2e-13)
    assert laws.min() > -1e-13
    laws = np.maximum(laws, 0.0)
    return np.array([rng.multinomial(REPEATS, row / row.sum()) for row in laws])


def measure(level):
    board = federation(level)
    dark, gates = inputs(board)
    g, h = map(int, gates[:2])
    visible = np.sort(
        np.r_[
            np.arange(g - g % 12, g - g % 12 + 12),
            np.arange(h - h % 12, h - h % 12 + 12),
        ]
    )
    route = outside_route(board, g, h)
    aux = [p for p in route if p not in visible]
    sites = np.sort(np.unique(np.r_[visible, route])).astype(int)
    reader = np.zeros(len(sites), complex)
    reader[np.searchsorted(sites, [g, h])] = [1 / np.sqrt(2), 1j / np.sqrt(2)]
    np.testing.assert_allclose(
        compile_wave(local_carrier(board, sites), reader), reader, atol=1e-12
    )
    first, second = product_ball(board, gates)
    kinetic, agreement, keys = operators(board, first, second)
    packet = seed_packet(keys, board.n, dark, gates)[:, :4]
    gen = -0.5j * (kinetic + agreement)
    packet = expm_multiply(gen, packet, traceA=gen.diagonal().sum())
    a, b = packet[:, :2], packet[:, 2:]
    arms = [
        joint_blocks(v, first, second, visible, aux, board)
        for v in (a, b, (a + 1j * b) / np.sqrt(2))
    ]
    mixed = [
        tuple((m + n) / 2 for m, n in zip(l, r, strict=True))
        for l, r in zip(arms[0][0], arms[1][0], strict=True)
    ]
    mixed_route = tuple(
        (m + n) / 2 for m, n in zip(arms[0][1], arms[1][1], strict=True)
    )
    result = {}
    for k, (name, (blocks, routed)) in enumerate(
        zip(("incoherent", "coherent"), [(mixed, mixed_route), arms[2]], strict=True)
    ):
        rng = np.random.default_rng([SEED, level, k, 63])
        mean, variation = (
            np.concatenate([np.diag(m).real for m in members])
            for members in zip(*blocks)
        )
        calibration = counts(mean, variation, [-1.0, 1.0], rng)
        cal = calibration / REPEATS
        threshold = 3 * np.sqrt(np.sum(cal * (1 - cal), axis=0) / REPEATS)
        keep = cal[1] - cal[0] > threshold
        held = counts(mean, variation, PHASE_STATES[1:, 1], rng)
        indices = np.searchsorted(visible, [g, h])
        rail_wave = np.zeros(len(visible), complex)
        rail_wave[indices] = [1 / np.sqrt(2), 1j / np.sqrt(2)]
        center, slope = [float((rail_wave.conj() @ m @ rail_wave).real) for m in routed]
        binary = counts(
            np.array([center, 1 - center]),
            np.array([slope, -slope]),
            PHASE_STATES[1:, 1],
            rng,
        )
        result[name] = dict(
            calibration_counts=calibration.tolist(),
            held_occupations=held.tolist(),
            selected_occupations=np.flatnonzero(keep).tolist(),
            occupation_contrast=contrast(held, keep),
            occupation_predicted_difference=float(0.25 * variation[keep].sum()),
            held_routed_reader_counts=binary.tolist(),
            routed_contrast=contrast(binary, [0]),
            routed_predicted_difference=0.25 * slope,
        )
    return dict(cells=board.cells, gateways=[g, h], auxiliary_ports=aux, studies=result)


def run():
    return dict(
        protocol_sha256=hashlib.sha256(
            (HERE / "JOINT_REPLICATION.md").read_bytes()
        ).hexdigest(),
        seed=SEED,
        shots_per_ensemble=REPEATS,
        boards=[measure(level) for level in range(3)],
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

        check_replay(json.loads((HERE / "joint_replication.json").read_text()), data)
    for board in data["boards"]:
        for name, s in board["studies"].items():
            print(board["cells"], name, s["occupation_contrast"], s["routed_contrast"])
