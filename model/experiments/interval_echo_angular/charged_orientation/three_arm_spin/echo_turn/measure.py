"""Spatial turns of the remounted, physically delivered three-head code."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from model.executable_model.symmetry import rotations
from ..apparatus import native_graph
from ..measure import BASE, ONE, X, Y, CY, build, pt
from ..rotation_domain import rotation_order
from .material import chart, compose, dressed, encode, rotate_program

HERE = Path(__file__).resolve().parent
IDENTITY = tuple(range(12))
SHOTS = 2**38


def sequence(turns, phase_first=False):
    frame = IDENTITY
    operation = np.eye(8, dtype=complex)
    for i, r in enumerate(turns):
        frame, step = dressed(frame, r, phase_first and i == 0)
        operation = step @ operation
    return frame, operation


def closed(turns, phase_first=False):
    frame, op = sequence(turns, phase_first)
    assert frame == IDENTITY
    return op


def turnaround():
    lab = chart()
    moves = rotations()
    index = {r: i for i, r in enumerate(moves)}

    def gateway(r):
        return {lab["base"] + r[p - lab["base"]] for p in lab["gateways"]}

    rigid = [i for i, r in enumerate(moves) if gateway(r) == lab["gateways"]]
    assert len(rigid) == 3 and [rotation_order(moves[i]) for i in rigid] == [1, 3, 3]
    half = [r for r in moves if rotation_order(r) == 2]
    pairs = [(x, y) for x, y in combinations(half, 2) if compose(x, y) == compose(y, x)]
    assert len(pairs) == 15
    errors = []
    for x, y in pairs:
        common = compose(x, y)
        a = closed([x, y, common])
        b = closed([y, x, common])
        errors.append(float(np.linalg.norm(a.conj().T @ b - np.eye(8), 2)))
    all_returns = [
        float(np.linalg.norm(closed([r] * rotation_order(r)) - np.eye(8), 2))
        for r in moves
    ]
    cycle = moves[rigid[1]]
    positive = closed([cycle] * 3)
    negative = closed([cycle] * 3, True)
    material_phase = np.diag([(-1) ** m.bit_count() for m in range(8)])
    np.testing.assert_allclose(positive, np.eye(8), atol=1e-11)
    np.testing.assert_allclose(negative, material_phase, atol=1e-11)
    # Independent mathematical spin-half control with phase-sensitive readout.
    hx, hy = -1j * X, -1j * Y
    commutator = (hx @ hy).conj().T @ (hy @ hx)
    np.testing.assert_array_equal(commutator, -np.eye(2))
    return (
        dict(
            proper_rotations=len(moves),
            gateway_preserving=rigid,
            addressed_internal_full_turns=len(all_returns),
            max_full_turn_operator_error=max(all_returns),
            perpendicular_pairs=[[index[x], index[y]] for x, y in pairs],
            max_half_turn_order_operator_error=max(errors),
            full_boundary_cycle=list(cycle),
            blank_swaps_per_boundary_step=len(rotate_program(cycle)),
            introduced_material_phase_signs=material_phase.diagonal().tolist(),
            spinor_control_dark_probability=1,
        ),
        positive,
        negative,
    )


def detected(slab, operation):
    rotated = np.einsum("jr,rsi->jsi", operation, slab)
    return np.einsum("ro,rsi->soi", BASE, rotated[ONE])


def bank(slab, operation, rng, admission=1.0):
    bright = detected(slab, (np.eye(8) + operation) / 2)
    dark = detected(slab, (np.eye(8) - operation) / 2)
    bare = float(np.linalg.norm(slab) ** 2 / 2)
    b, d = pt(bright)[0] / 9, pt(dark)[0] / 9
    bins = [b, d, bare / 3 - b - d, (admission - bare) / 3]
    if admission < 1:
        bins.append((1 - admission) / 3)
    raw = np.array([*bins, 2 / 3])
    assert np.all(raw >= 0)
    counts = rng.multinomial(SHOTS, raw / raw.sum())
    n = int(counts[0] + counts[1])
    k = int(counts[1])
    # Exact conditional interval for dark fraction among accepted launches;
    # simultaneous 99 percent over the two independent phase settings.
    tail = 0.01 / 4
    lo = beta.ppf(tail, k, n - k + 1) if k else 0.0
    hi = beta.isf(tail, k + 1, n - k) if k < n else 1.0
    return dict(
        starts=SHOTS,
        counts=counts.tolist(),
        probabilities=raw.tolist(),
        dark_fraction_99=[float(lo), float(hi)],
    )


def filtered_delivery(positive, negative):
    # Recompute the entire remote register, including every unread source word.
    # Take the finite-clock instrument and aim from the recorded experiment.
    earlier = json.loads((HERE.parent / "results.json").read_text())
    clock = earlier["clock"]
    ev, _, beta0, q, v, _, perms = build()
    n = clock["number_delays"]
    angles = (ev - clock["center"]) * clock["tick"] / (2 * np.pi)
    f = np.exp(-1j * np.pi * (n - 1) * angles) * np.sinc(n * angles) / np.sinc(angles)
    weighted = v @ (beta0 * f)
    admit = float(np.linalg.norm(weighted) ** 2)
    c = perms[CY]
    rows = np.arange(64, 128)
    first = q[rows] @ weighted
    second = ((q[c[rows]] - q[c[c[rows]]]) @ weighted) / np.sqrt(3)
    check = BASE[:, 1] @ ((q[c[ONE]] - q[c[c[ONE]]]) @ v) / np.sqrt(3)
    if np.dot(beta0, check) < 0:
        second = -second
    slab = np.stack((first, second), axis=-1).reshape(8, 8, 2)
    incoming = pt(detected(slab, np.eye(8)))
    outgoing = pt(detected(slab, positive))
    np.testing.assert_allclose(
        outgoing[0] / 9, earlier["filtered"]["ideal_total_yield"], rtol=1e-8, atol=1e-18
    )
    np.testing.assert_allclose(
        outgoing[1] / 9,
        earlier["filtered"]["original_PT_minimum"],
        rtol=1e-8,
        atol=1e-18,
    )
    np.testing.assert_allclose(incoming[2], outgoing[2], atol=1e-17)
    rng = np.random.default_rng(608434)
    return dict(
        energy=earlier["selected_energy"],
        clock=clock,
        admission=admit,
        decoded_yield_per_original=outgoing[0] / 9,
        PT_minimum_per_original=outgoing[1] / 9,
        conditional_PT=outgoing[1] / outgoing[0],
        Choi_change_norm=float(np.linalg.norm(incoming[2] - outgoing[2])),
        simulation_seed=608434,
        bin_order=[
            "bright",
            "dark",
            "remote_code_refusal",
            "no_remote_charge_after_clock",
            "clock_refusal",
            "source_refusal",
        ],
        boundary_full_turn=bank(slab, positive, rng, admit),
        boundary_full_turn_with_material_private_phase=bank(slab, negative, rng, admit),
    )


def run():
    lab = chart()
    result, positive, negative = turnaround()
    h = native_graph()
    source = np.zeros((h.shape[0], 2))
    source[ONE] = BASE
    evolution = expm_multiply((-32j) * h, source, traceA=(-32j) * h.diagonal().sum())
    slab = evolution[64:128].reshape(8, 8, 2)
    incoming = pt(detected(slab, np.eye(8)))
    outgoing = pt(detected(slab, positive))
    np.testing.assert_allclose(incoming[2], outgoing[2], atol=2e-16)
    rng = np.random.default_rng(608433)
    return dict(
        material=dict(
            terminal=lab["tip"],
            direction_heads=lab["marks"],
            staging=lab["storage"],
            probe_vacancies=lab["probes"],
            probe_compensators=lab["stock"],
            original_compensator=lab["donor"],
            phase_reference=lab["reference"],
            service_shutter=lab["service"],
            base_phase_route_swaps=encode(IDENTITY)[1],
            original_private_closures=len(lab["private"]),
            internal_silent_closures=len(lab["quiet"]),
        ),
        turns=result,
        free_delivery=dict(
            time=32,
            source_and_decoder_apertures=[1 / 3, 1 / 3],
            terminal_given_source=float(np.linalg.norm(slab) ** 2 / 2),
            decoded_yield_per_original=outgoing[0] / 9,
            PT_minimum_per_original=outgoing[1] / 9,
            conditional_PT=outgoing[1] / outgoing[0],
            Choi_change_norm=float(np.linalg.norm(incoming[2] - outgoing[2])),
            simulation_seed=608433,
            bin_order=[
                "bright",
                "dark",
                "remote_code_refusal",
                "no_remote_charge",
                "source_refusal",
            ],
            boundary_full_turn=bank(slab, positive, rng),
            boundary_full_turn_with_material_private_phase=bank(slab, negative, rng),
        ),
        filtered_delivery=filtered_delivery(positive, negative),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    path = parser.parse_args().output
    path.write_text(json.dumps(run(), indent=2) + "\n")
