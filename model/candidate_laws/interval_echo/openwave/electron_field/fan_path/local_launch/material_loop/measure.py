"""Replay material-counter, returning-word and stationary-current records at 1045."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from ....replay import compare
from ...fixture import FanPath, archived_fixture
from ...measure import conjugacy
from ..face_current.measure import pair_amplitudes
from ..measure import PAULI
from .triangle import (
    HUB,
    METER_SIZE,
    exact_leading_determinant,
    geometry,
    get_counters,
    side_analyzers,
    stationary_meter,
)


HERE = Path(__file__).resolve().parent
LABEL = [1045, 0, 334, 365]
PROP_ERROR = 3e-10
TIMES = [27.0625, 30.8125]


def partial_transpose(operator):
    return operator.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)


def moments(rho):
    return np.array(
        [[np.trace(rho @ np.kron(a, b)).real for b in PAULI] for a in PAULI]
    )


def source_records(path, archive):
    case = next(c for c in archive["cases"] if c["fan_count"] == 4)
    arm = next(a for a in case["stock_arms"] if a["preparation"] == "clock_bypassed")
    out = []
    for time in TIMES:
        flight = next(
            f
            for f in arm["flights"]
            if f["appointment"] == time and f["material_bridge_rate"] == 4
        )
        gate = next(
            p for p in flight["pair_readings"] if p["port_face_and_words"] == LABEL
        )
        wave = expm_multiply(-1j * time * path.h, path.source())
        samples, cell, hub = pair_amplitudes(path, wave, LABEL)
        flat = samples.reshape(-1, 4)
        rho = flat.T @ flat.conj() / 2
        t = moments(rho)
        np.testing.assert_allclose(
            t, gate["receiver_reference_pauli_moments"], atol=1e-13, rtol=0
        )
        out.append(
            (
                rho,
                t,
                {
                    "source_flight": time,
                    "rejected_before_terminal_fan": 1 - cell,
                    "terminal_fan_away_from_hub": cell - hub,
                    "wrong_hub_word": hub - t[0, 0],
                    "passed_on_original_launch": t[0, 0],
                    "input_pauli_error_upper": flight[
                        "observable_absolute_error_upper_per_pair"
                    ],
                    "input_pt_minimum_interval": gate[
                        "partial_transpose_minimum_interval"
                    ],
                },
            )
        )
    return out


def matched_axis(t):
    normalized = t / t[0, 0]
    cov = normalized[1:, 1:] - np.outer(normalized[1:, 0], normalized[0, 1:])
    n = np.cross(cov[0], cov[2])
    n /= np.linalg.norm(n)
    return n if cov[1] @ n >= 0 else -n


def cohorts(t, n):
    bank = {}
    for name, weight, unit in [
        ("plus", 0.5, n),
        ("minus", 0.5, -n),
        ("unread", 1, n * 0),
    ]:
        coeff = weight * (t[:, 0] + t[:, 1:] @ unit)
        rho = np.einsum("a,aij->ij", coeff / 2, PAULI)
        assert np.linalg.eigvalsh(rho)[0] > 0
        bank[name] = (rho, coeff[1:] / coeff[0])
    np.testing.assert_allclose(
        bank["plus"][0] + bank["minus"][0], bank["unread"][0], atol=1e-15
    )
    return bank


def complex_entries(matrix):
    return np.stack((matrix.real, matrix.imag), axis=-1).tolist()


def circuit_conditional(k, rho, gate):
    joint_map = np.kron(k, np.eye(2))
    result = joint_map @ rho @ joint_map.conj().T
    p = float(result.trace().real)
    least = float(np.linalg.svd(k, compute_uv=False)[-1])
    minimum = float(np.linalg.eigvalsh(partial_transpose(result))[0])
    hi_input = gate["input_pt_minimum_interval"][1]
    # Congruence under conditioning gives a robust negative witness even
    # when the returning field is close to a single word.
    upper = hi_input * max(0, least - PROP_ERROR) ** 2 if hi_input < 0 else None
    if upper is not None:
        assert minimum <= upper + 1e-15
    return {
        "word_map_complex_pairs": complex_entries(k),
        "per_original_launch": p,
        "least_singular_value": least,
        "pt_minimum_per_original": minimum,
        "pt_minimum_congruence_upper": upper,
        "return_word_source_pauli_moments_per_original": moments(result).tolist(),
    }


def counter_read(k, rho, gate, bank):
    result = {}
    for name, (local, _) in bank.items():
        table = np.einsum("mvi,ij,mvj->mv", k, local, k.conj()).real
        assert table.min() > -1e-15
        p = local.trace().real
        np.testing.assert_allclose(table.sum(), p, atol=2e-14, rtol=0)
        hub = table[:, HUB].sum(axis=1)
        result[name] = {
            "face_and_reference_per_original": p,
            "displaced_to_rim_per_original": table[:, 1:3].sum(),
            "zero_circuit_hub_per_original": hub[0],
            "other_winding_hub_per_original": hub[2:-1].sum(),
            "positive_hub_per_original": hub[1],
            "negative_hub_per_original": hub[-1],
            "joint_counter_port_word_per_original": table.tolist(),
            "oriented_single_circuit_score_given_face_and_reference": (hub[1] - hub[-1])
            / p,
        }
    score_key = "oriented_single_circuit_score_given_face_and_reference"
    diff = result["plus"][score_key] - result["minus"][score_key]
    effect = k[1, HUB].conj().T @ k[1, HUB] - k[-1, HUB].conj().T @ k[-1, HUB]
    coeff = np.array([np.trace(v @ effect).real / 2 for v in PAULI[1:]])
    delta_bloch = bank["plus"][1] - bank["minus"][1]
    np.testing.assert_allclose(diff, delta_bloch @ coeff, atol=2e-14, rtol=0)
    return {
        "cohorts": result,
        "matched_score_difference": diff,
        "score_xyz_coefficients": coeff.tolist(),
        "positive_return": circuit_conditional(k[1, HUB], rho, gate),
        "negative_return": circuit_conditional(k[-1, HUB], rho, gate),
    }


def read_stationary(meter, bank, source):
    first, both, j, free_wait, native_lines = meter
    analyzers = side_analyzers()
    out = {}
    for name, (local, _) in bank.items():
        after = both @ local @ both.conj().T
        p0, p1, p2 = (
            np.trace(a).real for a in (local, first @ local @ first.conj().T, after)
        )
        current = float(np.trace(j @ after).real)
        side_counts = []
        for i, dest, turn in analyzers:
            for additional_wait, evolved in (
                (0, after),
                (1, free_wait @ after @ free_wait.conj().T),
            ):
                final = np.diag(turn @ evolved @ turn.conj().T).real
                src, dst = final[[i, dest]]
                np.testing.assert_allclose(4 * (dst - src), current, atol=2e-14, rtol=0)
                side_counts.append(
                    {
                        "chain_edge": [i, dest],
                        "later_wait_in_pi": additional_wait,
                        "source_destination_other_per_original": [
                            src,
                            dst,
                            p2 - src - dst,
                        ],
                    }
                )
        native_counts = []
        for line in native_lines:
            selected = line @ local @ line.conj().T
            equilibrium = []
            for i, dest, turn in analyzers:
                final = np.diag(turn @ selected @ turn.conj().T).real
                equilibrium.append(final[[i, dest]].tolist())
                np.testing.assert_allclose(final[dest], final[i], atol=1e-14, rtol=0)
            native_counts.append(
                {
                    "line_per_original": selected.trace().real,
                    "side_detector_pairs": equilibrium,
                }
            )
        p_error = 4 * source["input_pauli_error_upper"] + PROP_ERROR
        # j has norm 4. For this linear expectation use the local trace distance.
        lo, hi = (current + s * 4 * p_error for s in (-1, 1))
        conditional = [
            min(a / b for a in (lo, hi) for b in (p2 - p_error, p2 + p_error)),
            max(a / b for a in (lo, hi) for b in (p2 - p_error, p2 + p_error)),
        ]
        out[name] = {
            "first_refused_after_face_per_original": p0 - p1,
            "second_refused_after_face_per_original": p1 - p2,
            "both_pass_per_original": p2,
            "spatial_current_per_original": current,
            "spatial_current_given_pass": current / p2,
            "spatial_current_given_pass_interval": conditional,
            "three_side_counts": side_counts,
            "energy_resolved_native_control": native_counts,
        }
    return out


def run():
    reads = HERE.parent / "word_reads.json"
    with archived_fixture() as archive:
        path = FanPath(archive, 4, 1, asymmetric=True)
        vacancy = FanPath(archive, 4, -1, asymmetric=True)
        charts = {"excess": geometry(path), "vacancy": geometry(vacancy)}
        gauges = conjugacy(path, vacancy)
        sources = source_records(path, json.loads(reads.read_text()))
    n = matched_axis(sources[-1][1])
    counters = get_counters()
    comparison, meter = stationary_meter()
    flights = []
    for rho, t, gate in sources:
        bank = cohorts(t, n)
        pmin = min(float(r.trace().real) for r, _ in bank.values())
        pad = gate["input_pauli_error_upper"]
        # Each steered Pauli expectation and its success probability have
        # error <= (1+sqrt(3))*pad/2. Also include propagator error in scores.
        ratio_error = 8 * pad / (pmin - 4 * pad) + 8 * PROP_ERROR
        if gate["source_flight"] == TIMES[-1]:
            diff = bank["plus"][1] - bank["minus"][1]
            np.testing.assert_allclose(diff[[0, 2]], 0, atol=2e-14)
            assert diff[1] > 0.1
        cases = []
        for (g, eighths), k in counters.items():
            case = counter_read(k, rho, gate, bank)
            if g == 0:
                for m in (-1, 1):
                    assert np.linalg.svd(k[m, HUB], compute_uv=False)[-1] < 2e-14
                assert abs(case["score_xyz_coefficients"][1]) < 2e-14
            cases.append({"ring_rate": g, "wait_in_pi_over_8": eighths, **case})
        flights.append(
            {
                **gate,
                "reference_axis_xyz": n.tolist(),
                "source_bloch_vectors": {
                    name: bloch.tolist() for name, (_, bloch) in bank.items()
                },
                "absolute_guard_each_conditional_score_or_bloch_component": ratio_error,
                "counter_readings": cases,
                "tuned_stationary_filter": read_stationary(meter, bank, gate),
            }
        )
    largest = flights[-1]
    scores = [
        abs(c["matched_score_difference"])
        for c in largest["counter_readings"]
        if c["ring_rate"] == 0.5 and c["wait_in_pi_over_8"]
    ]
    assert (
        max(scores)
        > 4 * largest["absolute_guard_each_conditional_score_or_bloch_component"]
    )
    assert (
        largest["tuned_stationary_filter"]["unread"][
            "spatial_current_given_pass_interval"
        ][0]
        * (
            largest["tuned_stationary_filter"]["unread"][
                "spatial_current_given_pass_interval"
            ][1]
        )
        > 0
    )
    return {
        "design_sha256": sha256((HERE / "DESIGN.md").read_bytes()).hexdigest(),
        "arrival_archive_sha256": sha256(reads.read_bytes()).hexdigest(),
        "register_residue_order": list(range(METER_SIZE)),
        "microscopic_charts": charts,
        "conjugate_source_receiver_gauges": gauges,
        "counter_transfer_operator_error_upper": PROP_ERROR,
        "single_return_alias_amplitude_upper": float(
            2 * Fraction(12 * 22, 7) ** 192 / factorial(192)
        ),
        "determinant_leading_coefficient_at_wait_power_8": exact_leading_determinant(),
        "stationary_comparison": comparison,
        "flights": flights,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write a new receipt")
    args = parser.parse_args()
    data = json.loads(json.dumps(run(), allow_nan=False))
    if args.output:
        with args.output.open("x") as f:
            f.write(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        compare(data, json.loads((HERE / "records.json").read_text()), "material_loop")
        print(
            "Verified both source flights, 24 counter banks, both signed charts and the stationary comparison."
        )


if __name__ == "__main__":
    main()
