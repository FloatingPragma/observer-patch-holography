"""Certify a finite coarse and coherently erased energy admission instrument."""

import argparse
from fractions import Fraction
import json
from math import ceil, floor
from pathlib import Path
from runpy import run_path
import sys

import numpy as np
from scipy.sparse.linalg import eigsh

HERE = Path(__file__).resolve().parent
FIELD = HERE.parent
BITS = 18
LENGTH = 1 << BITS
REPEATS = 41
TAU = Fraction(1, 64)
WINDOW = (Fraction(-1), Fraction(9, 2))
LOWER_CUT = Fraction(17, 4)
UPPER_CUT = Fraction(145, 32)


def pi_interval():
    """Enclose pi using Machin's identity and rational alternating series."""
    intervals = []
    for x in (Fraction(1, 5), Fraction(1, 239)):
        upper = sum((-1) ** k * x ** (2 * k + 1) / (2 * k + 1) for k in range(9))
        lower = upper - x**19 / 19
        intervals.append((lower, upper))
    (a, b), (c, d) = intervals
    return 16 * a - 4 * d, 16 * b - 4 * c


def circuit():
    low_pi, high_pi = pi_interval()
    # An output bin j represents estimate -2 pi signed(j)/(tau N).
    half_bin = TAU * LENGTH / 2
    minimum = ceil(WINDOW[0] * half_bin / low_pi)
    maximum = floor(WINDOW[1] * half_bin / low_pi)
    assert minimum == ceil(WINDOW[0] * half_bin / high_pi)
    assert maximum == floor(WINDOW[1] * half_bin / high_pi)
    margin = min(1, WINDOW[1] - LOWER_CUT, UPPER_CUT - WINDOW[1])
    bins_away = floor(LENGTH * TAU * margin / (2 * high_pi))
    assert bins_away >= 2
    one_read_error = Fraction(1, 2 * (bins_away - 1))
    majority = (REPEATS + 1) // 2
    y_error = 2**REPEATS * one_read_error**majority
    votes = maximum - minimum + 1
    counter_bits = REPEATS.bit_length()
    return {
        "phase_bits_per_register": BITS,
        "grid_points": LENGTH,
        "controlled_base_exposure": str(TAU),
        "energy_window": list(map(str, WINDOW)),
        "accepted_negative_of_signed_residue_from_to_inclusive": [minimum, maximum],
        "actual_unsigned_accepted_bin_intervals_inclusive": [
            [0, -minimum],
            [LENGTH - maximum, LENGTH - 1],
        ],
        "guaranteed_energy_margin_to_incorrect_estimate": str(margin),
        "minimum_whole_clock_bins_to_incorrect_estimate": bins_away,
        "one_read_misclassification_probability_upper": str(one_read_error),
        "coherent_repeats": REPEATS,
        "majority_votes_required": majority,
        "operator_error_upper_exact": str(y_error),
        "operator_error_upper": upper(y_error),
        "resources": {
            "phase_qubits": REPEATS * BITS,
            "vote_and_counter_and_two_binary_record_qubits": REPEATS + counter_bits + 2,
            "controlled_forward_and_inverse_power_calls": 2 * REPEATS * BITS,
            "longest_single_controlled_evolution": str(TAU * (LENGTH // 2)),
            "sum_of_absolute_controlled_evolution_exposures": str(
                2 * REPEATS * (LENGTH - 1) * TAU
            ),
            "phase_circuit_hadamards": 2 * REPEATS * 2 * BITS,
            "phase_circuit_controlled_phase_gates": 2
            * REPEATS
            * BITS
            * (BITS - 1)
            // 2,
            "phase_circuit_bit_reversal_swaps": 2 * REPEATS * (BITS // 2),
            "literal_vote_comparator_multi_controlled_X_upper": 2 * REPEATS * votes,
            "literal_vote_comparator_single_X_upper": 4 * REPEATS * votes * BITS,
            "majority_compute_and_uncompute_controlled_X_upper": 2
            * REPEATS
            * counter_bits
            + 2**counter_bits
            - majority,
            "clock_reset_binary_read_controlled_X": 1,
        },
    }, y_error


def upper(x):
    return float(np.nextafter(float(x), np.inf))


def lower(x):
    return float(np.nextafter(float(x), -np.inf))


def small_circuit_check():
    """Check yes/no and reset amplitudes with three explicit two-bit registers."""
    k = np.arange(4)
    inverse_fourier = np.exp(-2j * np.pi * np.outer(k, k) / 4) / 2
    four_hadamard = np.kron([[1, 1], [1, -1]], [[1, 1], [1, -1]]) / 2
    vote = np.indices((4,) * 3).reshape(3, -1)
    majority = np.count_nonzero(vote == 0, axis=0) >= 2
    results = []
    for phase in (0.0, 0.125, 0.25):
        qpe = inverse_fourier @ (
            np.exp(2j * np.pi * phase * k[:, None]) * four_hadamard
        )
        repeated = np.kron(qpe, np.kron(qpe, qpe))
        phase_read_state = repeated[:, 0]
        yes = repeated.conj().T @ (majority * phase_read_state)
        no = repeated.conj().T @ (~majority * phase_read_state)
        p = abs(qpe[0, 0]) ** 2
        f = 3 * p**2 * (1 - p) + p**3
        np.testing.assert_allclose([yes[0], no[0]], [f, 1 - f], atol=3e-15, rtol=0)
        failed_reset = float(np.sum(abs(yes[1:]) ** 2 + abs(no[1:]) ** 2))
        np.testing.assert_allclose(failed_reset, 2 * f * (1 - f), atol=3e-15, rtol=0)
        results.append(
            {
                "phase": phase,
                "yes_and_reset_probability": float(abs(yes[0]) ** 2),
                "no_and_reset_probability": float(abs(no[0]) ** 2),
                "reset_refusal_probability": failed_reset,
            }
        )
    return results


def spectral_domain(h, verified_inertia):
    doubled = (2 * h).astype(np.int64)
    assert (doubled != 2 * h).nnz == (doubled != doubled.T).nnz == 0
    off_diagonal = doubled.tocoo()
    assert not np.any(off_diagonal.data[off_diagonal.row != off_diagonal.col] > 0)
    assert np.min(doubled.sum(axis=1)) >= 0
    # The sign condition and nonnegative row sums make H positive semidefinite.
    bound = Fraction(int(np.max(abs(doubled).sum(axis=1))), 2)
    low_pi, _ = pi_interval()
    assert TAU * bound < low_pi / 2
    counts = [verified_inertia(h, cut) for cut in (LOWER_CUT, UPPER_CUT)]
    assert all(record["count_below"] == 2 for record in counts)
    return str(bound), counts


def finite_reads(certificate, y_error, sqrt_bounds):
    source = Fraction(certificate["source_minimum_singular_value_lower"])
    receiver = Fraction(certificate["vacuum_and_receiver_minimum_singular_value_lower"])
    joint = source * receiver
    assert 1000 * y_error < min(source, joint)
    local_negative_margin = Fraction(
        certificate[
            "all_waits_receiver_only_negative_eigenvalue_magnitude_per_launch_lower"
        ]
    )
    assert local_negative_margin > 2 * y_error
    admitted_error = 2 * y_error / source
    assert receiver > admitted_error
    separated = []
    for number in certificate["all_states_separated_probability_upper_by_distance"]:
        ideal = Fraction(number)
        amplitude = sqrt_bounds(ideal)[1]
        separated.append(upper(min(1, (amplitude + admitted_error) ** 2)))
    return {
        "coarse_yes_and_reset_probability_lower_all_source_inputs": lower(
            (source - y_error) ** 2
        ),
        "joint_all_other_fans_vacuum_and_receiving_doublet_probability_lower_all_inputs_and_waits": lower(
            (joint - y_error) ** 2
        ),
        "joint_receiver_probability_lower_conditioned_on_admission_all_inputs_and_waits": lower(
            (receiver - admitted_error) ** 2
        ),
        "joint_amplitude_error_fraction_of_ideal_worst_case": upper(y_error / joint),
        "admitted_normalized_vector_error_upper_all_source_inputs_and_waits": upper(
            admitted_error
        ),
        "admitted_separation_probability_upper_by_distance_all_source_inputs_and_waits": separated,
        "receiver_only_unnormalized_choi_trace_norm_error_upper": upper(2 * y_error),
        "receiver_only_unnormalized_choi_partial_transpose_minimum_eigenvalue_shift_upper": upper(
            2 * y_error
        ),
        "receiver_only_unnormalized_choi_negative_eigenvalue_magnitude_lower_all_waits": lower(
            local_negative_margin - 2 * y_error
        ),
    }


def run():
    sys.path.insert(0, str(FIELD.parents[4]))
    from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.channel_bound import (
        certify_channel,
        sqrt_bounds,
    )
    from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.fixture import (
        FanPath,
        archived_fixture,
    )
    from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.measure import (
        conjugacy,
    )
    from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.spectral import (
        separated_reference,
        verified_inertia,
    )

    design, error = circuit()
    rows = []
    with archived_fixture() as archive:
        for n in (3, 4):
            path = FanPath(archive, n, 1, asymmetric=True)
            negative = FanPath(archive, n, -1, asymmetric=True)
            conjugation = conjugacy(path, negative)
            h = path.h
            assert path.board.cells == 320
            bound, counts = spectral_domain(h, verified_inertia)
            _, bracket, parked = separated_reference(path)
            assert WINDOW[1] < bracket[0] < bracket[1] < UPPER_CUT
            rng = np.random.default_rng(923)
            energy, vectors = eigsh(
                h, k=4, which="SA", v0=rng.standard_normal(h.shape[0]), tol=1e-13
            )
            floor_energy = min(
                Fraction(block["charged_ground_dyadic_bracket"][0]) for block in parked
            )
            channel = certify_channel(
                path, h, energy[:2], vectors[:, :2], LOWER_CUT, floor_energy
            )
            rows.append(
                {
                    "board_cells": path.board.cells,
                    "fan_count": n,
                    "configurations_per_sign": h.shape[0],
                    "physical_hubs": [fan[0] for fan in path.cells],
                    "bridges": [list(link) for link in path.bridges],
                    "operator_norm_bound_for_phase_unwrapping": bound,
                    "charge_conjugacy_source_and_receiver_gauges": conjugation,
                    "absolute_energy_brackets_for_both_signs": counts,
                    "parked_separated_threshold_bracket": list(map(str, bracket)),
                    "ideal_transfer_and_separation_certificate": {
                        key: channel[key]
                        for key in (
                            "dyadic_energy_centers",
                            "eigenvalue_error_ceilings",
                            "eigenvector_frame_error_upper",
                            "source_minimum_singular_value_lower",
                            "vacuum_and_receiver_minimum_singular_value_lower",
                            "all_waits_joint_probability_lower_for_every_input",
                            "all_waits_receiver_only_negative_eigenvalue_magnitude_per_launch_lower",
                            "all_states_separated_probability_upper_by_distance",
                        )
                    },
                    "finite_coarse_filter_certificate": finite_reads(
                        channel, error, sqrt_bounds
                    ),
                }
            )
    return {
        "finite_instrument": design,
        "small_unitary_clock_check": small_circuit_check(),
        "paths": rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = run()
    if args.output:
        with args.output.open("x") as f:
            f.write(json.dumps(record, indent=2, allow_nan=False) + "\n")
    else:
        compare = run_path(str(FIELD / "replay.py"))["compare"]
        compare(
            record,
            json.loads((HERE / "coherent_filter.json").read_text()),
            "coherent_filter",
        )
        print("Verified the finite coherent filter and receiver bounds.")


if __name__ == "__main__":
    main()
