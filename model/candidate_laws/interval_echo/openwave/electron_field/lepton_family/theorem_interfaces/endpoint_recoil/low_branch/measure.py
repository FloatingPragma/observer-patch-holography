"""Enclose low interacting spectral support and a finite-clock receiver window."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import comb, isqrt
import json
from pathlib import Path

from .. import energy_bridge, measure as recoil


HERE = Path(__file__).resolve().parent
GROUND = "model/particle_primary_occupation_bound/results.json"
FOURIER = "model/particle_primary_joint_spectroscopy/PROOF.md"
FLOOR = -F(1455, 2)
GROUND_UPPER = -F(2895, 4)


def sqrt_upper(value):
    units = 2**40
    return F(isqrt(value.numerator * units**2 // value.denominator) + 1, units)


def low_witness(center, half_gap, variance_cap, turn, guard=2):
    """A negative quadratic on the low interval, bounded on the high tail."""
    floor = FLOOR - center
    k = turn / (4 * (2 * turn**2 - half_gap**2))

    def r(x):
        return x - turn - k * (x**2 - half_gap**2)

    assert turn > half_gap and 1 - 2 * k * (turn + guard) > 0
    assert r(turn + guard) < 0
    normalizer = -r(floor)
    assert normalizer > 0
    # Maximum of (x-turn)/(x²-d²), x>=turn, is
    # 1/(2(turn+sqrt(turn²-d²))). Use sqrt(turn²-d²)>=turn-d²/turn.
    tail = k * variance_cap
    assert tail < F(1, 2)
    lipschitz = (1 - 2 * k * floor - r(turn) / guard) / normalizer
    return {
        "window": [FLOOR, center + turn + guard],
        "first_edge_of_high_tail": center + turn,
        "linear_top_guard_width": guard,
        "polynomial_coefficient": k,
        "maximum_weighted_tail_ratio": k,
        "tail_contribution_upper": tail,
        "polynomial_normalizer": normalizer,
        "ideal_receiver_amplitude_lower": (F(1, 2) - tail) / normalizer,
        "filter_Lipschitz_per_energy_upper": lipschitz,
    }


def excited_witness(center, half_gap, variance_cap, turn=224, guard=4):
    """A positive quadratic clipped above the certified physical ground."""
    floor = FLOOR - center
    lower_guard_start = -F(1447, 2) - center
    lower_full = -723 - center
    a = (floor + lower_full) / 2
    k = (turn - a) / (2 * turn**2 - a**2 - half_gap**2)

    def r(x):
        return (x - a) * (1 - k * (x + a))

    assert GROUND_UPPER < center + lower_guard_start
    assert floor < a < lower_guard_start < lower_full < -half_gap < 0
    assert turn > half_gap and a < 0 and a**2 > half_gap**2
    assert k > 0 and 1 - k * (turn + guard + a) > 0
    assert turn**2 > half_gap**2
    # For x>=turn, (x-a+k(a²-d²))/(x²-d²) decreases from 2k to zero.
    # On the other excluded tail use its center root, interval length and
    # minimum of the positive weight x²-d².
    low_ratio = (lower_full - a) * (1 - k * (floor + a)) / (lower_full**2 - half_gap**2)
    tail_ratio = max(low_ratio, k)
    tail = tail_ratio * variance_cap
    assert tail < F(1, 2)
    normalizer = 1 / (4 * k) - a + k * a**2
    lower_slope = r(lower_full) / (lower_full - lower_guard_start)
    maximum_r_derivative = max(
        1 - 2 * k * lower_guard_start, abs(1 - 2 * k * (turn + guard))
    )
    lipschitz = (
        maximum_r_derivative + max(lower_slope, normalizer / guard)
    ) / normalizer
    return {
        "window": [center + lower_guard_start, center + turn + guard],
        "last_edge_of_low_tail": center + lower_full,
        "first_edge_of_high_tail": center + turn,
        "linear_bottom_guard_width": lower_full - lower_guard_start,
        "linear_top_guard_width": guard,
        "low_root_offset_from_half_sum": a,
        "polynomial_coefficient": k,
        "maximum_weighted_tail_ratio": tail_ratio,
        "tail_contribution_upper": tail,
        "polynomial_normalizer": normalizer,
        "ideal_receiver_amplitude_lower": (F(1, 2) - tail) / normalizer,
        "filter_Lipschitz_per_energy_upper": lipschitz,
    }


def clock_plan(radius):
    resolution = F(1, 4096)
    tick = 1 / (4 * radius)
    alias_margin = tick * resolution / 7
    target_modes = 32 / alias_margin
    modes = 1 << ((target_modes.numerator - 1) // target_modes.denominator).bit_length()
    assert modes >= target_modes
    tail = F(1, 2) / (modes * alias_margin) ** 2 + F(1, 2) / (modes * alias_margin)
    reads = 9
    majority_tail = sum(
        comb(reads, j) * tail**j * (1 - tail) ** (reads - j)
        for j in range(reads // 2 + 1, reads + 1)
    )
    width = modes.bit_length() - 1
    allowance = F(1, 2**17)
    return {
        "clock_tick": tick,
        "single_read_register_bits": width,
        "Fourier_register_dimension": modes,
        "median_registers": reads,
        "largest_controlled_native_wait": tick * modes / 2,
        "energy_resolution": resolution,
        "single_read_tail_bound": tail,
        "median_tail_bound": majority_tail,
        "complete_gate_operator_error_allowance": allowance,
        "per_controlled_clock_operator_error_sufficient": allowance
        / (4 * reads * width),
        "ancilla_circuit_operator_error_sufficient": allowance / 2,
    }


def compiled_window(witness, plan, wait):
    lo, hi = witness["window"]
    ideal = witness["ideal_receiver_amplitude_lower"]
    compiled = (
        plan["median_tail_bound"]
        + plan["energy_resolution"] * witness["filter_Lipschitz_per_energy_upper"]
        + plan["complete_gate_operator_error_allowance"]
    )
    flight = wait * (hi - lo) / 2
    assert ideal > flight + compiled
    return {
        **witness,
        "ground_relative_interaction_energy_enclosure": [
            max(0, lo - GROUND_UPPER),
            hi - FLOOR,
        ],
        "compiled_Kraus_operator_error_upper": compiled,
        "additional_native_hold": wait,
        "conditioned_two_spin_transfer_minimum_singular_value_lower": ideal
        - flight
        - compiled,
        "clock_accept_and_full_receiver_at_hold_zero_lower": (ideal - compiled) ** 2,
        "clock_accept_and_full_receiver_at_independent_hold_lower": (
            ideal - flight - compiled
        )
        ** 2,
        "probability_outside_exact_window_and_clock_accept_upper": compiled**2,
    }


def run():
    contract = recoil.read_archive(recoil.CONTRACT)
    ground = recoil.read_archive(GROUND)
    archive = json.loads(ground)
    member = json.loads(contract)["fixed_primary_member"]
    ids = {"member_id": member["member_id"], **member["ids"]}
    charged = recoil.run()
    evidence = recoil.HERE / "results.json"
    assert json.loads(json.dumps(charged, default=str)) == json.loads(
        evidence.read_text()
    )
    assert charged["original_source"] == ids
    upper = archive["selected_member_bound"]["result"]
    (lower,) = [
        r for r in archive["ground_floors"] if r["member_id"] == ids["member_id"]
    ]
    for record in (upper, lower):
        assert all(record[key] == val for key, val in ids.items())
        assert record["full_physical_carrier"]
    assert F(lower["physical_ground_lower"]) >= FLOOR
    assert F(upper["physical_ground_upper"]) <= GROUND_UPPER
    state = member["state"]
    radius = max(abs(F(state[key])) for key in ("spectral_lower", "spectral_upper"))
    plan = clock_plan(radius)
    gamma = recoil.read_archive(recoil.GAMMA)
    for hop in recoil.beta_hops(recoil.clifford(gamma)):
        # Distinct input spin columns can meet on a single-hop Fock word only
        # through the same row or column of a within-beta-pair block.
        for r in (slice(0, 2), slice(2, 4)):
            for s in (slice(0, 2), slice(2, 4)):
                block = hop[r, s]
                assert all(
                    sum(not zo.is_zero for zo in matrix.row(j)) <= 1
                    for matrix in (block, block.T)
                    for j in range(2)
                )
    rows = []
    for row in charged["clocked_rows"]:
        sending = row["absolute_source_clock_diagnostics"]
        receiving = row["absolute_receiver_clock_diagnostics"]
        center = (sending["mean"] + receiving["mean"]) / 2
        half_gap = abs(sending["mean"] - receiving["mean"]) / 2
        variance_cap = sqrt_upper(
            max(sending["unrecombined_sharp_spin_variances"])
            * max(receiving["unrecombined_sharp_spin_variances"])
        )
        compact = low_witness(center, half_gap, variance_cap, F(48))
        onset = low_witness(center, half_gap, variance_cap, F(42), F(1, 4))
        strict = excited_witness(center, half_gap, variance_cap)
        rows.append(
            {
                "path": row["path"],
                "moving_and_holder_hypercharge": row[
                    "native_moving_and_holder_hypercharge"
                ],
                "midpoint_clock_energy": center,
                "half_energy_detuning": half_gap,
                "sqrt_product_of_endpoint_variances_upper": variance_cap,
                "full_spin_matrices": {
                    "zeroth_sender_receiver": "0 times I2",
                    "mixed_first": "-1/2 times I2",
                    "mixed_second_about_midpoint": "0 times I2",
                    "sender_centered_second_diagonal": sending[
                        "unrecombined_sharp_spin_variances"
                    ],
                    "receiver_centered_second_diagonal": receiving[
                        "unrecombined_sharp_spin_variances"
                    ],
                },
                "exact_lower_spectral_window": compiled_window(
                    compact, plan, F(1, 2**18)
                ),
                "existence_of_common_energy_below": onset["first_edge_of_high_tail"],
                "existence_bound_uses_tighter_discontinuous_window_only": True,
                "exact_excited_only_window": compiled_window(strict, plan, F(1, 2**21)),
            }
        )
    return {
        "source_revision": recoil.SOURCE_REVISION,
        "source_identity": ids,
        "inputs_sha256": {
            recoil.CONTRACT: sha256(contract).hexdigest(),
            recoil.GAMMA: sha256(gamma).hexdigest(),
            GROUND: sha256(ground).hexdigest(),
            FOURIER: sha256(recoil.read_archive(FOURIER)).hexdigest(),
            "endpoint_recoil/results.json": sha256(evidence.read_bytes()).hexdigest(),
        },
        "native_clock_energy_zero_changed": False,
        "full_interacting_ground_energy_bracket": [FLOOR, GROUND_UPPER],
        "one_stock_original_admission_lower": charged[
            "starting_stock_original_admission_lower"
        ],
        "conditional_clock_and_receiver_bounds": rows,
        "clock_compilation_plan": plan,
        "exact_opposite_holder_spin_charge_conjugation": energy_bridge.CONJUGATION,
        "full_quantum_clock_circuit_executed": False,
        "isolated_line_measured": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = HERE / "results.json"
    actual = json.dumps(run(), default=str, sort_keys=True, indent=2) + "\n"
    if args.write:
        receipt.write_text(actual)
    else:
        assert receipt.read_text() == actual
        print("Verified full Wilson low-energy and excited-only receiver windows.")


if __name__ == "__main__":
    main()
