"""Resolve the supplied neutral-current coefficients using the Wilson endpoint."""

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
from subprocess import check_output

import sympy as sp

from ..charged_endpoints.cube import CELLS, LINKS, MASKS, Q6
from ..charged_endpoints.measure import (
    BOUNDS,
    CHARGE_READS,
    CONTRACT,
    GAMMA,
    ROOT,
    SOURCE_REVISION,
    clifford,
    matrices,
    read_archive,
)
from ..neutral_direction import verify as original_direction_check
from .wilson_probe import DOUBLET, verify_link_and_frame

HERE = Path(__file__).resolve().parent
OCCUPATION = "model/particle_primary_occupation_bound/results.json"
DICTIONARY = "code/defect_particle_classes/results.json"
DICTIONARY_REVISION = "3078a56fd1ed0a9c818a5966b2d1efb7e36e72bb"
NEUTRAL_LEAN = "Lean/Screen/NeutralCurrentDictionary.lean"
LOCAL_LEAN = HERE.parent / "EndpointNeutralDirection.lean"
CALIBRATIONS = {
    "hypercharge": (0, 1),
    "combined": (1, 1),
    "reversed_weak": (-1, 1),
    "weak_only": (1, 0),
    "sham": (0, 0),
}
RADIX = 5
TRANSFER_DIGITS = [0, 6, -3, 3, -6]
NATIVE_BLOCKS = (
    (3, -4, Fraction(4, 5)),
    (2, -3, Fraction(9, 20)),
    (6, 1, Fraction(21, 20)),
    (3, 2, Fraction(7, 10)),
    (1, 6, Fraction(3, 10)),
)


def sqrt_up(x):
    unit = 2**72
    radicand = x * unit**2
    ceil = isqrt(radicand.numerator // radicand.denominator)
    ceil += int(ceil**2 < radicand)
    result = Fraction(ceil, unit)
    assert result**2 >= x and (result - Fraction(1, unit)) ** 2 <= x
    return result


def allocation_loss(m, d, inverse_margin):
    q = 2 * (d + 1)
    w = inverse_margin
    active_defects = 4

    def lower(defects):
        return defects * m - w * q

    def shift(defects):
        bottom = lower(defects)
        return (d * w - bottom + sqrt_up((bottom + d * w) ** 2 + 4 * w**2)) / 2

    c = shift(active_defects)
    if lower(2) + c <= 0:
        active_defects = 2
        c = shift(active_defects)
    c = max(c, w * q - m, Fraction())
    assert c >= w * d and lower(2) + c > 0
    assert (c - w * d) * (lower(active_defects) + c) >= w**2
    return c


def control_admission(occupation, contract, success):
    """Bound the original local stock acceptance by the source energy budget."""
    source = contract["source"]
    assert source["side"] == 2 and Fraction(source["spacing"]) == 1
    assert Fraction(source["electric_coefficient"]) == 105
    assert source["interior_plaquettes"] == 6 and source["interior_hopping_bonds"] == 12
    assert (
        tuple(
            (
                row["dimension"],
                row["integer_Q6_charge"],
                Fraction(row["trace_metric_Casimir"]),
            )
            for row in occupation["local_bound"]["native_blocks"]
        )
        == NATIVE_BLOCKS
    )
    # Retain beta defect energy only on the colorless weak doublet.
    r = Fraction(13, 20)
    k = 48 * r
    s = k / 189
    electric = 105 * (1 - s)
    blocks = [(d, q, c, 1 - r if q == -3 else Fraction(1)) for d, q, c in NATIVE_BLOCKS]
    shift = sum(4 * d * (sqrt_up(m**2 + Fraction(1, 4)) - m) for d, _, _, m in blocks)
    worst_margin = electric * Fraction(3, 20) - shift
    assert worst_margin > 0
    losses = []
    for d, q, casimir, m in blocks:
        margin = electric * casimir / 2 - shift
        loss = allocation_loss(m, d, 1 / margin)
        losses.append({"dimension": d, "native_transfer": q, "loss": str(loss)})
    total_loss = 12 * sum(Fraction(x["loss"]) for x in losses) + Fraction(6, 7)
    source_upper = Fraction("-723.7678199426")
    assert source_upper >= Fraction(occupation["physical_ground_upper"])
    excess = Fraction(contract["state"]["target_ground_excess"])
    assert excess == Fraction(1, 16)
    assert (
        Fraction(occupation["selected_mean_source_energy_upper"])
        <= source_upper + excess
    )
    upper = source_upper + excess + 720 + total_loss
    # Sum the two selected empty/full doublet modes over all 48 addresses.
    # Each corresponding individual beta defect is charged exactly 3 times;
    # each nonconstant link counts 4 times. C >= (3/10)(I-P0) on every bond.
    directions = [(x, y) for a, b in LINKS for x, y in ((a, b), (b, a))]
    for sign in (1, -1):
        charged = Counter()
        addressed = Counter()
        for x, y in directions:
            for spin in (0, 1):
                charged.update([(x, sign, spin), (y, -sign, spin)])
                addressed[tuple(sorted((x, y)))] += 1
        assert set(charged.values()) == {3} and len(charged) == 32
        assert set(addressed.values()) == {4} and len(addressed) == 12
    assert k * Fraction(1, 16) == 3 * r
    assert k * Fraction(5, 18) == 105 * s / 2
    lower = 1 - upper / k
    assert lower > Fraction("0.78320864")
    assert Fraction(success) > Fraction("0.79446654")
    return {
        "original_stock_probability_lower_exact": str(lower),
        "simple_probability_lower": "0.78320864",
        "same_bound_each_reference_orientation_and_each_beta_sign": True,
        "full_source_polynomial_preparation_used": True,
        "uniform_spatial_and_spin_address_average": 48,
        "original_gate_refusal_probability_upper": str(1 - lower),
        "retained_weak_doublet_fraction": str(r),
        "retained_electric_fraction": str(s),
        "original_selected_energy_upper": str(source_upper + excess),
        "full_fock_schur_shift_upper": str(total_loss),
        "weak_doublet_beta_defect_budget_coefficient": str(3 * r),
        "electric_Casimir_budget_coefficient": str(105 * s / 2),
        "certified_nonconstant_link_margin_lower": str(worst_margin),
        "per_native_block_vacuum_resolvent_enclosures": losses,
    }


def light_ports(phase):
    # One balanced branch carries the signed transfer and the second
    # carries the transfer with a current pulse conjugation.
    path = sp.Matrix([1, phase]) / sp.sqrt(2)
    x_meter = sp.Matrix([[1, 1], [1, -1]]) / sp.sqrt(2)
    y_meter = sp.Matrix([[1, -sp.I], [1, sp.I]]) / sp.sqrt(2)
    result = {}
    for name, meter in (("cosine", x_meter), ("sine", y_meter)):
        assert meter.H * meter == sp.eye(2)
        ports = [sp.simplify(a * sp.conjugate(a)) for a in meter * path]
        assert sum(ports) == 1
        result[name] = ports
    return result


def probability_law(response, admission_floor):
    phase = sp.expand_complex(sp.exp(sp.I * sp.pi * response / 2))
    assert sp.simplify(phase * phase.conjugate()) == 1
    ports = light_ports(phase)
    readout = sp.simplify(
        2 * ports["cosine"][0] - 1 + sp.I * (2 * ports["sine"][0] - 1)
    )
    assert sp.simplify(readout - phase) == 0
    floor = sp.Rational(admission_floor)
    return {
        "signed_response_reconstructed_from_phase": str(response),
        "unit_phase": str(phase),
        "conditioned_plus_minus": {
            axis: list(map(str, row)) for axis, row in ports.items()
        },
        "original_plus_minus_weight_algebraic_lower": {
            axis: list(map(str, (sp.simplify(p * floor) for p in row)))
            for axis, row in ports.items()
        },
    }


def transfer_geometry(archive):
    shapes = archive["source_probe_charge_records"]
    (failure,) = [row for row in shapes if row.get("failure")]
    assert failure["charge_shift"] == [0] * 8
    successes = [r for r in shapes if not r.get("failure")]
    seen = Counter()
    joint = []
    for record in successes:
        x, y = tuple(record["target"]), tuple(record["origin"])
        native = record["charge_shift"]
        assert tuple(sorted((x, y))) in LINKS
        assert native == [6 * (int(z == x) - int(z == y)) for z in CELLS]
        seen.update([(x, y, record["matched_beta_component"])])
        joint.append(
            {
                "target": list(x),
                "source": list(y),
                "matched_beta_component": record["matched_beta_component"],
                "spinor_native_transfer": native,
                "doublet_particle_native_transfer": [-q // 2 for q in native],
                "doublet_hole_native_transfer": [q // 2 for q in native],
            }
        )
    assert len(seen) == 48 and len(successes) == 48
    assert set(seen.values()) == {1}
    for q in TRANSFER_DIGITS:
        assert TRANSFER_DIGITS[q % RADIX] == q
    assert archive["charge_menu"]["difference_meter_radix"] == RADIX
    for row in [*joint, {"spinor_native_transfer": [0] * 8}]:
        assert all(
            TRANSFER_DIGITS[q % RADIX] == q
            for name, shifts in row.items()
            if name.endswith("native_transfer")
            for q in shifts
        )
    return joint


def calibrated_records(table, sign, native_floor, control_floor):
    def response(mask, a, b):
        entry = table[mask]
        return sign * (a * sp.Rational(entry["T3"]) + b * sp.Rational(entry["Y"]))

    output = []
    for label, (a, b) in CALIBRATIONS.items():
        s = response(24, a, b)
        u = response(DOUBLET[0], a, b)
        lower = response(DOUBLET[1], a, b)
        assert (sign * (s + 2 * u), sign * s) == (a, b)
        assert lower == -s - u
        output.append(
            {
                "meter_setting": label,
                "source_endpoint_at_named_receiver": probability_law(s, native_floor),
                "other_address_or_failure_phase": probability_law(sp.S(0), Fraction(0)),
                "independent_upper_control": probability_law(u, control_floor),
                "held_independent_lower_control": probability_law(lower, control_floor),
            }
        )
    return output


def run():
    files = (
        CONTRACT,
        BOUNDS,
        CHARGE_READS,
        GAMMA,
        OCCUPATION,
        NEUTRAL_LEAN,
    )
    frozen = {path: read_archive(path) for path in files}
    dictionary_file = check_output(
        ["git", "show", f"{DICTIONARY_REVISION}:{DICTIONARY}"], cwd=ROOT
    )
    contract = json.loads(frozen[CONTRACT])["fixed_primary_member"]
    band, charged, occupation = (
        json.loads(frozen[path])["selected_member_bound"]["result"]
        for path in (BOUNDS, CHARGE_READS, OCCUPATION)
    )
    identity = {"member_id": contract["member_id"], **contract["ids"]}
    for value in (band, charged, occupation):
        assert all(value[k] == expected for k, expected in identity.items())
    dictionary = {r["mask"]: r for r in json.loads(dictionary_file)["patterns"]}
    assert MASKS[Q6] == 24
    assert (dictionary[24]["n2"], dictionary[24]["n3"]) == (2, 0)
    assert all((dictionary[m]["n2"], dictionary[m]["n3"]) == (1, 3) for m in DOUBLET)
    weights = {m: dictionary[m]["six_Y"] for m in MASKS}
    assert sorted(Counter(weights.values()).items()) == sorted(
        (q, d) for d, q, _ in NATIVE_BLOCKS
    )
    # Four Dirac copies; one colorless and three colored weak doublets
    # per copy. Keep every local weak mode in the common current pulse.
    weak_doublet_copies = 4 * sum(dictionary[m]["n2"] == 1 for m in MASKS) // 2
    assert weak_doublet_copies == 16
    spinor = matrices(clifford(frozen[GAMMA]))
    legacy = original_direction_check()
    assert legacy["endpoint_probes_bound_to_original_transfer_alphabet"] == 48
    assert legacy["held_charged_doublet_response"] == "-s-u"
    geometry = transfer_geometry(charged)
    source_floor = Fraction(band["original_success_weight_lower"])
    energy_floor = Fraction(band["reported_transfer_original_weight_lower"])
    assert source_floor > Fraction("0.79446654") and energy_floor > Fraction(1, 4)
    assert energy_floor == Fraction(
        charged["reported_joint_opposite_charge_transfer_original_weight_lower"]
    )
    floor = control_admission(occupation, contract, source_floor)
    control_floor = Fraction(floor["original_stock_probability_lower_exact"])
    # These probabilities use complete original starts. The same local stock
    # event precedes both current settings, so the common weight also multiplies
    # their total variation on the recorded sine port.
    exact_control_floor = sp.Rational(control_floor)
    upper_difference = exact_control_floor * sp.sqrt(2) / 4
    lower_difference = exact_control_floor * (2 - sp.sqrt(2)) / 4
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "immutable_source_revision": SOURCE_REVISION,
        "immutable_source_hashes": {
            path: sha256(data).hexdigest() for path, data in frozen.items()
        },
        "exterior_dictionary_revision": DICTIONARY_REVISION,
        "exterior_dictionary_sha256": sha256(dictionary_file).hexdigest(),
        "endpoint_theorem_sha256": sha256(LOCAL_LEAN.read_bytes()).hexdigest(),
        "source_identity": identity,
        "same_spinor_modes": spinor,
        "actual_directed_bonds_spin_and_native_transfer_decoder": geometry,
        "expanded_radix_five_decoder": TRANSFER_DIGITS,
        "bigger_meter_or_changed_source_required_for_decoder": False,
        "reference_setting_phase_exposure_in_pi_units": "1/2",
        "address_uniformly_chosen_before_string_on_original_bank": 8,
        "named_endpoint_each_sign_per_original_probability_lower": str(
            source_floor / 8
        ),
        "same_energy_interval": band["reported_transfer_interval"],
        "named_endpoint_each_sign_energy_joint_original_weight_lower": str(
            energy_floor / 8
        ),
        "native_string_total_failure_probability_upper": str(1 - source_floor),
        "success_off_named_receiver_original_weight_lower": str(3 * source_floor / 4),
        "sum_phase_ports_is_same_original_charge_energy_and_idle_reference_channel": True,
        "same_full_polynomial_source_upper_and_lower_gate": floor,
        "pure_Y_vs_combined_total_variation_lower_per_original_start_by_control": {
            "upper": str(upper_difference),
            "held_lower": str(lower_difference),
        },
        "gauge_and_weak_checks": verify_link_and_frame(weights, weak_doublet_copies),
        "source_using_original_clock_denominations": {
            "endpoint_onsite_and_electric_center": "87/4",
            "independent_control_onsite_plus_electric_center": "237/8",
            "control_energy_full_histogram_executed": False,
        },
        "indistinguishability_for_all_original_singlet_input_records": "(a,b) ↦ (a′,b)",
        "signed_banks": [
            {
                "receiver_particle_hole_sign": sign,
                "spinor_native_receiver_increment": sign * 6,
                "doublet_native_receiver_increment": -sign * 3,
                "same_current_settings": calibrated_records(
                    dictionary, sign, source_floor / 8, control_floor
                ),
            }
            for sign in (1, -1)
        ],
        "full_canonical_source_histogram_executed": False,
        "higgs_dynamics_or_native_neutral_mode_identified": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as out:
            out.write(
                json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
            )
    else:
        assert result == json.loads((HERE / "records.json").read_text())
        print(
            "Replayed the signed Wilson currents, full-source stock bound and frame-erasure controls."
        )


if __name__ == "__main__":
    main()
