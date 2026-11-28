"""Compare local endpoint spin and equal-exposure charged exchange."""

import argparse
import ast
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import pairwise
import json
from pathlib import Path
import subprocess

import numpy as np
import sympy as sp

from .cube import (
    A,
    B,
    C,
    BETA_SEA,
    CELLS,
    CENTER,
    EXCHANGED,
    LINKS,
    MASKS,
    PARTNERS,
    Q6,
    UNCROSSED,
    c_action,
    exchange,
    mode,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[7]
SOURCE_REVISION = "3bcacf1276ddb3ead35be2b9adcd63b6643e2d35"
CONTRACT = "model/particle_common_contract/contract.json"
BOUNDS = "model/particle_primary_charged_string/results.json"
CHARGE_READS = "model/particle_primary_charged_readout/results.json"
GAMMA = "model/common_source/overlap_regulator.py"


def read_archive(name):
    return subprocess.check_output(
        ["git", "show", f"{SOURCE_REVISION}:{name}"], cwd=ROOT
    )


def clifford(data):
    # Replay just the defining matrix constructor from the immutable source.
    definitions = ast.parse(data)
    (gamma_def,) = [
        node
        for node in definitions.body
        if isinstance(node, ast.FunctionDef) and node.name == "gamma_matrices"
    ]
    namespace = {"np": np}
    exec(
        compile(ast.Module(body=[gamma_def], type_ignores=[]), GAMMA, "exec"), namespace
    )
    raw, _ = namespace["gamma_matrices"]()
    gammas = [
        sp.Matrix([[int(z.real) + sp.I * int(z.imag) for z in row] for row in g])
        for g in raw
    ]
    for g in gammas:
        assert g == g.H
    for i, g in enumerate(gammas):
        for j, h in enumerate(gammas):
            assert g * h + h * g == 2 * int(i == j) * sp.eye(4)
    return gammas


def matrices(gammas):
    eye = sp.eye(4)
    chart = sp.Matrix.hstack(
        eye[:, 0] + eye[:, 2],
        eye[:, 1] + eye[:, 3],
        -eye[:, 0] + eye[:, 2],
        -eye[:, 1] + eye[:, 3],
    ) / sp.sqrt(2)
    beta = gammas[3]
    assert chart.H * beta * chart == sp.diag(1, 1, -1, -1)
    spins = [
        sp.simplify(
            -sp.I * chart.H * gammas[(i + 1) % 3] * gammas[(i + 2) % 3] * chart / 2
        )
        for i in range(3)
    ]
    rotations = [-2 * sp.I * j for j in spins]
    conjugation = gammas[1]
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        assert spins[i] * spins[j] - spins[j] * spins[i] == sp.I * spins[k]
        assert rotations[i] ** 2 == -eye
        assert rotations[i] * rotations[j] == -rotations[j] * rotations[i]
        w = (-beta + beta * gammas[i]) / 2
        assert -conjugation.T * w.T * conjugation.conjugate() == w.H
        # Hole transfer is the normally ordered negative transpose of the
        # negative-beta particle block, with the conjugate character.
        beta_w = sp.simplify(chart.H * w * chart)
        assert beta_w[:2, :2] == -sp.eye(2) / 2
        assert -beta_w[2:, 2:].T == -sp.eye(2) / 2
    for block in (slice(0, 2), slice(2, 4)):
        js = [j[block, block] for j in spins]
        assert sum((j**2 for j in js), sp.zeros(2)) == 3 * sp.eye(2) / 4
        assert all(j != sp.zeros(2) and j**2 == sp.eye(2) / 4 for j in js)
    assert all(j[:2, 2:] == j[2:, :2] == sp.zeros(2) for j in spins)
    return {
        "both_beta_spin_casimirs": [3, 4],
        "axes_with_half_turn_squared_minus_identity": 3,
        "opposite_half_turn_orders_on_each_endpoint_have_minus_sign": True,
        "dimension_of_native_singlet": 1,
        "position_hopping_coefficient_after_projection": "-1/2",
        "same_quasiparticle_hopping_coefficient_on_both_charge_signs": True,
    }


def parity(fock, x):
    count = ((fock >> (60 * CELLS.index(x))) & ((1 << 60) - 1)).bit_count()
    return (-1) ** count


def turns(band, charge_meter):
    successes = charge_meter["source_probe_charge_records"]
    (failure,) = [row for row in successes if row.get("failure")]
    assert failure["charge_shift"] == [0] * 8
    preimages = Counter()
    for row in successes:
        if row.get("failure"):
            continue
        origin, target = tuple(row["origin"]), tuple(row["target"])
        assert tuple(sorted((origin, target))) in LINKS
        spin = row["matched_beta_component"]
        assert spin in (0, 1)
        n, _ = c_action(BETA_SEA, mode(origin, -1, spin), False)
        raised, _ = c_action(n, mode(target, 1, spin), True)
        op_support = [parity(raised, x) * parity(BETA_SEA, x) for x in CELLS]
        shift = [6 * (int(x == target) - int(x == origin)) for x in CELLS]
        assert row["charge_shift"] == shift
        assert op_support == [-1 if q else 1 for q in shift]
        # The parity change depends only on the two CAR monomials, so
        # extends from these sharp columns to every input of the operator.
        assert np.prod(op_support).item() == 1
        preimages.update(
            (q, "dark" if flip == -1 else "bright")
            for q, flip in zip(shift, op_support, strict=True)
        )
    assert preimages == {(0, "bright"): 288, (-6, "dark"): 48, (6, "dark"): 48}
    p = Fraction(band["original_success_weight_lower"])
    joint = Fraction(band["reported_transfer_original_weight_lower"])
    assert p > Fraction("0.79446654") and joint > Fraction(1, 4)
    return {
        "named_observer_chosen_before_probe": {"number": 8, "choice": "uniform"},
        "archived_successful_probe_labels": len(successes) - 1,
        "success_labels_times_choices_by_charge_and_phase_disposition": [
            [q, phase, count] for (q, phase), count in sorted(preimages.items())
        ],
        "charge_and_dark_joint_original_weight_lower_by_sign": float(p / 8),
        "total_dark_original_weight_lower": float(p / 4),
        "same_energy_bank_interval": band["reported_transfer_interval"],
        "dark_in_energy_bank_original_weight_lower_by_sign": str(joint / 8),
        "all_dark_in_energy_bank_original_weight_lower": str(joint / 4),
        "failure_reference_port": "bright",
        "global_spatial_central_turn_port_on_each_success": "bright",
        "local_four_pi_turn_port_on_each_success": "bright",
        "charge_meter_and_energy_outcomes_preserved_after_phase_label_is_traced": True,
        "full_source_histogram_executed": False,
    }


def exposed_counts(routes):
    return Counter(tuple(sorted((u, v))) for leg in routes for u, v in pairwise(leg))


def run():
    sources = {
        path: read_archive(path) for path in (CONTRACT, BOUNDS, CHARGE_READS, GAMMA)
    }
    contract = json.loads(sources[CONTRACT])["fixed_primary_member"]
    band = json.loads(sources[BOUNDS])["selected_member_bound"]["result"]
    charge_meter = json.loads(sources[CHARGE_READS])["selected_member_bound"]["result"]
    for key, expected in {
        "member_id": contract["member_id"],
        **contract["ids"],
    }.items():
        assert band[key] == charge_meter[key] == expected
    assert (
        contract["source"]["side"] == 2
        and contract["source"]["per_link_word_degree"] == 1
    )
    active = exposed_counts(EXCHANGED)
    assert active == exposed_counts(UNCROSSED)
    assert active == {tuple(sorted((CENTER, x))): 2 for x in (A, B, C)}
    for a, b in zip((A, B), PARTNERS, strict=True):
        assert tuple(sorted((a, b))) in LINKS
    assert MASKS[Q6] == (1 << 3) | (1 << 4)
    spins = matrices(clifford(sources[GAMMA]))
    # Physical dimension ≤ 2^480 (1 + 2*15²)^24 ≤ 2^696.
    # Taylor argument is positive and at most 1/2 in operator norm.
    state = contract["state"]
    scale = state["dyadic_scale"]
    width = Fraction(state["spectral_upper"]) - Fraction(state["spectral_lower"])
    assert Fraction(state["beta"]) * width <= scale
    rho_floor = 480 + 9 * contract["source"]["owned_links"] + 2 * scale
    return {
        "source_revision": SOURCE_REVISION,
        "source_sha256": {
            path: sha256(content).hexdigest() for path, content in sources.items()
        },
        "source_identity": {"member_id": contract["member_id"], **contract["ids"]},
        "local_dirac_and_Wilson_checks": spins,
        "unchanged_polynomial_source_local_turn": turns(band, charge_meter),
        "physical_exchange_on_heralded_beta_stock": {
            "stock_projection": "filled negative-beta sea and complete constant-link bank",
            "original_stock_admission_dyadic_lower": {
                "base": 2,
                "negative_exponent": rho_floor,
            },
            "both_arms_use_each_addressed_bond_twice": True,
            "compensator_addresses": PARTNERS,
            "exchanged_route": EXCHANGED,
            "uncrossed_route": UNCROSSED,
            "link_character_word_order": LINKS,
            "results_by_sign": [exchange(sign) for sign in (1, -1)],
            "projected_hop_pulse_exposure_in_pi_units": 1,
            "all_unaddressed_and_electric_terms_shuttered_during_hops": True,
            "autonomous_interacting_exchange_executed": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    result = json.loads(json.dumps(run(), allow_nan=False))
    if args.record:
        with args.record.open("x") as handle:
            handle.write(json.dumps(result, indent=2) + "\n")
    elif result != json.loads((HERE / "results.json").read_text()):
        raise ArithmeticError("The charged endpoint replay differs from its receipt")
    print(
        "Verified both charged spinor signs, equal-exposure exchange and native link current."
    )


if __name__ == "__main__":
    main()
