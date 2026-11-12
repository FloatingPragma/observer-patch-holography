"""Exact charge lattice of the declared five-level exterior code.

Fresh implementation of the equations cited in source_refs. No result receipt
or measured particle table enters the calculation. SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import product
import json
from math import comb


SOURCE_REFS = [
    "model/defect_particle_classes/model.py",
    "model/defect_particle_classes/PROOF.md",
    "model/defect_integer_lift/PROOF.md",
    "model/native_quantum_gauss/PROOF.md",
]


def enumerate_code(hypercharge=None, electric_scale=F(1)):
    """Occupations in Lambda(C^3+ C^2), with Q=T3+Y in declared units."""
    y_levels = tuple(hypercharge or (F(-1, 3),)*3+(F(1, 2),)*2)
    if len(y_levels) != 5 or any(type(x) not in (int, F) for x in y_levels):
        raise ValueError("five exact rational hypercharges required")
    if type(electric_scale) not in (int, F):
        raise ValueError("exact electric normalization required")
    rows = []
    for mask in range(32):
        bits = tuple((mask >> i) & 1 for i in range(5))
        n3, n2 = sum(bits[:3]), sum(bits[3:])
        y = sum((b*w for b, w in zip(bits, y_levels)), F())
        t3 = F(bits[3]-bits[4], 2)
        rows.append({"mask": mask, "occupation": bits, "n3": n3, "n2": n2,
                     "Y": y, "T3": t3, "Q": electric_scale*(t3+y),
                     "k": (3*n2-2*n3) % 6,
                     "colour_singlet": n3 in (0, 3)})
    return rows


def verify_charge_lattice(rows):
    """Check the stated code and its exact electric and centre consequences."""
    if len(rows) != 32 or {r["mask"] for r in rows} != set(range(32)):
        raise ValueError("the complete occupation code is required")
    if any(r["Q"] != F(r["occupation"][3])-F(r["n3"], 3) for r in rows):
        raise ValueError("declared electric normalization or weight law changed")
    if any(r["Y"] != F(3*r["n2"]-2*r["n3"], 6) for r in rows):
        raise ValueError("hypercharge weight law changed")
    if any((F(r["n3"], 3)+F(r["n2"], 2)+r["Y"]).denominator != 1 for r in rows):
        raise ValueError("the Z6 kernel is not trivial")
    if any((3*r["Q"]).denominator != 1 for r in rows):
        raise ValueError("electric charge is outside one-third integer units")
    if any(r["colour_singlet"] and r["Q"].denominator != 1 for r in rows):
        raise ValueError("a colour-singlet code weight has fractional charge")
    if any(abs(r["Q"]) == F(1, 2) for r in rows):
        raise ValueError("the code contains a half-charge weight")
    combinatorial = Counter()
    for n3, n2 in product(range(4), range(3)):
        combinatorial[(3*n2-2*n3) % 6] += comb(3, n3)*comb(2, n2)
    if Counter(r["k"] for r in rows) != combinatorial:
        raise ValueError("centre multiplicities disagree with the exterior dimensions")
    return True


def _rejects(rows):
    try:
        verify_charge_lattice(rows)
    except ValueError:
        return True
    return False


def run():
    rows = enumerate_code()
    verify_charge_lattice(rows)
    half_classes = {str(k): [j for j in range(6) if (2*j-k) % 6 == 0]
                    for k in range(6)}
    singlets = [r for r in rows if r["colour_singlet"]]
    # This exhausts pairs, including differently charged daughters. It does
    # not assume charge conservation is sufficient for a physical transition.
    electron = next(r for r in rows if r["Q"] == -1 and r["n2"] == 0)
    algebraic_splits = [(a, b) for a, b in product(rows, repeat=2)
                        if a["Q"]+b["Q"] == electron["Q"]
                        and (a["k"]+b["k"]-electron["k"]) % 6 == 0]
    fractional_split = next((a, b) for a, b in algebraic_splits
                            if a["Q"] == F(-2, 3) and b["Q"] == F(-1, 3))
    checks = {
        "complete_32_occupation_patterns": len(rows) == 32,
        "exact_Z6_kernel_on_all_patterns": all(
            (F(r["n3"], 3)+F(r["n2"], 2)+r["Y"]).denominator == 1 for r in rows),
        "all_charges_in_one_third_integer_lattice": all((3*r["Q"]).denominator == 1 for r in rows),
        "all_eight_colour_singlet_weights_integral": len(singlets) == 8 and all(r["Q"].denominator == 1 for r in singlets),
        "half_unit_charge_absent_from_entire_code": all(abs(r["Q"]) != F(1, 2) for r in rows),
        "equal_class_halves_of_odd_k_absent": all(not half_classes[str(k)] for k in (1, 3, 5)),
        "electron_k0_has_equal_class_solutions": half_classes["0"] == [0, 3],
        "unequal_fractional_split_is_algebraically_allowed": bool(fractional_split),
    }
    mutations = {
        "half_electric_normalization_rejected": _rejects(enumerate_code(electric_scale=F(1, 2))),
        "changed_weak_hypercharge_rejected": _rejects(enumerate_code((F(-1, 3),)*3+(F(1, 4),)*2)),
        "missing_code_weight_rejected": _rejects(rows[:-1]),
    }
    if not all(checks.values()) or not all(mutations.values()):
        raise ArithmeticError("charge calculation or a negative control failed")
    serial = lambda r: {k: str(v) if isinstance(v, F) else v for k, v in r.items()}
    return {
        "schema": "oph.openwave.charge.v1", "computational_pass": True,
        "criterion_rows": [1, 12], "full_criterion_met": False,
        "checks": checks, "mutations": mutations, "source_refs": SOURCE_REFS,
        "patterns": [serial(r) for r in rows],
        "centre_multiplicities": [sum(r["k"] == k for r in rows) for k in range(6)],
        "electric_weights": {str(q): n for q, n in sorted(Counter(r["Q"] for r in rows).items())},
        "colour_singlet_masks": [r["mask"] for r in singlets],
        "equal_class_halves": half_classes,
        "electron_charge_and_class": {"Q": "-1", "k": electron["k"]},
        "algebraically_allowed_coloured_daughters": [serial(r) for r in fractional_split],
        "parameter_ancestry": {
            "declared": ["five internal levels decomposed as three colour plus two weak levels",
                         "exterior occupation code", "hypercharges (-1/3,-1/3,-1/3,1/2,1/2)",
                         "T3=(b_s1-b_s2)/2", "electric generator Q=T3+Y and its unit normalization"],
            "empirical_calibrations_used": [], "observed_particle_names_used_as_inputs": False,
        },
        "scope": {
            "licenses": "Exact charge-lattice and centre statements on the declared 32-weight code, plus the tensor-product implication: total colour triality zero implies integral electric weight.",
            "does_not_license": "No integrated winding charge, energetic protection against fractionation, isolation/confinement theorem, physical quark identification or electric coupling value is computed.",
            "electron_modular_no_halving_argument_refuted": True,
            "charge_conservation_does_not_establish_a_transition": True,
            "official_OpenWave_cell_assigned": False,
        },
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
