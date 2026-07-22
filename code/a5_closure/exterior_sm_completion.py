#!/usr/bin/env python3
"""Exact exterior Standard-Model completion certificate.

This module checks the finite representation, hypercharge, anomaly, Yukawa,
weak-multiplicity, and face-phase arithmetic used by the conditional
screen-to-Standard-Model theorem. It does not construct the PSTRC carrier,
identify coefficient directions with physical currents, select the
trace-balanced global form, attach the face module to physical families, or
exclude extra light sectors.
"""
from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def q(value: int, denominator: int = 1) -> Fraction:
    return Fraction(value, denominator)


# V = C + W with hypercharges fixed by total trace zero:
# 3*(-1/3) + 2*(1/2) = 0.
carrier = {
    "C": {"dimension": 3, "su3": "3", "su2": "1", "hypercharge": q(-1, 3)},
    "W": {"dimension": 2, "su3": "1", "su2": "2", "hypercharge": q(1, 2)},
}
assert 3 * carrier["C"]["hypercharge"] + 2 * carrier["W"]["hypercharge"] == 0


# The non-vacuum even exterior package M1 = Lambda^2 V + Lambda^4 V.
# Each row is (SU3 representation, SU2 representation, Y, complex dimension).
fields = {
    "Q": {"origin": "C tensor W", "su3": "3", "su2": "2", "Y": q(1, 6), "dimension": 6},
    "u_c": {"origin": "Lambda^2 C", "su3": "3bar", "su2": "1", "Y": q(-2, 3), "dimension": 3},
    "e_c": {"origin": "Lambda^2 W", "su3": "1", "su2": "1", "Y": q(1), "dimension": 1},
    "d_c": {
        "origin": "Lambda^2 C tensor Lambda^2 W",
        "su3": "3bar",
        "su2": "1",
        "Y": q(1, 3),
        "dimension": 3,
    },
    "L": {"origin": "Lambda^3 C tensor W", "su3": "1", "su2": "2", "Y": q(-1, 2), "dimension": 2},
}
assert sum(row["dimension"] for row in fields.values()) == 15
assert sum(
    row["dimension"]
    for row in fields.values()
    if row["origin"] in {"C tensor W", "Lambda^2 C", "Lambda^2 W"}
) == 10
assert fields["d_c"]["dimension"] + fields["L"]["dimension"] == 5

field_signatures = {
    (row["su3"], row["su2"], row["Y"])
    for row in fields.values()
}
conjugate_su3 = {"1": "1", "3": "3bar", "3bar": "3"}
dual_signatures = {
    (conjugate_su3[row["su3"]], row["su2"], -row["Y"])
    for row in fields.values()
}
assert field_signatures.isdisjoint(dual_signatures)


# Gauge-invariant one-Higgs Yukawa channels. H = W has Y=+1/2.
scalars = {
    "H": {"su3": "1", "su2": "2", "Y": carrier["W"]["hypercharge"]},
    "Hdag": {"su3": "1", "su2": "2", "Y": -carrier["W"]["hypercharge"]},
}
yukawa_factors = {
    "Q_H_u_c": (fields["Q"], scalars["H"], fields["u_c"]),
    "Q_Hdag_d_c": (fields["Q"], scalars["Hdag"], fields["d_c"]),
    "L_Hdag_e_c": (fields["L"], scalars["Hdag"], fields["e_c"]),
}


def singlet_multiplicity(group: str, factors: tuple[dict, ...]) -> int:
    reps = [row[group] for row in factors if row[group] != "1"]
    if not reps:
        return 1
    if group == "su3" and sorted(reps) == ["3", "3bar"]:
        return 1
    if group == "su2" and reps == ["2", "2"]:
        return 1
    return 0


yukawa_charge_sums = {
    name: sum((row["Y"] for row in factors), q(0))
    for name, factors in yukawa_factors.items()
}
yukawa_invariant_multiplicities = {
    name: {
        "su3": singlet_multiplicity("su3", factors),
        "su2": singlet_multiplicity("su2", factors),
    }
    for name, factors in yukawa_factors.items()
}
assert set(yukawa_charge_sums.values()) == {q(0)}
assert all(row == {"su3": 1, "su2": 1} for row in yukawa_invariant_multiplicities.values())


# Exact one-generation anomaly arithmetic in left-handed Weyl convention.
anomalies = {
    "SU3_cubed": q(2) - q(1) - q(1),
    "SU3_squared_U1": q(2) * q(1, 2) * fields["Q"]["Y"]
    + q(1, 2) * fields["u_c"]["Y"]
    + q(1, 2) * fields["d_c"]["Y"],
    "SU2_squared_U1": q(3) * q(1, 2) * fields["Q"]["Y"]
    + q(1, 2) * fields["L"]["Y"],
    "gravity_squared_U1": sum(row["dimension"] * row["Y"] for row in fields.values()),
    "U1_cubed": sum(row["dimension"] * row["Y"] ** 3 for row in fields.values()),
}
assert set(anomalies.values()) == {q(0)}


# SU(2)-doublet multiplicity: three color copies of Q plus one lepton copy.
weak_doublet_multiplicity = 3 + 1
assert weak_doublet_multiplicity == 4
assert weak_doublet_multiplicity % 2 == 0


# Face-stabilizer induction. For A5 irreps (1,3,3',4,5), the character on
# the order-three class is (1,0,0,1,-1). Frobenius reciprocity gives the
# multiplicity of either non-trivial C3 character as (dim-chi(3A))/3.
a5_face_phase = {
    "1": {"dimension": 1, "chi_3A": 1},
    "3": {"dimension": 3, "chi_3A": 0},
    "3prime": {"dimension": 3, "chi_3A": 0},
    "4": {"dimension": 4, "chi_3A": 1},
    "5": {"dimension": 5, "chi_3A": -1},
}
for row in a5_face_phase.values():
    assert (row["dimension"] - row["chi_3A"]) % 3 == 0
    row["omega_multiplicity"] = (row["dimension"] - row["chi_3A"]) // 3
induced_decomposition = {
    name: row["omega_multiplicity"]
    for name, row in a5_face_phase.items()
    if row["omega_multiplicity"]
}
assert induced_decomposition == {"3": 1, "3prime": 1, "4": 1, "5": 2}
assert sum(
    a5_face_phase[name]["dimension"] * mult
    for name, mult in induced_decomposition.items()
) == 20


# Abstract control for the proposed 24-slot deck reduction. This proves only
# what follows from a declared free action, not that the OPH register carries it.
abstract_deck_orbits = [list(range(start, start + 6)) for start in range(0, 24, 6)]
assert len(abstract_deck_orbits) == 4
assert sorted(x for orbit in abstract_deck_orbits for x in orbit) == list(range(24))


def frac(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


payload = {
    "schema": "OPH conditional exterior Standard-Model completion certificate v1",
    "carrier": {
        name: {**row, "hypercharge": frac(row["hypercharge"])}
        for name, row in carrier.items()
    },
    "matter_package": {
        "definition": "M1 = Lambda^2(C+W) + Lambda^4(C+W)",
        "complex_dimension": 15,
        "fields": {
            name: {**row, "Y": frac(row["Y"])}
            for name, row in fields.items()
        },
        "chirality": "M1* = Lambda^3 V + Lambda^1 V has no irreducible summand in common with M1",
        "global_descent": "automatic once V is a genuine representation of S(U(3)xU(2))",
    },
    "yukawa_channels": {
        name: {
            "hypercharge_sum": frac(value),
            "su3_singlet_multiplicity": yukawa_invariant_multiplicities[name]["su3"],
            "su2_singlet_multiplicity": yukawa_invariant_multiplicities[name]["su2"],
            "invariant_multiplicity_before_families": 1,
        }
        for name, value in yukawa_charge_sums.items()
    },
    "anomalies": {
        "scope": "listed four-dimensional perturbative anomalies",
        "coefficients": {name: frac(value) for name, value in anomalies.items()},
        "su2_witten_parity": (
            "PASS: SU(2) Witten parity is even, with four left-handed weak "
            "doublets modulo color multiplicity"
        ),
    },
    "weak_load": {
        "doublet_multiplicity_per_generation": weak_doublet_multiplicity,
        "decomposition": "(C tensor W) + (Lambda^3 C tensor W) = (C + complex line) tensor W",
        "beta_EW_if_defined_as_per_generation_doublet_multiplicity": 4,
        "three_generation_doublets": 12,
        "reversibly_oriented_three_generation_count": 24,
        "physical_gate": "PORT-LOAD-TRACE and family normalization",
    },
    "face_phase": {
        "induced_module": "Ind_C3^A5(omega)",
        "decomposition": induced_decomposition,
        "minimal_irreps": ["3", "3prime"],
        "physical_gate": "A5-FAMILY-ATTACHMENT",
        "yukawa_boundary": (
            "An exact unbroken A5 action on the family multiplicity space restricts "
            "Yukawa tensors to A5-invariant pairings. General family matrices require "
            "a source-derived breaking, hiding, or forgetting mechanism."
        ),
    },
    "abstract_deck_control": {
        "premise": "a declared free Z6 action on 24 slots",
        "orbit_count": len(abstract_deck_orbits),
        "invariant_function_dimension": len(abstract_deck_orbits),
        "physical_gates": ["physical free deck action", "PORT-WEAK-INTERTWINER"],
    },
    "conditional_closure": {
        "proved_mathematics": [
            "one-generation exterior decomposition and exact hypercharges",
            "chirality at representation level",
            "three one-Higgs Yukawa invariant lines",
            "listed perturbative anomaly cancellation and the SU(2) Witten parity condition",
            "four weak-doublet copies per generation",
            "minimal three-dimensional face-phase extension",
            "four invariant orbits for any declared free Z6 action on 24 slots",
        ],
        "remaining_producers": [
            "certified echosahedral carrier lineage (UD12 and RP-A5 are closed there, but not for arbitrary OPH carriers)",
            "PORT-CURRENT-INNER physical source binding remains open; port_current_inner_certificate.py verifies the conditional algebraic construction for a declared charged-double-triplet representation and four signed coefficients",
            "PORT-SPIN-LIFT physical source binding remains open; super_tannakian_matter_lift_certificate.py verifies the conditional algebraic matter lift for the declared matter-lift contracts: the non-split double cover, the projector-realized fifteen-state module, kernel emission, and MAR nonemptiness",
            "BLOCK-DETERMINANT-BALANCE",
            "AXIS-CENTER-DESCENT and the physical line spectrum (consuming the emitted action-kernel data: generator, deck relation, residual order six)",
            "physical selection of the non-vacuum even exterior package beyond the declared selection contract",
            "physical selection of H=W as the scalar doublet beyond the declared one-scalar contract",
            "A5-FAMILY-ATTACHMENT",
            "A5 family-symmetry breaking, hiding, or forgetting before general Yukawa matrices",
            "PORT-WEAK-INTERTWINER and PORT-LOAD-TRACE",
            "continuum spin/QFT realization",
        ],
        "no_extra_sector_boundary": (
            "M1 is the non-vacuum even exterior package, not the full even Clifford module. "
            "The full even exterior module also contains Lambda^0 V, and M1 alone is not "
            "Clifford-stable. Removing that singlet and "
            "excluding other anomaly-free light sectors requires MGFC/MAR or an "
            "observer-visible discriminator."
        ),
        "physical_status": (
            "conditional structural theorem; no mass, coupling, pole, or "
            "unconditional screen-forcing claim"
        ),
    },
}


if __name__ == "__main__":
    out = HERE / "exterior_sm_completion.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
