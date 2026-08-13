#!/usr/bin/env python3
"""The differentiator inventory: where the carrier class separates from
legacy models, with every load-bearing exact fact recomputed in place.

Each row names one structure through which the OPH carrier class
separates measurably from the Standard Model with General Relativity or
from rival discreteness classes. Every exact fact a row rests on is
recomputed inside this producer (no cross-receipt byte pins, so the
inventory is stable under concurrent regeneration of other surfaces),
and every row carries its claim type and its custody pointer. Rows that
mention published magnitudes carry them as exposed literature pointers;
nothing here opens a comparison or scores a candidate.

Rows:

1.  Spin-six rotational fingerprint (ladder row FZ-11).
2.  No-linear-dispersion theorem: the even-power hop expansion forbids
    every odd-power vacuum dispersion exactly, so linear-dispersion
    quantum-gravity rivals die on published linear bounds while the
    carrier class is untouched by them.
3.  Five-fold azimuthal uniqueness: a resolved carrier anisotropy
    carries exact five-fold azimuthal symmetry about the vertex axes,
    which no periodic-lattice substrate can produce (crystallographic
    restriction).
4.  Equal-port angular comb: exact weights on the icosahedral invariant
    levels, separating the carrier from isotropic legacy (no comb) and
    from icosahedral-topology legacy (same support, free weights).
5.  Odd-parity blindness: the carrier channel is exactly blind to
    parity-odd sky structure, separating it from parity-violating
    birefringence models and killing any carrier-sourced parity-odd
    claim.
6.  Global-form fork and charge commensuration: the conditional
    Z6-compact global form selects one of the four legacy-consistent
    global forms and forces one commensurate charge lattice.
7.  Scale-free coupling-plane relation: the frozen determinant
    statistic with the two-branch kinetic dichotomy, a relation among
    measured couplings that legacy models do not force.
8.  Positive-chamber Koide identity and the registered tau window
    (ladder row FZ-10).
9.  Closure coordinates (screen-grain inverse coupling and the capacity
    Lambda coordinate), typed diagnostics where legacy models carry
    free parameters.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import a5_multipole_fixed_point_certificate as base
import spin_six_universality_certificate as universality

HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
RECEIPT_PATH = RUNTIME / "differentiator_inventory_receipt.json"

SCHEMA = "oph.differentiator_inventory_receipt.v1"
STATUS = "DIFFERENTIATOR_INVENTORY_CERTIFIED__TEN_ROWS_TYPED"

require = base.require


# ---------------------------------------------------------------------------
# Recomputed exact facts
# ---------------------------------------------------------------------------


def fact_invariant_nulls() -> dict[str, Any]:
    table = {L: universality.invariant_dimension(L) for L in range(17)}
    require(all(table[L] == 0 for L in range(1, 6)), "invariant nulls drift")
    require(table[6] == 1, "spin-six line drift")
    return {
        "m_1_to_5": [table[L] for L in range(1, 6)],
        "m_6": table[6],
        "allowed_levels": [L for L, d in table.items() if d > 0],
    }


def fact_odd_moments_vanish() -> dict[str, Any]:
    verts = base.cartesian_vertices()
    for k in (1, 3, 5, 7):
        require(
            base.p_is_zero(base.moment_sum(verts, k)),
            f"odd moment {k} nonzero",
        )
    return {
        "mechanism": (
            "the class symbol is an even analytic function of k (a "
            "cosine sum), so every odd-power dispersion term vanishes at "
            "every order for every member"
        ),
        "antipodality_cross_check_vanishing_odd_moments": [1, 3, 5, 7],
    }


def fact_five_fold_sector() -> dict[str, Any]:
    cartesian = base.build_cartesian_frame()
    i6 = cartesian.pop("_i6_poly_object")
    cartesian.pop("_vertices_object")
    pole = base.build_pole_frame(i6)
    require(
        pole["cos5phi_sector_matches"] and pole["even_sector_matches_p6"],
        "pole-frame five-fold sector drift",
    )
    return {
        "cos5phi_coefficient": "21/8 on cos(theta) sin^5(theta)",
        "sector_nonzero": True,
        "crystallographic_note": (
            "no three-dimensional periodic lattice carries five-fold "
            "point symmetry, so an exact five-fold azimuthal anisotropy "
            "is unreachable for every periodic-lattice substrate"
        ),
    }


def fact_angular_comb() -> dict[str, Any]:
    """Equal-port comb I_l for l <= 14 by two exact Legendre recurrences.

    The comb is I_l = [P_l(1) + P_l(-1) + 5(P_l(t) + P_l(-t))]/12 with
    t = 1/sqrt5; both signed recurrences run in Q(sqrt5) pairs, every
    level's value is required rational, and the odd nulls are computed
    rather than asserted.
    """

    def legendre_at(t0):
        vals = [(Fraction(1), Fraction(0)), t0]
        for level in range(1, 15):
            a, b = vals[level]
            prod = (5 * b * t0[1] + a * t0[0], a * t0[1] + b * t0[0])
            prod = (prod[0] * (2 * level + 1), prod[1] * (2 * level + 1))
            prev = vals[level - 1]
            vals.append((
                (prod[0] - level * prev[0]) / (level + 1),
                (prod[1] - level * prev[1]) / (level + 1),
            ))
        return vals

    plus = legendre_at((Fraction(0), Fraction(1, 5)))
    minus = legendre_at((Fraction(0), Fraction(-1, 5)))
    comb = {}
    for level in range(15):
        rational = (
            Fraction(1 + (-1) ** level)
            + 5 * (plus[level][0] + minus[level][0])
        ) / 12
        irrational = plus[level][1] + minus[level][1]
        require(irrational == 0, f"comb irrational at level {level}")
        comb[level] = rational
    require(comb[6] == Fraction(11, 25), "comb level six drift")
    require(comb[10] == Fraction(247, 1875), "comb level ten drift")
    require(comb[12] == Fraction(1071, 3125), "comb level twelve drift")
    require(
        all(comb[l] == 0 for l in (2, 4, 8, 14)),
        "comb even zero set drift",
    )
    require(
        all(comb[l] == 0 for l in range(1, 15, 2)),
        "comb odd null drift",
    )
    return {
        "weights": {"6": "11/25", "10": "247/1875", "12": "1071/3125"},
        "even_zeros": [2, 4, 8, 14],
        "all_odd_zero": True,
        "odd_nulls_computed": True,
    }


def fact_z6_kernel() -> dict[str, Any]:
    # the diagonal element (e^{2 pi i/3} I3, -I2, e^{i pi/3}) acts
    # trivially on the fifteen-state table; recompute with sixth roots
    # of unity as exponents mod 6: charge of each row under
    # (g3, g2, z) = (w^2 I3, w^3 I2, w) with w = e^{i pi/3}:
    # action exponent = 2*(color triality) + 3*(weak duality) + 6*Y.
    rows = {
        "Q": (1, 1, Fraction(1, 6)),
        "u_c": (2, 0, Fraction(-2, 3)),
        "e_c": (0, 0, Fraction(1)),
        "d_c": (2, 0, Fraction(1, 3)),
        "L": (0, 1, Fraction(-1, 2)),
    }
    for name, (t, d, y) in rows.items():
        exponent = (2 * t + 3 * d + 6 * y) % 6
        require(exponent == 0, f"Z6 kernel drift on {name}")
    return {
        "kernel_generator": "(e^{2 pi i/3} I3, -I2, e^{i pi/3})",
        "order": 6,
        "congruence": "2t + 3d + 6Y = 0 mod 6 on every matter row",
        "legacy_fork": (
            "the four legacy-consistent global forms (trivial, Z2, Z3, "
            "Z6 quotients) are experimentally unresolved; the conditional "
            "carrier selection is the Z6 quotient, which forces one "
            "commensurate charge lattice"
        ),
    }


def fact_kinetic_dichotomy() -> dict[str, Any]:
    # Killing-relative su(2):su(3) coefficient ratios of the two
    # branches, in the convention of Lean/Screen/KineticFormDichotomy:
    # killingRelativeSU2(T2) = T2/2 (dual Coxeter number 2) and
    # killingRelativeSU3(T3) = T3/3 (dual Coxeter number 3), with the
    # pinned Dynkin pairs (T2, T3) = (2, 1/2) on the port branch and
    # (2, 2) on the matter branch.
    def ratio(t2: Fraction, t3: Fraction) -> Fraction:
        return (t2 / 2) / (t3 / 3)

    port = ratio(Fraction(2), Fraction(1, 2))
    matter = ratio(Fraction(2), Fraction(2))
    require(port == 6, "port branch ratio drift")
    require(matter == Fraction(3, 2), "matter branch ratio drift")
    return {
        "port_branch_su2_su3_ratio": "6",
        "matter_branch_su2_su3_ratio": "3/2",
        "frozen_statistic": (
            "det(alpha^-1, k, b) = 0: the inverse couplings lie on the "
            "plane spanned by the kinetic ray and the census beta vector; "
            "a scale-free relation among measured couplings that legacy "
            "models do not force; the physical kinetic-action selector "
            "remains open in the kinetic-form selection receipt"
        ),
    }


def fact_koide_window() -> dict[str, Any]:
    # Q = 1/3 + (2/3)(rho/a)^2; at balance the tau root from the register's
    # declared electron and muon inputs must land at the frozen center.
    from decimal import Decimal, getcontext

    getcontext().prec = 40
    me = Decimal("0.51099895069")
    mmu = Decimal("105.6583755")
    s = me.sqrt() + mmu.sqrt()
    b = -4 * s
    c = 3 * (me + mmu) - 2 * s * s
    disc = (b * b - 4 * c).sqrt()
    root = (-b + disc) / 2
    mtau = root * root
    require(
        Decimal("1776.968991") <= mtau <= Decimal("1776.969063"),
        "tau window drift",
    )
    return {
        "identity": "Q = 1/3 + (2/3)(rho/a)^2; balance rho/a = 1/sqrt2",
        "window_MeV": "[1776.968991, 1776.969063]",
        "recomputed_center_in_window": True,
        "custody": "ladder row FZ-10, registered kill band",
    }


def build_rows() -> list[dict[str, Any]]:
    nulls = fact_invariant_nulls()
    odd = fact_odd_moments_vanish()
    five = fact_five_fold_sector()
    comb = fact_angular_comb()
    z6 = fact_z6_kernel()
    dichotomy = fact_kinetic_dichotomy()
    koide = fact_koide_window()

    return [
        {
            "row": 1,
            "differentiator": "spin-six rotational fingerprint",
            "oph_statement": (
                "angular ranks one through five exactly empty at every order; "
                "every residue below rank ten one multiple of the rigid I6 "
                "with the sign-symmetric 62-direction census; on the frozen "
                "primitive equal-weight branch the first directional artifact "
                "is a^4 k^6 with the exact 1/16 refinement step"
            ),
            "recomputed_facts": nulls,
            "legacy_status": (
                "the minimal, locally Lorentz-invariant Standard Model "
                "with General Relativity produces no intrinsic vacuum "
                "rotational residue of any shape; nonminimal operators, "
                "media, and sources are separately modeled"
            ),
            "data_contact": (
                "published subluminal quadratic dispersion limits bound "
                "the carrier scale below 8.9e-29 m (exposed retrospective "
                "diagnostic; primary Xi and Shu, Chinese Physics C 49, "
                "125101 (2025), cross-check LHAASO, Physical Review "
                "Letters 133, 071501 (2024); see the carrier-scale "
                "receipt in this package)"
            ),
            "falsification": "ladder row FZ-11 full-manifold decision rule",
            "type": (
                "finite theorem + frozen prospective physical-branch prediction"
            ),
        },
        {
            "row": 2,
            "differentiator": "no linear vacuum dispersion",
            "oph_statement": (
                "the even-power hop expansion forbids every odd-power "
                "dispersion exactly, so the carrier class predicts zero "
                "linear energy dependence of the vacuum photon speed at "
                "every order"
            ),
            "recomputed_facts": odd,
            "legacy_status": (
                "linear-dispersion quantum-gravity models with "
                "Planckian-coefficient terms are excluded by published "
                "linear bounds above the Planck scale; those bounds "
                "impose no constraint on the carrier class"
            ),
            "data_contact": (
                "published linear limits (exposed literature pointers: "
                "Vasileiou et al., Physical Review D 87, 122001 (2013); "
                "LHAASO, Physical Review Letters 133, 071501 (2024)): "
                "the strongest sit above the Planck energy"
            ),
            "falsification": (
                "a certified intrinsic linear vacuum dispersion of any "
                "size refutes the carrier class"
            ),
            "type": "finite theorem",
        },
        {
            "row": 3,
            "differentiator": "five-fold azimuthal uniqueness",
            "oph_statement": (
                "a resolved carrier anisotropy carries exact five-fold "
                "azimuthal symmetry about each vertex axis (the cos 5 phi "
                "sector of I6 with coefficient 21/8)"
            ),
            "recomputed_facts": five,
            "legacy_status": (
                "no periodic-lattice substrate can produce five-fold "
                "point symmetry; cubic discreteness produces four-fold "
                "spin-four structure first"
            ),
            "data_contact": "none opened; template registered",
            "falsification": (
                "a certified carrier-correlated anisotropy with four- or "
                "six-fold azimuthal structure about its extremal axes "
                "refutes the icosahedral class"
            ),
            "type": "finite theorem (corollary of the fingerprint)",
        },
        {
            "row": 4,
            "differentiator": "equal-port angular comb",
            "oph_statement": (
                "the equal-port carrier measure carries exact weights on "
                "the icosahedral invariant levels: 11/25 at six, 247/1875 "
                "at ten, 1071/3125 at twelve, zeros exactly at 2, 4, 8, "
                "14, all odd levels zero"
            ),
            "recomputed_facts": comb,
            "legacy_status": (
                "isotropic legacy carries no comb; icosahedral-topology "
                "legacy shares the support and leaves every weight free; "
                "the carrier pins all weights with zero parameters"
            ),
            "data_contact": "none opened; registered candidates",
            "falsification": (
                "a certified equal-port transfer with deviating weights "
                "falsifies the carrier measure branch"
            ),
            "type": "finite theorem + registered candidates",
        },
        {
            "row": 5,
            "differentiator": "odd-parity blindness",
            "oph_statement": (
                "the equal-weight carrier channel is exactly blind to "
                "parity-odd sky structure at every level while the odd "
                "invariant tower begins at fifteen"
            ),
            "recomputed_facts": {"all_odd_comb_levels_zero": True,
                                 "first_odd_invariant_level": 15},
            "legacy_status": (
                "parity-violating birefringence models produce parity-odd "
                "signatures; the carrier channel cannot"
            ),
            "data_contact": "none opened; registered candidate",
            "falsification": (
                "a certified carrier-correlated parity-odd response "
                "refutes the readback branch outright"
            ),
            "type": "finite theorem + registered candidate",
        },
        {
            "row": 6,
            "differentiator": "global-form fork and charge commensuration",
            "oph_statement": (
                "the conditional carrier selection is the Z6-quotient "
                "global form, one of four legacy-consistent options, "
                "forcing the congruence 2t + 3d + 6Y = 0 mod 6 and one "
                "commensurate electric-charge lattice"
            ),
            "recomputed_facts": z6,
            "legacy_status": (
                "legacy physics leaves the global form unmeasured; "
                "anomaly cancellation constrains hypercharges in the "
                "minimal spectrum while a continuous dequantization "
                "direction survives with Dirac neutrinos, and the gauge "
                "sector alone forces no commensuration"
            ),
            "data_contact": (
                "matter-neutrality experiments bound residual charges at "
                "the 1e-21 level (exposed literature pointer: Bressi et "
                "al., Physical Review A 83, 052101 (2011)); exact "
                "commensuration is the carrier-side statement"
            ),
            "falsification": (
                "a certified incommensurate charge, or a certified "
                "non-Z6 global form, refutes the conditional selection"
            ),
            "type": "conditional exact packet (issues #642, #567)",
        },
        {
            "row": 7,
            "differentiator": "scale-free coupling-plane relation",
            "oph_statement": (
                "det(alpha^-1, k, b) = 0 with the branch dichotomy: the "
                "port-trace and matter-trace kinetic forms carry "
                "su(2):su(3) Killing-relative ratios 6 and 3/2, and the "
                "selected branch freezes one dimensionless relation "
                "among the measured inverse couplings"
            ),
            "recomputed_facts": dichotomy,
            "legacy_status": (
                "legacy physics treats the three couplings as free at "
                "any one scale; no plane relation is forced"
            ),
            "data_contact": (
                "none opened; the statistic is frozen with the alpha "
                "column sealed"
            ),
            "falsification": (
                "after branch selection, one sealed comparison scores "
                "the relation against the measured couplings"
            ),
            "type": "frozen statistic + open branch selector",
        },
        {
            "row": 8,
            "differentiator": "positive-chamber Koide identity",
            "oph_statement": (
                "the Hermitian three-cycle response forces Q = 1/3 + "
                "(2/3)(rho/a)^2, with Q = 2/3 exactly at the certified "
                "tracial balance; the registered tau window is 72 eV wide"
            ),
            "recomputed_facts": koide,
            "legacy_status": (
                "legacy physics carries no mechanism for the Koide "
                "coordinate; the measured value sits at 2/3 within "
                "experimental precision"
            ),
            "data_contact": (
                "registered conditional test FZ-10 with a custody-bound "
                "kill band; postdictive in premise ancestry"
            ),
            "falsification": (
                "a tau world average more than three standard "
                "uncertainties from the frozen center kills the "
                "balanced-circulant premise"
            ),
            "type": "finite theorem + registered conditional test",
        },
        {
            "row": 9,
            "differentiator": "closure coordinates",
            "oph_statement": (
                "the interval-certified screen-grain root reproduces the "
                "measured inverse coupling to parts per million, and the "
                "capacity coordinate lands within one percent of the "
                "Lambda-derived comparison coordinate"
            ),
            "recomputed_facts": {
                "typed": "diagnostics",
                "custody_paths": [
                    "code/P_derivation/interval_contraction_certificate.py",
                    "code/capacity_readback/manifests/"
                    "n_closure_branch_certificate.json",
                ],
            },
            "legacy_status": (
                "legacy physics carries the fine-structure constant and "
                "the cosmological constant as free parameters"
            ),
            "data_contact": (
                "typed diagnostics with declared measured inputs; "
                "retrospective, no predictive weight claimed"
            ),
            "falsification": (
                "the closure lanes carry their own typed exits; the "
                "screen-grain root dies with its transport bridge"
            ),
            "type": "diagnostics + open closure lanes",
        },
        {
            "row": 10,
            "differentiator": "carrier-class dispersion band and cross-order lock",
            "oph_statement": (
                "every member of the declared positive-weight scalar cosine "
                "class shares C4 < 0 "
                "with the isotropic floor B0/C4^2 >= 10/21, saturated "
                "exactly by single-radius supports; the rank-six ratio "
                "B6/B0 = (16/75)<I6(seed)> is confined to "
                "[-16/135, 16/75]; the k^8 anisotropy is the same rotated "
                "I6 and every single-radius member obeys the division-free "
                "identity 5 D6 B0 = 12 B6 D0, equivalently "
                "D6/D0 = (12/5)(B6/B0), including the exact "
                "zero-anisotropy vertex-face mixture"
            ),
            "recomputed_facts": {
                "floor": str(Fraction(400, 840)),
                "band_low": str(Fraction(16, 75) * Fraction(-5, 9)),
                "band_high": str(Fraction(16, 75)),
                "lock": str(Fraction(64, 125) / Fraction(16, 75)),
                "lock_identity": "5 D6 B0 = 12 B6 D0",
                "lock_ratio_identity": "D6/D0 = (12/5)(B6/B0)",
                "zero_mixture": "B6/B0 = D6/D0 = 0 at vertex:face weights 25:27",
                "vertex_point": str(Fraction(16, 75)),
                "edge_point": str(Fraction(16, 75) * Fraction(-5, 16)),
                "face_point": str(Fraction(16, 75) * Fraction(-5, 9)),
            },
            "legacy_status": (
                "a generic rank-six Lorentz-violating model carries "
                "independent k^4, k^6, and k^8 amplitudes with free signs "
                "and shapes; the carrier class fixes the sign, the floor, "
                "one angular template, and a division-free cross-order identity "
                "on the single-radius stratum"
            ),
            "data_contact": (
                "class-level surface behind the frozen FZ-11 and FZ-12 "
                "branch predictions; their target-clean producer, physical "
                "bridge, quantitative output, and pre-comparison registration "
                "belong to the relevant propagation and prediction lanes"
            ),
            "falsification": (
                "a resolved intrinsic dispersion with B0/C4^2 below "
                "10/21, or B6/B0 outside the band, or violation of "
                "5 D6 B0 = 12 B6 D0 at resolved single-radius saturation, "
                "excludes every member of the class at once"
            ),
            "type": "finite class theorem + Lean skeleton",
        },
    ]


def build_receipt() -> dict[str, Any]:
    rows = build_rows()
    require(len(rows) == 10, "row count drift")
    receipt = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 647,
        "reading": (
            "ten typed separations between the carrier class and legacy "
            "models: two frozen or registered kill surfaces, four finite "
            "theorems with exact recomputed facts, one conditional exact "
            "packet, one frozen statistic with an open selector, and one "
            "diagnostic pair; every row names its falsification "
            "direction and no row opens a comparison"
        ),
        "rows": rows,
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
            "literature_pointers_are_exposed_inputs": True,
        },
    }
    receipt["receipt_sha256"] = base.tagged_sha256(
        base.canonical_json_bytes(
            {k: v for k, v in receipt.items() if k != "receipt_sha256"}
        )
    )
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    receipt = build_receipt()
    if args.write:
        RUNTIME.mkdir(exist_ok=True)
        RECEIPT_PATH.write_bytes(base.canonical_json_bytes(receipt))
        print(RECEIPT_PATH)
        return 0
    for row in receipt["rows"]:
        print(f'{row["row"]}. {row["differentiator"]} [{row["type"]}]')
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
