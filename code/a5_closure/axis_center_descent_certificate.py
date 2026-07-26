#!/usr/bin/env python3
"""Exact certificate for GitHub issue #567: AXIS-CENTER-DESCENT.

The input is the hash-pinned #314 matter-lift receipt, whose kernel emission
carries the common action kernel on the simply connected cover as data. This
certificate consumes that data and derives, rather than assumes:

* the common kernel of the central action on every realized matter tensor,
  computed by exhaustive enumeration over the 36 candidate central elements
  rather than inferred from Lie type;
* the maximal effective image
  (SU(3) x SU(2) x U(1)) / Z6 of the declared tensor representation,
  together with its character lattice and dual cocharacter lattice by exact
  integer/rational arithmetic;
* the primitive correlated cocharacter, in color, weak, and electron-Dirac
  units under q = 6Y with
  h = (omega_3, -I_2, exp(i pi/3)); this is lattice arithmetic, not a
  physical choice of genuine Wilson/'t Hooft lines or a theta-period theorem;
* the exact four-way non-identifiability result: quotienting by
  1, Z2, Z3, or Z6 gives four different character lattices while every
  declared local tensor descends through all four;
* the algebraic weak-center/U(1) relation h^3 = (1,-I_2,-1), which acts
  trivially on every declared weight.  It is not fermion parity.  Indeed no
  element among the 36 central candidates acts by -1 on all five matter
  multiplets, so no spacetime Spin attachment follows from this calculation;
* weight-level refinement invariance under the sixty carrier rotations and
  the artifact persistence maps;
* countermodels: an adjoint-only tensor set leaves all four global forms
  admissible, and a fractionally charged extra tensor shrinks the kernel.

The source packet contains no physical loop/deck holonomy, tangent-frame Spin
lift, genuine-line category or UV polarization, principal-bundle/instanton
sector, or refinement map on those structures.  The physical #567 gate
therefore fails closed.

Every arithmetic decision is exact integer or rational arithmetic; no
floating point appears in a proof step.
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import echosahedral_selector_certificate as e565  # noqa: E402

CertificateError = e565.CertificateError
require = e565.require
sha256_json = e565.sha256_json
load_json = e565.load_json
write_json = e565.write_json

SCHEMA = "oph.axis_center_descent_manifest.v3"
RECEIPT_SCHEMA = "oph.axis_center_descent_receipt.v3"
NEGATIVE_SCHEMA = "oph.axis_center_descent_negative_controls.v3"

# Realized weight table under the corpus convention q = 6Y: every realized
# matter field and both carrier blocks, each with (q, triality, duality). The
# scalar is deliberately unnecessary: its weak-block weight duplicates an
# already present character, and scalar existence/economy is not source-bound.
REALIZED_WEIGHTS: dict[str, tuple[int, int, int]] = {
    "Q": (1, 1, 1),
    "u_c": (-4, 2, 0),
    "e_c": (6, 0, 0),
    "d_c": (2, 2, 0),
    "L": (-3, 0, 1),
    "carrier_color": (-2, 1, 0),
    "carrier_weak": (3, 0, 1),
}

MATTER_LABELS = ("Q", "u_c", "e_c", "d_c", "L")

ADJOINT_WEIGHTS: dict[str, tuple[int, int, int]] = {
    "gluon": (0, 0, 0),
    "weak_boson": (0, 0, 0),
    "hypercharge_boson": (0, 0, 0),
}


def central_phase_sixths(element: tuple[int, int, int], weight: tuple[int, int, int]) -> int:
    """The phase of h^element on a weight, in sixths of a turn modulo six.

    The central element (k, l, r) acts on a tensor with data (q, t, d) by
    exp(2 pi i (k t / 3 + l d / 2 + r q / 6)); in sixths of a turn this is
    2 k t + 3 l d + r q modulo six.
    """

    k, l, r = element
    q, t, d = weight
    return (2 * k * t + 3 * l * d + r * q) % 6


def common_kernel(weights: Mapping[str, tuple[int, int, int]]) -> list[tuple[int, int, int]]:
    """Exhaustive enumeration of the common central kernel."""

    kernel = []
    for k in range(3):
        for l in range(2):
            for r in range(6):
                if all(
                    central_phase_sixths((k, l, r), weight) == 0
                    for weight in weights.values()
                ):
                    kernel.append((k, l, r))
    return kernel


def cyclic_generated(generator: tuple[int, int, int]) -> set[tuple[int, int, int]]:
    k, l, r = generator
    return {((k * n) % 3, (l * n) % 2, (r * n) % 6) for n in range(6)}


def central_power(
    generator: tuple[int, int, int], exponent: int
) -> tuple[int, int, int]:
    """A power of a central parameter tuple."""

    k, l, r = generator
    return ((k * exponent) % 3, (l * exponent) % 2, (r * exponent) % 6)


def common_scalar_phase_elements(
    weights: Mapping[str, tuple[int, int, int]], phase_sixths: int
) -> list[tuple[int, int, int]]:
    """Central elements acting with one prescribed scalar phase on all weights."""

    return [
        (k, l, r)
        for k in range(3)
        for l in range(2)
        for r in range(6)
        if all(
            central_phase_sixths((k, l, r), weight) == phase_sixths % 6
            for weight in weights.values()
        )
    ]


def smith_normal_form(matrix: list[list[int]]) -> list[list[int]]:
    """Smith normal form of an integer matrix (small sizes only)."""

    m = [row[:] for row in matrix]
    rows, cols = len(m), len(m[0])

    def find_pivot(start: int) -> tuple[int, int] | None:
        best = None
        for i in range(start, rows):
            for j in range(start, cols):
                if m[i][j] != 0 and (best is None or abs(m[i][j]) < abs(m[best[0]][best[1]])):
                    best = (i, j)
        return best

    for position in range(min(rows, cols)):
        while True:
            pivot = find_pivot(position)
            if pivot is None:
                return m
            i, j = pivot
            m[position], m[i] = m[i], m[position]
            for row in m:
                row[position], row[j] = row[j], row[position]
            value = m[position][position]
            reduced = True
            for i in range(position + 1, rows):
                quotient = m[i][position] // value
                if quotient:
                    for j in range(cols):
                        m[i][j] -= quotient * m[position][j]
                if m[i][position]:
                    reduced = False
            for j in range(position + 1, cols):
                quotient = m[position][j] // value
                if quotient:
                    for i in range(rows):
                        m[i][j] -= quotient * m[i][position]
                if m[position][j]:
                    reduced = False
            if reduced:
                if m[position][position] < 0:
                    for j in range(cols):
                        m[position][j] = -m[position][j]
                break
    # Divisibility normalization: adjacent diagonal pairs (a, b) with a not
    # dividing b become (gcd(a, b), lcm(a, b)).
    from math import gcd

    size = min(rows, cols)
    changed = True
    while changed:
        changed = False
        for index in range(size - 1):
            a, b = m[index][index], m[index + 1][index + 1]
            if a != 0 and b % a != 0:
                g = gcd(a, b)
                m[index][index] = g
                m[index + 1][index + 1] = a * b // g
                changed = True
    return m


def subgroups_of_z6() -> list[dict[str, Any]]:
    """All subgroups of the cyclic centre Z6, enumerated by element order."""

    elements = list(range(6))
    subgroups = []
    for divisor in (1, 2, 3, 6):
        generator = 6 // divisor
        members = sorted({(generator * n) % 6 for n in range(divisor)})
        subgroups.append({"order": divisor, "members": members})
    # Verify these are exactly the subsets closed under addition mod six.
    closed = []
    for mask in range(1 << 6):
        subset = [x for x in elements if mask & (1 << x)]
        if 0 not in subset or not subset:
            continue
        if all((a + b) % 6 in subset for a in subset for b in subset):
            closed.append(sorted(subset))
    require(
        sorted(closed) == sorted(s["members"] for s in subgroups),
        "GLOBAL_FORM_MENU",
        "the subgroup enumeration of the centre drifted",
    )
    return subgroups


def quotient_candidate_menu(
    weights: Mapping[str, tuple[int, int, int]],
    generator: tuple[int, int, int],
) -> list[dict[str, Any]]:
    """The four diagonal quotient choices and their local-data fingerprints.

    A subgroup of order ``d`` is generated by ``h^(6/d)``.  Its quotient
    character residues satisfy ``2t + 3d_w + q = 0 (mod d)``.  Every subgroup
    lies in the common kernel, hence every declared local weight descends
    through every one of these quotients.  The different residue counts show
    that their global character data is nevertheless different.
    """

    rows: list[dict[str, Any]] = []
    for order in (1, 2, 3, 6):
        exponent = 6 // order
        subgroup_generator = central_power(generator, exponent)
        subgroup_elements = {
            central_power(generator, exponent * n) for n in range(order)
        }
        require(
            all(
                central_phase_sixths(element, weight) == 0
                for element in subgroup_elements
                for weight in weights.values()
            ),
            "GLOBAL_FORM_MENU",
            f"a realized tensor does not descend through the order-{order} subgroup",
        )
        character_residues = [
            (t, d, q)
            for t in range(3)
            for d in range(2)
            for q in range(6)
            if (2 * t + 3 * d + q) % order == 0
        ]
        require(
            len(character_residues) == 36 // order,
            "GLOBAL_FORM_MENU",
            f"the order-{order} quotient character-residue count drifted",
        )
        rows.append(
            {
                "quotient_subgroup_order": order,
                "quotient_label": (
                    "SU(3) x SU(2) x U(1)"
                    if order == 1
                    else f"(SU(3) x SU(2) x U(1)) / Z{order}"
                ),
                "generator_power_of_h": 0 if order == 1 else exponent,
                "generator_central_parameters": list(subgroup_generator),
                "character_residue_class_count": len(character_residues),
                "all_declared_local_tensors_descend": True,
            }
        )
    require(
        len({row["character_residue_class_count"] for row in rows}) == 4,
        "GLOBAL_FORM_MENU",
        "the four quotient candidates were not distinguished by global character data",
    )
    return rows


def load_matter_receipt(manifest: Mapping[str, Any], base_dir: Path) -> dict[str, Any]:
    path_raw = manifest.get("matter_receipt_path")
    require(isinstance(path_raw, str), "UPSTREAM_REFERENCE", "matter_receipt_path is missing")
    path = Path(path_raw)
    if not path.is_absolute():
        path = base_dir / path
    receipt = load_json(path)
    require(
        manifest.get("matter_receipt_sha256") == sha256_json(receipt),
        "UPSTREAM_HASH",
        "the #314 matter receipt hash does not match the declared pin",
    )
    require(
        receipt.get("issue") == 314
        and str(receipt.get("schema", "")).startswith("oph.super_tannakian_matter_receipt"),
        "UPSTREAM_RECEIPT",
        "the pinned receipt is not a #314 matter-lift receipt",
    )
    require(
        receipt.get("conditional_algebraic_gate", {}).get("passed") is True,
        "UPSTREAM_RECEIPT",
        "the pinned matter receipt must record a passing conditional algebraic gate",
    )
    kernel_emission = receipt.get("kernel_emission")
    require(
        isinstance(kernel_emission, Mapping)
        and kernel_emission.get("global_quotient_assumed") is False,
        "UPSTREAM_RECEIPT",
        "the pinned receipt must emit the kernel without assuming the quotient",
    )
    return receipt


def validate_manifest(manifest: Mapping[str, Any]) -> None:
    require(manifest.get("schema") == SCHEMA, "SCHEMA", f"expected {SCHEMA}")
    require(
        manifest.get("hypercharge_convention") == "q = 6Y",
        "CONVENTION",
        "the manifest must declare the integer hypercharge convention q = 6Y",
    )
    for key in (
        "measured_coupling",
        "mass_target",
        "monopole_dynamics",
        "claim_physical_global_form",
        "block_map",
        "physical_deck_loop",
        "spacetime_spin_attachment",
        "genuine_line_category",
        "uv_line_polarization",
        "instanton_sector",
        "theta_periodicity",
    ):
        require(key not in manifest, "FORBIDDEN_DEPENDENCY", f"forbidden manifest key {key}")


def rederive_weights(receipt: Mapping[str, Any]) -> dict[str, tuple[int, int, int]]:
    """Re-derive the integer weight table from the pinned matter receipt."""

    package = receipt.get("realized_package", {})
    fields = package.get("fields", {})
    require(
        set(fields) == {"Q", "u_c", "e_c", "d_c", "L"},
        "WEIGHT_TABLE",
        "the realized package does not carry the five matter fields",
    )
    triality = {"Q": 1, "u_c": 2, "e_c": 0, "d_c": 2, "L": 0}
    duality = {"Q": 1, "u_c": 0, "e_c": 0, "d_c": 0, "L": 1}
    weights: dict[str, tuple[int, int, int]] = {}
    for label, data in fields.items():
        charge = Fraction(data.get("charge"))
        q = 6 * charge
        require(
            q.denominator == 1,
            "WEIGHT_TABLE",
            f"field {label} has a non-integral q = 6Y charge",
        )
        weights[label] = (int(q), triality[label], duality[label])
    carrier = receipt.get("matter_carrier", {})
    charges = carrier.get("block_trace_charges", {})
    q_color = 6 * Fraction(charges.get("color_block"))
    q_weak = 6 * Fraction(charges.get("weak_block"))
    require(
        q_color.denominator == 1 and q_weak.denominator == 1,
        "WEIGHT_TABLE",
        "the carrier block charges are not integral in q = 6Y",
    )
    weights["carrier_color"] = (int(q_color), 1, 0)
    weights["carrier_weak"] = (int(q_weak), 0, 1)
    require(
        weights == REALIZED_WEIGHTS,
        "WEIGHT_TABLE",
        f"the re-derived weight table drifted: {weights}",
    )
    return weights


def certificate_payload(manifest: Mapping[str, Any], base_dir: Path | None = None) -> dict[str, Any]:
    base = base_dir or MODULE_DIR
    validate_manifest(manifest)
    receipt = load_matter_receipt(manifest, base)
    weights = rederive_weights(receipt)

    # --- The kernel on every realized matter tensor -------------------------
    kernel = common_kernel(weights)
    require(len(kernel) == 6, "KERNEL", f"expected a six-element kernel, got {len(kernel)}")
    generator = (1, 1, 1)
    require(
        set(kernel) == cyclic_generated(generator),
        "KERNEL",
        "the kernel is not the cyclic group generated by (omega_3, -I_2, zeta_6)",
    )
    emission = receipt["kernel_emission"]
    emitted_generator = emission.get("kernel_generator", {})
    require(
        emitted_generator.get("su3_center_power") == 1
        and emitted_generator.get("su2_center_power") == 1
        and Fraction(emitted_generator.get("u1_phase_turns")) == Fraction(1, 6),
        "KERNEL",
        "the recomputed kernel generator does not match the emitted kernel data",
    )
    # Tensor additivity: the kernel condition is additive in (q, t, d), so
    # triviality on the verified weight list extends to every realized tensor.
    for first in weights.values():
        for second in weights.values():
            summed = (
                first[0] + second[0],
                (first[1] + second[1]) % 3,
                (first[2] + second[2]) % 2,
            )
            for element in kernel:
                require(
                    central_phase_sixths(element, summed) == 0,
                    "KERNEL",
                    "the kernel fails on a tensor product of realized weights",
                )

    # --- Character and cocharacter lattices by Smith reduction ---------------
    # Characters of the cover torus centre data are labeled (t, d, q) in Z^3;
    # characters of the quotient are those trivial on the kernel:
    # 2t + 3d + q = 0 mod 6. The explicit basis (1,0,-2), (0,1,-3), (0,0,6)
    # of that sublattice has Smith invariants (1, 1, 6): the quotient
    # character lattice has index six with cyclic quotient Z6, and dually the
    # quotient cocharacter lattice contains the cover lattice with index six.
    reduced = smith_normal_form([[1, 0, -2], [0, 1, -3], [0, 0, 6]])
    invariants = [reduced[i][i] for i in range(3)]
    require(
        invariants == [1, 1, 6],
        "LATTICE",
        f"the character sublattice has Smith invariants {invariants}",
    )
    allowed_classes = [
        (t, d, q)
        for t in range(3)
        for d in range(2)
        for q in range(6)
        if (2 * t + 3 * d + q) % 6 == 0
    ]
    require(
        len(allowed_classes) == 6,
        "LATTICE",
        "the maximal-quotient character residue table does not have six elements",
    )
    for label, weight in weights.items():
        q, t, d = weight
        require(
            (2 * t + 3 * d + q) % 6 == 0,
            "LATTICE",
            f"realized weight {label} is not a quotient character",
        )

    # --- Dual cocharacter lattice and its primitive correlated class ----------
    # The quotient character lattice is the kernel of (t, d, q) -> 2t + 3d + q
    # modulo six inside Z^3; the vectors (1, 0, -2), (0, 1, -3), (0, 0, 6)
    # form a basis (each satisfies the constraint and their determinant is
    # six, the index of the sublattice). A magnetic cocharacter (m, n, g)
    # with color centre flux m/3, weak flux n/2, and U(1) charge g pairs
    # integrally with the whole lattice exactly when it pairs integrally
    # with this basis.
    basis = [(1, 0, -2), (0, 1, -3), (0, 0, 6)]
    for t_basis, d_basis, q_basis in basis:
        require(
            (2 * t_basis + 3 * d_basis + q_basis) % 6 == 0,
            "LATTICE",
            "a declared character-lattice basis vector violates the constraint",
        )
    determinant = basis[0][0] * (basis[1][1] * basis[2][2] - basis[1][2] * basis[2][1])
    require(determinant == 6, "LATTICE", f"the basis determinant is {determinant}, not six")

    def admissible_magnetic(m: int, n: int, g_sixths: Fraction) -> bool:
        pair_one = Fraction(m, 3) - 2 * g_sixths
        pair_two = Fraction(n, 2) - 3 * g_sixths
        pair_three = 6 * g_sixths
        return (
            pair_one.denominator == 1
            and pair_two.denominator == 1
            and pair_three.denominator == 1
        )

    cocharacter_window = [
        (m, n, Fraction(g, 6))
        for m in range(3)
        for n in range(2)
        for g in range(0, 36)
        if admissible_magnetic(m, n, Fraction(g, 6))
    ]
    nontrivial = [row for row in cocharacter_window if row != (0, 0, Fraction(0))]
    minimal = min(nontrivial, key=lambda row: (row[2], row[0], row[1]))
    require(
        minimal == (1, 1, Fraction(1, 6)),
        "COCHARACTER_LATTICE",
        f"the primitive correlated cocharacter drifted: {minimal}",
    )
    pure_u1 = [
        row
        for row in cocharacter_window
        if row[0] == 0 and row[1] == 0 and row[2] != 0
    ]
    require(
        min(row[2] for row in pure_u1) == 1,
        "COCHARACTER_LATTICE",
        "the primitive pure-U(1) cocharacter is not one electron-Dirac quantum",
    )
    # If B is the displayed character basis, its dual lattice is B^{-1} Z^3.
    # These three exact columns generate it.
    dual_basis = [
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(1, 3), Fraction(1, 2), Fraction(1, 6)),
    ]
    for cocharacter in dual_basis:
        for character in basis:
            pairing = sum(
                Fraction(character[index]) * cocharacter[index]
                for index in range(3)
            )
            require(
                pairing.denominator == 1,
                "COCHARACTER_LATTICE",
                "the declared dual-basis pairing is not integral",
            )

    # --- Weak-centre relation; no spacetime Spin conclusion ------------------
    # h^3 = (1, -I_2, -1) is a gauge-kernel relation.  It acts as +1, not
    # fermion parity, on every realized tensor.  Exhaustive enumeration also
    # shows that no central element acts by -1 on all five matter multiplets.
    weak_center_u1_relation = central_power(generator, 3)
    require(
        weak_center_u1_relation == (0, 1, 3)
        and weak_center_u1_relation in set(kernel),
        "CENTRAL_RELATION",
        "h^3 is not the expected weak-centre/U(1) kernel relation",
    )
    matter_weights = {label: weights[label] for label in MATTER_LABELS}
    fermion_minus_one_candidates = common_scalar_phase_elements(matter_weights, 3)
    require(
        fermion_minus_one_candidates == [],
        "FERMION_PARITY",
        "a central gauge element unexpectedly acts by -1 on every matter multiplet",
    )

    # Refinement invariance: the weight data (q, t, d) is a class function of
    # the realized fields; the sixty carrier rotations and the artifact
    # persistence maps permute states inside each field block (verified in the
    # pinned #314 receipt), so the kernel computation is refinement invariant.
    physical_maps = receipt.get("refinement", {}).get("physical_maps", [])
    require(
        len(physical_maps) >= 2 and all(row.get("intertwined") for row in physical_maps),
        "REFINEMENT_INVARIANCE",
        "the pinned receipt does not carry intertwined physical refinement maps",
    )

    # --- Four locally indistinguishable global forms -------------------------
    menu = subgroups_of_z6()
    quotient_candidates = quotient_candidate_menu(weights, generator)
    adjoint_kernel = common_kernel(ADJOINT_WEIGHTS)
    require(
        len(adjoint_kernel) == 36,
        "GLOBAL_FORM_MENU",
        "the adjoint-only kernel is not the full centre",
    )
    # All four subgroups of the common diagonal kernel are compatible with
    # the realized local tensors.  The full-kernel quotient is the maximal
    # effective image, but selecting it as the physical gauge group requires
    # additional global/line data.
    fractional_weights = dict(REALIZED_WEIGHTS)
    fractional_weights["fractional_singlet"] = (1, 0, 0)
    fractional_kernel = common_kernel(fractional_weights)
    require(
        len(fractional_kernel) == 1,
        "GLOBAL_FORM_MENU",
        "the fractionally charged countermodel does not shrink the kernel to the identity",
    )

    return {
        "schema": RECEIPT_SCHEMA,
        "issue": 567,
        "manifest_sha256": sha256_json(manifest),
        "matter_receipt_sha256": sha256_json(receipt),
        "convention": {
            "hypercharge": "q = 6Y",
            "kernel_generator": "h = (omega_3 I_3, -I_2, exp(i pi/3))",
            "overall_charge_sign": (
                "the displayed generator uses the even-Weyl convention; charge "
                "conjugation inverts it and leaves the generated Z6 subgroup unchanged"
            ),
        },
        "realized_weight_table": {
            label: {"q": q, "triality": t, "duality": d}
            for label, (q, t, d) in sorted(weights.items())
        },
        "kernel_on_realized_tensors": {
            "candidates_enumerated": 36,
            "kernel_order": 6,
            "kernel_elements": [list(row) for row in sorted(kernel)],
            "cyclic_generator": [1, 1, 1],
            "matches_emitted_kernel_data": True,
            "tensor_additivity_checked_pairs": len(weights) ** 2,
            "computed_not_inferred": "the kernel is the exhaustive common stabilizer of the realized weight table, not a Lie-type inference",
        },
        "maximal_effective_image": {
            "group": "(SU(3) x SU(2) x U(1)) / Z6",
            "status": (
                "maximal effective image of the declared tensor representation; "
                "not a selected physical global gauge group"
            ),
            "character_lattice_smith_invariants": invariants,
            "character_residue_classes": [list(row) for row in allowed_classes],
            "character_residue_class_count": 6,
        },
        "dual_cocharacter_lattice": {
            "basis_coordinates_color_weak_u1": [
                [str(value) for value in row] for row in dual_basis
            ],
            "primitive_correlated_cocharacter": {
                "color_coweight": "1/3",
                "weak_coweight": "1/2",
                "u1_charge_electron_dirac_units": "1/6",
            },
            "primitive_pure_u1_cocharacter": "one electron-Dirac quantum",
            "interpretation_boundary": (
                "the dual lattice and its finite central-flux residues are exact; "
                "a physical spectrum of genuine Wilson, 't Hooft, or dyonic lines "
                "requires a line category and mutually-local UV polarization"
            ),
            "theta_periodicity_status": (
                "not derived: this packet contains no principal-bundle sectors, "
                "topological action normalization, or instanton-charge lattice"
            ),
            "monopole_dynamics_not_inferred": True,
        },
        "algebraic_weak_center_u1_relation": {
            "h_cubed": list(weak_center_u1_relation),
            "phase_on_every_declared_weight_sixths": 0,
            "universal_fermion_minus_one_candidates": [
                list(row) for row in fermion_minus_one_candidates
            ],
            "conclusion": (
                "h^3 is a gauge-kernel relation acting trivially on the declared "
                "weights; it is not fermion parity, and the finite center table "
                "does not attach the matter object to a spacetime Spin bundle"
            ),
        },
        "weight_level_refinement_invariance": {
            "carrier_rotations": 60,
            "artifact_persistence_maps": len(physical_maps),
            "conclusion": (
                "the central weight table and common-kernel computation are "
                "invariant under the declared state permutations"
            ),
            "physical_loop_or_bundle_refinement_naturality_derived": False,
        },
        "global_form_nonidentifiability": {
            "subgroup_menu": menu,
            "quotient_candidates": quotient_candidates,
            "subgroup_menu_scope": (
                "exactly the four subgroups of the selected diagonal cyclic Z6; "
                "not every subgroup of the full adjoint centre Z3 x Z2 x Z6"
            ),
            "adjoint_only_kernel_order": 36,
            "conclusion_menu": (
                "the diagonal quotient chain has four members: 1, Z2, Z3, and Z6; "
                "all declared local tensors descend through every member"
            ),
            "adjoint_only_scope": (
                "adjoint tensors leave the full 36-element finite central slice "
                "invisible and therefore do not select the diagonal chain"
            ),
            "maximal_effective_quotient": "Z6",
            "physical_global_form_selected": False,
            "selection_boundary": (
                "quotienting by the full common kernel is the maximal effective "
                "local action; choosing it as the physical gauge group requires "
                "source-derived global/line data"
            ),
            "fractional_singlet_countermodel_kernel_order": 1,
            "conclusion_countermodels": (
                "a fractionally charged extra tensor shrinks the common kernel to "
                "the identity, while removing matter enlarges it; local weights "
                "determine the maximal kernel but do not choose which subgroup is gauged"
            ),
        },
        "two_z6_constructions": {
            "six_axis_lattice_quotient": (
                "Lean/Screen/Z6Exact.lean: the abstract coefficient-lattice "
                "quotient Lambda_plus/(Lambda_1 (+) Lambda_5)"
            ),
            "tensor_action_kernel": (
                "Lean/Screen/TraceBalancedKernel.lean: the kernel of "
                "2k + 3l + r on the central parameters"
            ),
            "algebraic_intertwiner": (
                "Lean/Screen/Z6Descent.lean: an exact group isomorphism matching "
                "chosen generators and inversion"
            ),
            "physical_loop_intertwiner_derived": False,
            "boundary": (
                "an isomorphism between two abstract cyclic groups is not a "
                "source-derived port-loop holonomy or deck action"
            ),
        },
        "claim_boundary": {
            "proves": (
                "the common central kernel on the declared tensors, the maximal "
                "effective quotient and its character/cocharacter lattices, the "
                "algebraic intertwiner of the two abstract Z6 groups, weight-level "
                "refinement invariance, and exact four-way global-form non-identifiability"
            ),
            "status": "derived_conditionally_on_unsourced_upstream_response_and_matter_branch",
            "does_not_close": [
                "physical selection/attachment of this global form (#567)",
                "the independent #599 response producer and physical #314 matter source gate",
                "a source-derived physical deck/loop class rather than declared central weight data",
                "a genuine-line category, mutually-local UV polarization, and completeness of the physical line spectrum",
                "spacetime Spin/fermion attachment (h^3 is a trivial gauge relation, not fermion parity)",
                "refinement-natural transport of loops, bundles, Spin data, or line categories",
                "theta periodicity or instanton-sector quantization",
                "monopole or dyon dynamics (only cocharacter-lattice arithmetic is derived)",
                "family attachment and any three-family claim (#569)",
                "laboratory measurement of any line or flux",
                "continuum quantum field theory",
            ],
        },
        "conditional_algebraic_gate": {
            "kernel_enumerated": True,
            "character_and_cocharacter_lattices_computed": True,
            "maximal_effective_image_computed": True,
            "four_global_forms_locally_indistinguishable": True,
            "weak_center_relation_not_misidentified_as_spin": True,
            "passed": True,
        },
        "physical_global_form_gate": {
            "upstream_response_physically_source_bound": False,
            "upstream_matter_physically_source_bound": False,
            "source_derived_deck_loop_class": False,
            "spacetime_spin_attachment": False,
            "genuine_line_category_selected": False,
            "uv_mutual_locality_polarization_selected": False,
            "refinement_natural_loop_bundle_transport": False,
            "instanton_sector_and_action_normalization": False,
            "theta_periodicity_derived": False,
            "laboratory_global_form_attachment": False,
            "passed": False,
        },
        "acceptance_criteria_status": {
            "abstract_z6_intertwiner": True,
            "source_derived_loop_deck_intertwiner": False,
            "kernel_on_every_declared_tensor": True,
            "maximal_effective_image_character_lattice": True,
            "dual_cocharacter_lattice": True,
            "physical_global_quotient_selected": False,
            "genuine_wilson_tHooft_dyonic_line_category_selected": False,
            "theta_periodicity_from_instantiated_topological_sectors": False,
            "primitive_correlated_cocharacter_in_electron_dirac_and_color_weak_units": True,
            "monopole_dynamics_not_inferred": True,
            "algebraic_weak_center_u1_relation": True,
            "spacetime_spin_fermion_attachment": False,
            "weight_level_refinement_invariance": True,
            "loop_bundle_line_refinement_naturality": False,
            "larger_cover_and_smaller_quotient_countermodels": True,
            "four_compatible_cover_level_global_forms": True,
            "issue_567_closeable": False,
        },
        "issue_closure_condition": {
            "conditional_algebraic_gate_passed": True,
            "physical_global_form_gate_passed": False,
            "met_locally": False,
            "remaining_producer": (
                "after #599/#314 source closure: a source-derived loop/deck "
                "holonomy, spacetime Spin attachment, genuine-line category with "
                "UV polarization, refinement-natural transport of that data, and "
                "principal-bundle/instanton normalization sufficient for theta periodicity"
            ),
        },
        "verifier_command": (
            "python3 code/a5_closure/axis_center_descent_certificate.py verify "
            "--manifest code/a5_closure/manifests/axis_center_descent_reference.json "
            "--receipt code/a5_closure/receipts/axis_center_descent_reference.receipt.json"
        ),
    }


def negative_control_cases(manifest: Mapping[str, Any]) -> list[tuple[str, dict[str, Any], str]]:
    cases: list[tuple[str, dict[str, Any], str]] = []

    wrong_convention = copy.deepcopy(dict(manifest))
    wrong_convention["hypercharge_convention"] = "q = Y"
    cases.append(("wrong_hypercharge_convention", wrong_convention, "CONVENTION"))

    physical_promotion = copy.deepcopy(dict(manifest))
    physical_promotion["claim_physical_global_form"] = True
    cases.append(("unsupported_physical_promotion", physical_promotion, "FORBIDDEN_DEPENDENCY"))

    wrong_pin = copy.deepcopy(dict(manifest))
    wrong_pin["matter_receipt_sha256"] = "0" * 64
    cases.append(("wrong_matter_receipt_pin", wrong_pin, "UPSTREAM_HASH"))

    forbidden = copy.deepcopy(dict(manifest))
    forbidden["monopole_dynamics"] = {"target": "monopole mass"}
    cases.append(("monopole_dynamics_injection", forbidden, "FORBIDDEN_DEPENDENCY"))

    for name, key in (
        ("declared_deck_loop_injection", "physical_deck_loop"),
        ("declared_spacetime_spin_injection", "spacetime_spin_attachment"),
        ("declared_line_category_injection", "genuine_line_category"),
        ("declared_uv_polarization_injection", "uv_line_polarization"),
        ("declared_instanton_sector_injection", "instanton_sector"),
        ("declared_theta_period_injection", "theta_periodicity"),
    ):
        mutant = copy.deepcopy(dict(manifest))
        mutant[key] = {"declared_without_source_receipt": True}
        cases.append((name, mutant, "FORBIDDEN_DEPENDENCY"))

    return cases


def negative_control_payload(manifest: Mapping[str, Any], base_dir: Path | None = None) -> dict[str, Any]:
    results = []
    for name, mutant, expected_code in negative_control_cases(manifest):
        actual_code = "ACCEPTED"
        try:
            certificate_payload(mutant, base_dir)
        except CertificateError as exc:
            actual_code = exc.code
        require(
            actual_code == expected_code,
            "NEGATIVE_CONTROL_FAILED",
            f"{name}: expected {expected_code}, got {actual_code}",
        )
        results.append({"name": name, "expected_error": expected_code, "actual_error": actual_code, "passed": True})
    return {
        "schema": NEGATIVE_SCHEMA,
        "issue": 567,
        "manifest_sha256": sha256_json(manifest),
        "finite_controls": results,
        "countermodel_witnesses": {
            "four_local_data_completions": (
                "the same declared local tensor table descends through quotients "
                "by 1, Z2, Z3, and Z6, whose character residue counts are "
                "respectively 36, 18, 12, and 6"
            ),
            "adjoint_only": (
                "the adjoint-only table has the full 36-element finite central "
                "slice in its kernel"
            ),
            "fractional_singlet": "kernel order 1: a larger cover is forced by a fractionally charged tensor",
            "promotion": "a manifest cannot promote the conditional kernel quotient to a physically attached global form",
            "fermion_parity": (
                "no central element in the 36-element table acts by -1 on all "
                "five matter multiplets; h^3 acts trivially"
            ),
        },
    }


def verify_receipt(manifest: Mapping[str, Any], receipt: Mapping[str, Any], base_dir: Path | None = None) -> None:
    expected = certificate_payload(manifest, base_dir)
    require(receipt == expected, "RECEIPT_MISMATCH", "receipt is stale, malformed, or tampered")


def default_paths() -> tuple[Path, Path, Path]:
    return (
        MODULE_DIR / "manifests" / "axis_center_descent_reference.json",
        MODULE_DIR / "receipts" / "axis_center_descent_reference.receipt.json",
        MODULE_DIR / "negative_controls" / "issue_567_negative_controls.json",
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    certify = sub.add_parser("certify")
    certify.add_argument("--manifest", type=Path, required=True)
    certify.add_argument("--output", type=Path, required=True)
    verify = sub.add_parser("verify")
    verify.add_argument("--manifest", type=Path, required=True)
    verify.add_argument("--receipt", type=Path, required=True)
    negative = sub.add_parser("negative-controls")
    negative.add_argument("--manifest", type=Path, required=True)
    negative.add_argument("--output", type=Path, required=True)
    all_cmd = sub.add_parser("all")
    all_cmd.add_argument("--manifest", type=Path)
    args = parser.parse_args(argv)
    if args.command == "certify":
        manifest = load_json(args.manifest)
        receipt = certificate_payload(manifest, args.manifest.resolve().parent.parent)
        write_json(args.output, receipt)
        print(json.dumps({"status": "PASS", "receipt": str(args.output)}, indent=2))
    elif args.command == "verify":
        manifest = load_json(args.manifest)
        receipt = load_json(args.receipt)
        verify_receipt(manifest, receipt, args.manifest.resolve().parent.parent)
        print(json.dumps({"status": "PASS", "receipt": str(args.receipt)}, indent=2))
    elif args.command == "negative-controls":
        manifest = load_json(args.manifest)
        payload = negative_control_payload(manifest, args.manifest.resolve().parent.parent)
        write_json(args.output, payload)
        print(json.dumps({"status": "PASS"}, indent=2))
    else:
        default_manifest, default_receipt, default_negative = default_paths()
        manifest_path = args.manifest or default_manifest
        manifest = load_json(manifest_path)
        write_json(default_receipt, certificate_payload(manifest))
        write_json(default_negative, negative_control_payload(manifest))
        print(
            json.dumps(
                {"status": "PASS", "receipt": str(default_receipt), "negative_controls": str(default_negative)},
                indent=2,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
