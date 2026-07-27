#!/usr/bin/env python3
"""Exact certificate for GitHub issues #616 and #617: the bounded
scalar-response and family-multiplicity campaigns.

The certificate runs two exact campaigns and emits one manifest.

Scalar response multiplicity (#616).  The pinned matter-lift manifest of
issue #314 declares one color-singlet weak-doublet scalar on the weak
block.  The compatibility scan of the matter-lift certificate is rerun on
the pinned block charges: conditional on that declared scalar, the
admissible primitive charge is the conjugate pair q_S = +/-3 in the
q = 6Y normalization and the compatible invariant Yukawa channel set has
exactly three elements.  A countermodel battery then shows that the
source grammar does not determine scalar existence or count: the empty
configuration, the duplicate two-doublet configuration, and the
one-inert-doublet configuration each pass every grammar-visible check
(charge compatibility per coupled copy, channel arithmetic, chiral
anomaly freedom, Witten parity, and the absence of any clause that
counts scalar copies).  The verdict records scalar existence as not
source determined and scalar multiplicity as independence limited; the
one-doublet branch stays a declared completion under the registry row
scalar_existence_and_multiplicity (class conditional_open_interface,
owned by issue #616).

Family multiplicity window (#617).  On the separately declared one-Higgs
branch, Proposition 6.9 of the paper gives the conditional generation
window 3 <= N_g <= 5 from two declared clauses of Definition 6.1a: the
intrinsic CKM CP-capability clause (v) and the weak-sector ultraviolet
completability clause (vi).  The certificate verifies both edges with
exact integer and rational arithmetic in the paper's conventions: the
physical CP phase count (N_g-1)(N_g-2)/2 vanishes for N_g <= 2 and is
one at N_g = 3, and the one-loop coefficient b_SU2 = 22/3
- (1/3) N_g (N_c+1) - 1/6 with N_c = 3 and one complex Higgs doublet
gives 1/2 at N_g = 5 and -5/6 at N_g = 6.  Inside the window, N_g = 3,
4, 5 all satisfy both clauses; the A_5 character table over Q(sqrt5) is
verified orthonormal and complete with irreducible dimensions 1, 3, 3,
4, 5 (no two-dimensional irreducible exists), so a three-slot carrier
exists and k independent family copies stay grammar-admissible at every
in-window count.  The window is exact, the count inside it is not
source-selected, and N_g = 3 stays a declared completion under the
registry row family_attachment_and_multiplicity.

Fail-closed controls: a claimed CP phase at N_g = 2, a claimed
asymptotically free N_g = 6 weak sector, and a claimed two-dimensional
A_5 irreducible are each rejected by the exact arithmetic with typed
error codes.  Every decision is exact in Q or Q(sqrt5); no floating
point appears in a proof step.
"""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import echosahedral_selector_certificate as e565  # noqa: E402
import super_tannakian_matter_lift_certificate as m314  # noqa: E402

SCHEMA = "oph.multiplicity_window_certificate.v1"
MANIFEST_PATH = MODULE_DIR / "manifests" / "multiplicity_window_reference.json"
MATTER_MANIFEST_PATH = MODULE_DIR / "manifests" / "super_tannakian_matter_reference.json"

CertificateError = e565.CertificateError
require = e565.require
sha256_json = e565.sha256_json
load_json = e565.load_json
write_json = e565.write_json

F5 = m314.F5


def frac(value: Fraction) -> str:
    """Deterministic string form of an exact rational."""

    return str(Fraction(value))


# ---------------------------------------------------------------------------
# Part A (#616): the scalar compatibility block on the pinned charges
# ---------------------------------------------------------------------------


def scalar_pair_scan(y_color: Fraction, y_weak: Fraction) -> dict[int, list[list[str]]]:
    """Channel sets per scalar charge, mirroring the matter-lift enumeration.

    For every unordered pair of realized matter fields satisfying the
    triality and even-doublet selection rules, charge conservation solves
    the scalar charge exactly in the q = 6Y normalization.  The map sends
    each admissible charge to its sorted channel list.
    """

    package = m314.GENERAL_PACKAGE
    q = {field: 6 * m314._charge_at(field, y_color, y_weak) for field in package}
    fields = sorted(package)
    scan: dict[int, list[list[str]]] = {}
    for i, left in enumerate(fields):
        for right in fields[i:]:
            triality = (package[left]["triality"] + package[right]["triality"]) % 3
            duality = (package[left]["duality"] + 1 + package[right]["duality"]) % 2
            if triality != 0 or duality != 0:
                continue
            pair_charge = q[left] + q[right]
            require(
                pair_charge.denominator == 1,
                "SCALAR_SCAN",
                "a compatible matter pair has non-integral q = 6Y charge",
            )
            for scalar, q_s in (("S", -pair_charge), ("Sbar", pair_charge)):
                if q_s == 0:
                    continue
                scan.setdefault(int(q_s), []).append([left, scalar, right])
    return {charge: sorted(rows) for charge, rows in scan.items()}


def pinned_scalar_block() -> dict[str, Any]:
    """Rerun the scalar compatibility block on the pinned #314 manifest values."""

    matter_manifest = load_json(MATTER_MANIFEST_PATH)
    exterior = matter_manifest.get("exterior_matter_contract")
    require(isinstance(exterior, Mapping), "PINNED_MANIFEST", "exterior_matter_contract is missing")
    charges = exterior.get("block_trace_charges", {})
    y_color = Fraction(str(charges.get("color_block")))
    y_weak = Fraction(str(charges.get("weak_block")))
    require(
        y_color == Fraction(-1, 3) and y_weak == Fraction(1, 2),
        "PINNED_MANIFEST",
        f"pinned block charges drifted: ({y_color}, {y_weak})",
    )
    require(
        exterior.get("one_scalar") == "weak_block" and exterior.get("extra_scalars") == [],
        "PINNED_MANIFEST",
        "the pinned manifest must declare exactly the one weak-block scalar",
    )
    declared_channels = [tuple(str(x) for x in row) for row in exterior.get("yukawa_channels", [])]

    block = m314.scalar_and_channel_selection_certificate(y_color, y_weak, declared_channels)
    require(
        block["admissible_scalar_charges"] == [3, -3],
        "SCALAR_BLOCK",
        "the admissible scalar charges are not the conjugate pair (3, -3)",
    )
    require(
        block["representative_scalar_charge"] == 3,
        "SCALAR_BLOCK",
        "the declared weak-block representative charge is not +3",
    )
    channels = block["derived_channels_for_declared_representative"]
    require(
        len(channels) == 3,
        "SCALAR_BLOCK",
        f"the compatible channel set has {len(channels)} elements, expected three",
    )

    scan = scalar_pair_scan(y_color, y_weak)
    require(
        sorted(scan) == [-3, 3] and scan[3] == sorted(channels),
        "SCALAR_BLOCK",
        "the local pair scan does not reproduce the matter-lift channel set",
    )

    return {
        "pinned_matter_lift_manifest": {
            "path": "manifests/super_tannakian_matter_reference.json",
            "sha256": sha256_json(matter_manifest),
        },
        "block_trace_charges": {"color_block": frac(y_color), "weak_block": frac(y_weak)},
        "charge_normalization": "q = 6Y integer normalization",
        "admissible_scalar_charges": [3, -3],
        "conjugation_relation": block["conjugation_relation"],
        "representative_scalar_charge": 3,
        "invariant_yukawa_channels": sorted(channels),
        "channel_count": 3,
        "conditional_on": (
            "one declared color-singlet weak-doublet scalar; existence, "
            "multiplicity, and economy are outside the scan"
        ),
        "matches_pinned_manifest_channels": True,
    }


GRAMMAR_CLAUSES = (
    "charge_sum_zero_in_q6Y_per_coupled_channel",
    "triality_zero_mod_three_per_coupled_channel",
    "even_total_doublet_count_per_coupled_channel",
    "chiral_anomaly_freedom",
    "witten_parity_even",
)


def fermion_anomaly_forms(y_color: Fraction, y_weak: Fraction) -> dict[str, Any]:
    """Exact chiral anomaly forms of the realized one-family fermion package."""

    package = m314.GENERAL_PACKAGE
    half = Fraction(1, 2)
    grav = Fraction(0)
    su3 = Fraction(0)
    su2 = Fraction(0)
    u1_cubed = Fraction(0)
    weyl_doublets = 0
    for field, row in package.items():
        q_f = m314._charge_at(field, y_color, y_weak)
        grav += row["dim"] * q_f
        su3 += row["su3_fund"] * half * q_f
        su2 += row["su2_doublets"] * half * q_f
        u1_cubed += row["dim"] * q_f**3
        weyl_doublets += row["su2_doublets"]
    return {
        "grav": grav,
        "su3_sq_u1": su3,
        "su2_sq_u1": su2,
        "u1_cubed": u1_cubed,
        "weyl_doublets": weyl_doublets,
    }


def scalar_configuration_record(
    name: str,
    copies: Sequence[Mapping[str, Any]],
    scan: Mapping[int, list[list[str]]],
    y_color: Fraction,
    y_weak: Fraction,
    reading: str,
) -> dict[str, Any]:
    """Evaluate every grammar-visible check on one scalar configuration.

    Scalar copies are bosonic and non-chiral, so each copy enters every
    chiral anomaly form and the Witten doublet count with weight zero.
    The check battery records the exact value of each clause together
    with its pass flag; a countermodel demonstrates non-determination
    exactly when every clause passes.
    """

    forms = fermion_anomaly_forms(y_color, y_weak)
    checks: list[dict[str, Any]] = []

    copy_rows: list[dict[str, Any]] = []
    for copy in copies:
        charge = int(copy["charge_q6Y"])
        coupled = copy["yukawa_coupling"] == "unit"
        channels = scan.get(charge, []) if coupled else []
        require(
            charge in (3, -3),
            "SCALAR_CONFIGURATION",
            f"configuration {name} carries an incompatible scalar charge {charge}",
        )
        copy_rows.append(
            {
                "label": str(copy["label"]),
                "color": "1",
                "weak": 2,
                "charge_q6Y": charge,
                "yukawa_coupling": str(copy["yukawa_coupling"]),
                "invariant_channel_count": len(channels),
            }
        )
    checks.append(
        {
            "check": "charge_compatibility_per_copy",
            "value": [row["charge_q6Y"] for row in copy_rows],
            "pass": all(row["charge_q6Y"] in (3, -3) for row in copy_rows),
        }
    )

    coupled_channel_counts = [row["invariant_channel_count"] for row in copy_rows]
    total_channels = sum(coupled_channel_counts)
    checks.append(
        {
            "check": "channel_structure",
            "value": {
                "per_copy": coupled_channel_counts,
                "total": total_channels,
                "channels_scale_with_coupled_copies": True,
            },
            "pass": all(
                count in (0, 3) for count in coupled_channel_counts
            ),
        }
    )

    for clause in (
        "charge_sum_zero_in_q6Y_per_coupled_channel",
        "triality_zero_mod_three_per_coupled_channel",
        "even_total_doublet_count_per_coupled_channel",
    ):
        # Every coupled channel comes from the pair scan, which enforces the
        # clause by construction; an empty channel list satisfies it vacuously.
        checks.append(
            {
                "check": clause,
                "value": {"coupled_channels": total_channels},
                "pass": True,
            }
        )

    scalar_chiral_weight = 0
    anomalies = {
        key: frac(forms[key] + scalar_chiral_weight * len(copies))
        for key in ("grav", "su3_sq_u1", "su2_sq_u1", "u1_cubed")
    }
    checks.append(
        {
            "check": "chiral_anomaly_freedom",
            "value": {
                "fermion_sector": anomalies,
                "scalar_chiral_weight_per_copy": scalar_chiral_weight,
            },
            "pass": all(forms[key] == 0 for key in ("grav", "su3_sq_u1", "su2_sq_u1", "u1_cubed")),
        }
    )

    witten = forms["weyl_doublets"] % 2
    checks.append(
        {
            "check": "witten_parity_even",
            "value": {
                "weyl_fermion_doublets": forms["weyl_doublets"],
                "scalar_doublet_contribution": 0,
                "parity": witten,
            },
            "pass": witten == 0,
        }
    )

    checks.append(
        {
            "check": "no_clause_counts_scalar_copies",
            "value": {
                "grammar_clauses": list(GRAMMAR_CLAUSES),
                "clauses_referencing_scalar_count": [],
                "scalar_copies": len(copies),
            },
            "pass": True,
        }
    )

    require(
        all(row["pass"] for row in checks),
        "SCALAR_CONFIGURATION",
        f"configuration {name} fails a grammar-visible check",
    )
    return {
        "name": name,
        "scalar_copies": copy_rows,
        "copy_count": len(copies),
        "checks": checks,
        "all_checks_pass": True,
        "reading": reading,
    }


def scalar_countermodel_battery(
    y_color: Fraction, y_weak: Fraction
) -> dict[str, Any]:
    """The three non-determination countermodels of issue #616."""

    scan = scalar_pair_scan(y_color, y_weak)
    configurations = [
        scalar_configuration_record(
            "n0_no_scalar",
            [],
            scan,
            y_color,
            y_weak,
            (
                "no grammar clause fails at zero scalars: scalar existence is "
                "not forced by the source grammar, and the absent Yukawa "
                "channels violate no A1/A2/A3-visible constraint"
            ),
        ),
        scalar_configuration_record(
            "n2_duplicate_identical_charge",
            [
                {"label": "S1", "charge_q6Y": 3, "yukawa_coupling": "unit"},
                {"label": "S2", "charge_q6Y": 3, "yukawa_coupling": "unit"},
            ],
            scan,
            y_color,
            y_weak,
            (
                "two identical charge-3 doublets pass charge compatibility and "
                "channel structure per copy; the channel set scales and the "
                "non-chiral scalars contribute nothing to any anomaly form"
            ),
        ),
        scalar_configuration_record(
            "n2_one_inert",
            [
                {"label": "S1", "charge_q6Y": 3, "yukawa_coupling": "unit"},
                {"label": "S_inert", "charge_q6Y": 3, "yukawa_coupling": "zero"},
            ],
            scan,
            y_color,
            y_weak,
            (
                "one coupled and one inert doublet pass every visible check; "
                "no clause forces a nonzero Yukawa coupling on any copy"
            ),
        ),
    ]
    return {
        "grammar_visible_checks": list(GRAMMAR_CLAUSES)
        + ["charge_compatibility_per_copy", "channel_structure", "no_clause_counts_scalar_copies"],
        "configurations": configurations,
        "every_configuration_passes_every_check": True,
        "conclusion": (
            "three mutually incompatible scalar contents pass every "
            "grammar-visible check, so the source grammar determines neither "
            "scalar existence nor scalar multiplicity"
        ),
    }


def scalar_verdict() -> dict[str, Any]:
    return {
        "scalar_existence": "not_source_determined",
        "scalar_multiplicity": "independence_limited",
        "countermodels": [
            "n0_no_scalar",
            "n2_duplicate_identical_charge",
            "n2_one_inert",
        ],
        "declared_completion": (
            "one color-singlet weak-doublet scalar with q_S = +/-3 stays a "
            "declared completion of the matter branch"
        ),
        "registry_row": {
            "id": "scalar_existence_and_multiplicity",
            "class": "conditional_open_interface",
            "area": "matter",
            "owner_issue": 616,
        },
    }


# ---------------------------------------------------------------------------
# Part B (#617): the family multiplicity window
# ---------------------------------------------------------------------------


def ckm_parameter_count(n_g: int) -> dict[str, int]:
    """Exact CKM parameter count for n_g generations.

    An n x n unitary mixing matrix carries n^2 real parameters; 2n - 1
    independent rephasings remove all but (n-1)^2, of which n(n-1)/2 are
    angles and (n-1)(n-2)/2 are physical CP phases.  The identity
    (n-1)^2 - n(n-1)/2 = (n-1)(n-2)/2 is checked per value.
    """

    angles = Fraction(n_g * (n_g - 1), 2)
    phases = Fraction((n_g - 1) * (n_g - 2), 2)
    physical = Fraction((n_g - 1) ** 2)
    require(
        angles.denominator == 1 and phases.denominator == 1,
        "CP_COUNT",
        f"non-integral CKM parameter count at N_g = {n_g}",
    )
    require(
        Fraction(n_g * n_g) - Fraction(2 * n_g - 1) == physical,
        "CP_COUNT",
        f"rephasing count fails at N_g = {n_g}",
    )
    require(
        physical - angles == phases,
        "CP_COUNT",
        f"angle/phase split fails at N_g = {n_g}",
    )
    return {"angles": int(angles), "phases": int(phases)}


def cp_capability_certificate() -> dict[str, Any]:
    """The CP-capability lower edge in the paper's exact convention."""

    table = {n: ckm_parameter_count(n) for n in range(2, 7)}
    require(
        [table[n]["angles"] for n in range(2, 7)] == [1, 3, 6, 10, 15],
        "CP_COUNT",
        "the CKM angle count table drifted",
    )
    require(
        [table[n]["phases"] for n in range(2, 7)] == [0, 1, 3, 6, 10],
        "CP_COUNT",
        "the CKM phase count table drifted",
    )
    require(
        table[2]["phases"] == 0 and table[3]["phases"] == 1,
        "CP_COUNT",
        "the lower-edge phase counts drifted",
    )
    excluded = [n for n in range(1, 7) if Fraction((n - 1) * (n - 2), 2) == 0]
    require(excluded == [1, 2], "CP_COUNT", "the CP clause exclusion set drifted")
    return {
        "formula": "#(CP phases) = (N_g - 1)(N_g - 2) / 2",
        "count_table": {
            str(n): {"angles": table[n]["angles"], "phases": table[n]["phases"]}
            for n in range(2, 7)
        },
        "excluded_counts": [1, 2],
        "lower_edge": 3,
        "paper_convention": {
            "source": "paper/tex_fragments/PAPER.tex, Proposition 6.9, Step 1",
            "clause": (
                "intrinsic CKM CP capability, clause (v) of Definition 6.1a "
                "on the declared one-Higgs quark branch"
            ),
            "paper_formula": "#(CP phases) = (N_g - 1)(N_g - 2)/2",
            "matches_certificate_formula": True,
        },
        "conclusion": (
            "the declared CP-capability clause excludes N_g <= 2: a "
            "two-generation CKM matrix carries zero physical phases and the "
            "first phase appears at N_g = 3"
        ),
    }


def b_su2(n_g: int, n_c: int = 3) -> Fraction:
    """One-loop SU(2) coefficient in the paper's convention.

    The paper writes b_SU2 = 22/3 - (1/3) N_g (N_c + 1) - 1/6, where the
    final -1/6 is one complex Higgs doublet (T = 1/2 scalar contributing
    (1/3)(1/2)).  With N_c = 3 this equals 22/3 - (4/3) N_g - 1/6: the
    N_g term counts four weak-doublet Weyl fermions per generation at
    T = 1/2, each contributing (2/3)(1/2).
    """

    paper_form = Fraction(22, 3) - Fraction(1, 3) * n_g * (n_c + 1) - Fraction(1, 6)
    direct_form = Fraction(22, 3) - Fraction(4, 3) * n_g - Fraction(1, 6)
    require(
        n_c != 3 or paper_form == direct_form,
        "UV_CONVENTION",
        f"the paper form and the direct form disagree at N_g = {n_g}",
    )
    require(
        (paper_form > 0) == (2 * n_g * (n_c + 1) < 43),
        "UV_CONVENTION",
        "the paper's inequality N_g (N_c + 1) < 43/2 is not equivalent to b > 0",
    )
    return paper_form


def su2_uv_certificate() -> dict[str, Any]:
    """The SU(2) ultraviolet upper edge in the paper's exact convention."""

    values = {n: b_su2(n) for n in range(1, 9)}
    require(
        values[5] == Fraction(1, 2) and values[6] == Fraction(-5, 6),
        "UV_SIGN",
        f"the edge coefficients drifted: b(5) = {values[5]}, b(6) = {values[6]}",
    )
    require(
        all(values[n] - values[n + 1] == Fraction(4, 3) for n in range(1, 8)),
        "UV_SIGN",
        "b_SU2 does not decrease by 4/3 per generation",
    )
    excluded = [n for n in range(1, 9) if values[n] <= 0]
    require(excluded == [6, 7, 8], "UV_SIGN", "the UV exclusion set drifted")
    return {
        "formula": "b_SU2 = 22/3 - (1/3) N_g (N_c + 1) - 1/6, N_c = 3",
        "higgs_term": "-1/6 from one complex Higgs doublet (T = 1/2 scalar)",
        "value_table": {str(n): frac(values[n]) for n in range(1, 9)},
        "b_at_5": frac(values[5]),
        "b_at_6": frac(values[6]),
        "excluded_counts_from_6": True,
        "upper_edge": 5,
        "paper_convention": {
            "source": "paper/tex_fragments/PAPER.tex, Proposition 6.9, Step 2",
            "clause": (
                "weak-sector UV completability, clause (vi) of Definition 6.1a "
                "on the same declared one-Higgs branch"
            ),
            "paper_formula": "b_SU2 = 22/3 - (1/3) N_g (N_c + 1) - 1/6",
            "paper_inequality": "asymptotic freedom b_SU2 > 0, i.e. N_g (N_c + 1) < 43/2",
            "includes_higgs_term": True,
            "matches_certificate_formula": True,
        },
        "conclusion": (
            "the declared UV clause excludes N_g >= 6: b_SU2 = 1/2 > 0 at "
            "N_g = 5 and b_SU2 = -5/6 < 0 at N_g = 6"
        ),
    }


A5_CLASS_LABELS = ("e", "(2,2)", "3-cycles", "5-cycles-a", "5-cycles-b")
A5_CLASS_SIZES = (1, 15, 20, 12, 12)


def _a5_character_table() -> dict[str, tuple[F5, ...]]:
    phi = F5(Fraction(1, 2), Fraction(1, 2))
    psi = F5(Fraction(1, 2), Fraction(-1, 2))
    return {
        "1": (F5(1), F5(1), F5(1), F5(1), F5(1)),
        "3": (F5(3), F5(-1), F5(0), phi, psi),
        "3p": (F5(3), F5(-1), F5(0), psi, phi),
        "4": (F5(4), F5(0), F5(1), F5(-1), F5(-1)),
        "5": (F5(5), F5(1), F5(-1), F5(0), F5(0)),
    }


def a5_irreducible_dimensions() -> list[int]:
    """The complete A_5 irreducible dimension list from exact character arithmetic."""

    require(sum(A5_CLASS_SIZES) == 60, "A5_CHARACTERS", "the class sizes do not sum to 60")
    table = _a5_character_table()
    labels = sorted(table)
    for i, left in enumerate(labels):
        for right in labels[i:]:
            inner = F5(0)
            for k in range(5):
                inner = inner + F5(A5_CLASS_SIZES[k]) * table[left][k] * table[right][k]
            expected = F5(60) if left == right else F5(0)
            require(
                inner == expected,
                "A5_CHARACTERS",
                f"character rows {left}, {right} are not orthonormal: {inner!r}",
            )
    dims: list[int] = []
    for label in ("1", "3", "3p", "4", "5"):
        value = table[label][0]
        require(
            value.b == 0 and value.a.denominator == 1,
            "A5_CHARACTERS",
            f"irreducible dimension of {label} is not a rational integer",
        )
        dims.append(int(value.a))
    require(sorted(dims) == [1, 3, 3, 4, 5], "A5_CHARACTERS", f"dimension list drifted: {dims}")
    require(
        sum(d * d for d in dims) == 60,
        "A5_CHARACTERS",
        "the squared dimensions do not sum to the group order",
    )
    return dims


def a5_carrier_certificate() -> dict[str, Any]:
    """The screen multiplicity space: A_5 carriers by exact character arithmetic."""

    dims = a5_irreducible_dimensions()
    two_dim = _decompositions_of(dims, 2)
    require(
        two_dim == [(2, 0, 0, 0, 0)],
        "A5_IRREP_DIMENSION",
        f"unexpected two-dimensional decompositions: {two_dim}",
    )
    # The single solution stacks two trivial lines, so every two-dimensional
    # A_5 module is reducible and no two-dimensional irreducible exists.
    nontrivial = sorted(d for d in dims if d > 1)
    require(nontrivial[0] == 3, "A5_IRREP_DIMENSION", "the smallest nontrivial dimension drifted")
    return {
        "class_labels": list(A5_CLASS_LABELS),
        "class_sizes": list(A5_CLASS_SIZES),
        "character_table_orthonormal": True,
        "table_complete": "five irreducibles on five classes; squared dimensions sum to 60",
        "irreducible_dimensions": [1, 3, 3, 4, 5],
        "two_dimensional_modules": [
            {"multiplicities": list(two_dim[0]), "irreducible": False, "content": "1 + 1"}
        ],
        "no_two_dimensional_irreducible": True,
        "smallest_nontrivial_carrier_dimension": 3,
        "three_slot_carrier_exists": True,
        "lean_references": [
            "Lean/Screen/A5CharacterField.lean (character table, Galois stability)",
            "Lean/Screen/A5AngularMultiplets.lean (character_table_orthonormal)",
        ],
    }


def in_window_non_selection() -> dict[str, Any]:
    """Verify that every in-window count satisfies both clauses and that the
    grammar admits k independent family copies at every in-window count."""

    forms = fermion_anomaly_forms(Fraction(-1, 3), Fraction(1, 2))
    rows: list[dict[str, Any]] = []
    for n in (3, 4, 5):
        phases = ckm_parameter_count(n)["phases"]
        coefficient = b_su2(n)
        require(
            phases >= 1 and coefficient > 0,
            "WINDOW_MEMBERSHIP",
            f"N_g = {n} fails a window clause",
        )
        scaled = {
            key: n * forms[key] for key in ("grav", "su3_sq_u1", "su2_sq_u1", "u1_cubed")
        }
        require(
            all(value == 0 for value in scaled.values()),
            "WINDOW_MEMBERSHIP",
            f"the anomaly forms fail to vanish at {n} family copies",
        )
        witten = (n * forms["weyl_doublets"]) % 2
        require(witten == 0, "WINDOW_MEMBERSHIP", f"odd Witten parity at {n} family copies")
        rows.append(
            {
                "n_g": n,
                "cp_phases": phases,
                "b_su2": frac(coefficient),
                "cp_clause_pass": True,
                "uv_clause_pass": True,
                "k_copies_anomaly_forms": {key: frac(value) for key, value in scaled.items()},
                "k_copies_witten_parity": witten,
                "k_copies_grammar_admissible": True,
            }
        )
    return {
        "window": [3, 4, 5],
        "members": rows,
        "copy_counting_clause_exists": False,
        "copy_admissibility": (
            "anomaly freedom and Witten parity are per-family and scale "
            "linearly, so k independent copies of any family carrier stay "
            "grammar-admissible at every in-window count; no visible "
            "constraint counts copies"
        ),
        "window_is_exact": True,
        "count_inside_window_source_selected": False,
    }


def family_verdict() -> dict[str, Any]:
    return {
        "window": [3, 4, 5],
        "window_is_exact": True,
        "count_inside_window": "not_source_selected",
        "declared_completion": (
            "N_g = 3 stays a declared completion inside the conditional window"
        ),
        "registry_row": {
            "id": "family_attachment_and_multiplicity",
            "class": "conditional_open_interface",
            "area": "matter",
            "campaign_issue": 617,
        },
    }


# ---------------------------------------------------------------------------
# Fail-closed controls
# ---------------------------------------------------------------------------


def check_cp_phase_claim(n_g: int, claimed_phases: int) -> None:
    phases = ckm_parameter_count(n_g)["phases"]
    require(
        claimed_phases == phases,
        "CP_COUNT_REJECTION",
        f"an N_g = {n_g} CKM matrix carries exactly {phases} physical phases, "
        f"the claimed count {claimed_phases} is rejected",
    )


def check_asymptotic_freedom_claim(n_g: int) -> None:
    coefficient = b_su2(n_g)
    require(
        coefficient > 0,
        "UV_SIGN_REJECTION",
        f"b_SU2 = {coefficient} at N_g = {n_g}; the asymptotic-freedom claim is rejected",
    )


def check_a5_irreducible_dimension_claim(dimension: int) -> None:
    dims = a5_irreducible_dimensions()
    solutions = [
        multiplicities
        for multiplicities in _decompositions_of(dims, dimension)
        if sum(multiplicities) == 1
    ]
    require(
        len(solutions) > 0,
        "A5_IRREP_REJECTION",
        f"no A5 irreducible of dimension {dimension} exists; the complete "
        f"dimension list is {sorted(dims)}",
    )


def _decompositions_of(dims: Sequence[int], total: int) -> list[tuple[int, ...]]:
    solutions: list[tuple[int, ...]] = []

    def extend(prefix: tuple[int, ...], remaining: int) -> None:
        index = len(prefix)
        if index == len(dims):
            if remaining == 0:
                solutions.append(prefix)
            return
        for count in range(remaining // dims[index] + 1):
            extend(prefix + (count,), remaining - count * dims[index])

    extend((), total)
    return solutions


def _expect_rejection(action: Callable[[], None], code: str, claim: str) -> dict[str, Any]:
    try:
        action()
    except CertificateError as error:
        require(
            error.code == code,
            "CONTROL_CODE",
            f"control {claim} raised {error.code}, expected {code}",
        )
        return {"claim": claim, "rejected": True, "error_code": code, "detail": error.message}
    raise CertificateError("CONTROL_NOT_REJECTED", f"control claim was accepted: {claim}")


def fail_closed_controls() -> list[dict[str, Any]]:
    controls = [
        _expect_rejection(
            lambda: check_cp_phase_claim(2, 1),
            "CP_COUNT_REJECTION",
            "one CP phase in a two-generation CKM matrix",
        ),
        _expect_rejection(
            lambda: check_asymptotic_freedom_claim(6),
            "UV_SIGN_REJECTION",
            "an asymptotically free SU(2) sector at N_g = 6",
        ),
        _expect_rejection(
            lambda: check_a5_irreducible_dimension_claim(2),
            "A5_IRREP_REJECTION",
            "a two-dimensional A5 irreducible representation",
        ),
    ]
    # The same checkers accept the true values, so the rejections are sharp.
    check_cp_phase_claim(3, 1)
    check_asymptotic_freedom_claim(5)
    check_a5_irreducible_dimension_claim(3)
    return controls


# ---------------------------------------------------------------------------
# Payload, manifest, verification, CLI
# ---------------------------------------------------------------------------


def certificate_payload() -> dict[str, Any]:
    y_color = Fraction(-1, 3)
    y_weak = Fraction(1, 2)
    cp = cp_capability_certificate()
    uv = su2_uv_certificate()
    return {
        "schema": SCHEMA,
        "issues": [616, 617],
        "description": (
            "bounded scalar-response and family-multiplicity campaigns: the "
            "scalar compatibility block on the pinned matter-lift charges "
            "with its non-determination countermodel battery, and the exact "
            "conditional generation window 3 <= N_g <= 5 with in-window "
            "non-selection and fail-closed controls"
        ),
        "scalar_response_multiplicity": {
            "compatibility_block": pinned_scalar_block(),
            "countermodel_battery": scalar_countermodel_battery(y_color, y_weak),
            "verdict": scalar_verdict(),
        },
        "family_multiplicity_window": {
            "cp_capability_lower_edge": cp,
            "su2_ultraviolet_upper_edge": uv,
            "screen_carrier": a5_carrier_certificate(),
            "in_window_non_selection": in_window_non_selection(),
            "verdict": family_verdict(),
        },
        "fail_closed_controls": fail_closed_controls(),
        "arithmetic": "exact integer, rational, and Q(sqrt5) arithmetic; no floating point",
    }


def build_manifest() -> dict[str, Any]:
    body = certificate_payload()
    manifest = dict(body)
    manifest["manifest_sha256"] = "sha256:" + sha256_json(body)
    return manifest


def verify_manifest(path: Path) -> dict[str, Any]:
    stored = load_json(path)
    require(stored.get("schema") == SCHEMA, "SCHEMA", f"expected {SCHEMA}")
    body = {key: value for key, value in stored.items() if key != "manifest_sha256"}
    require(
        stored.get("manifest_sha256") == "sha256:" + sha256_json(body),
        "MANIFEST_HASH",
        "the manifest self-hash does not recompute",
    )
    recomputed = certificate_payload()
    require(
        body == recomputed,
        "MANIFEST_DRIFT",
        "the stored manifest does not equal the recomputed certificate payload",
    )
    return stored


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--write",
        action="store_true",
        help="write the reference manifest from the recomputed payload",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=MANIFEST_PATH,
        help="manifest path (default: manifests/multiplicity_window_reference.json)",
    )
    args = parser.parse_args(argv)

    if args.write or not args.manifest.exists():
        manifest = build_manifest()
        write_json(args.manifest, manifest)
        print(f"wrote {args.manifest}")
    else:
        manifest = verify_manifest(args.manifest)
        print(f"verified {args.manifest}")

    scalar = manifest["scalar_response_multiplicity"]["verdict"]
    family = manifest["family_multiplicity_window"]["verdict"]
    print(
        "scalar_existence="
        + scalar["scalar_existence"]
        + " scalar_multiplicity="
        + scalar["scalar_multiplicity"]
    )
    print(
        "window="
        + str(family["window"])
        + " count_inside_window="
        + family["count_inside_window"]
    )
    print("manifest_sha256=" + manifest["manifest_sha256"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
