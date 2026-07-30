#!/usr/bin/env python3
"""Spectral-ledger certificate for GitHub issues #612 and #609: exterior
matter-menu completeness with an exact beyond-algebra boundary.

The declared matter carrier of the pinned #314 matter-lift packet is the
exterior algebra over the five-dimensional module V = C (+) W with the
derived primitive block charges (q_C, q_W) = (-2, 3) in q = 6Y units.
Inside that declared response algebra this certificate establishes:

* the complete isotypic menu: the 32-dimensional Fock space decomposes into
  twelve summands (the vacuum line, the ten nontrivial current-module
  summands, and the top line) with every quantum datum the interface
  grammar sees (dimension, exact charge block, fermionic parity, triality,
  duality, Gauss/vacuum class); the twelve diagonal projectors are
  orthogonal idempotents whose sum is the identity on the 32-dimensional
  space, so menu completeness inside the algebra is exact, and the pair
  (charge, parity) separates the twelve summands, so isotypic factorization
  through the menu is unique;

* the exhaustive 1024-subset cross-reference: every subset of the ten
  nontrivial summands carries an exclusion class (empty, vectorlike and
  non-chiral, anomalous with its exact failing traces, or survivor), the
  survivor count is exactly two, and the survivor masks equal the
  kernel-checked masks of Lean/Screen/ExteriorSelection.lean under the
  identical component indexing;

* explicit off-menu controls: an anomaly-free vectorlike pair excluded by
  the chirality clause alone; higher representations absent from the
  exterior algebra and therefore beyond the declared menu; the sterile
  rank-one trivial-charge countermodel, which every current observable of
  the declared grammar annihilates exactly, so source data cannot exclude
  it; scalar duplicates and inert doublets recorded as owned by issue
  #616; general source-invisible direct sums covered by additivity over
  the sterile witness;

* the declared light/heavy threshold of #609: the retained light sector is
  the scan-selected chiral anomaly-free conjugate pair, every other
  in-algebra configuration carries its exact exclusion clause, the
  Gauss/vacuum lines are the declared heavy/topological side, and a
  physical decoupling or pole identification justification stays a
  separate open physical interface.

Fail-closed controls: a claimed third scan survivor, a tampered charge
assignment, the sterile summand claimed as source-visible, and a dropped
menu summand are each rejected with a typed error code.  Every decision is
exact integer or rational arithmetic; no floating point appears.
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable, Mapping

MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import echosahedral_selector_certificate as e565  # noqa: E402
import super_tannakian_matter_lift_certificate as m314  # noqa: E402

SCHEMA = "oph.matter_menu_spectral_ledger_certificate.v1"
ISSUES = [612, 609]

CertificateError = e565.CertificateError
require = e565.require
sha256_json = e565.sha256_json
load_json = e565.load_json
write_json = e565.write_json
frac_text = m314.frac_text

PINNED_MANIFEST_RELPATH = "manifests/super_tannakian_matter_reference.json"
PINNED_RECEIPT_RELPATH = "receipts/super_tannakian_matter_reference.receipt.json"
LEAN_RELPATH = "Lean/Screen/ExteriorSelection.lean"
LEAN_EVEN_MASK = 796
LEAN_ODD_MASK = 227
ANOMALY_NAMES = ("grav", "su3", "su2", "u1_cubed")

# Component index -> (color exterior degree, weak exterior degree).  The
# component order is exactly m314._scan_components, which equals the Fin 10
# order of Lean/Screen/ExteriorSelection.lean.
COMPONENT_DEGREES: dict[int, tuple[int, int]] = {
    0: (1, 0),
    1: (0, 1),
    2: (2, 0),
    3: (1, 1),
    4: (0, 2),
    5: (1, 2),
    6: (2, 1),
    7: (3, 0),
    8: (2, 2),
    9: (3, 1),
}
DEGREES_TO_COMPONENT = {v: k for k, v in COMPONENT_DEGREES.items()}

# Field correspondence with the pinned #314 rank-fifteen block ledger: the
# even-parity components carry the field labels of the matter-lift receipt
# and the odd-parity components carry their charge conjugates.
FIELD_LABELS: dict[int, str] = {
    2: "u_c",
    3: "Q",
    4: "e_c",
    8: "d_c",
    9: "L",
    0: "d_c_conjugate",
    1: "L_conjugate",
    5: "u_c_conjugate",
    6: "Q_conjugate",
    7: "e_c_conjugate",
}
PINNED_FIELD_DIMENSIONS = {"Q": 6, "u_c": 3, "e_c": 1, "d_c": 3, "L": 2}

CLASS_CODES = {
    "empty": "E",
    "vectorlike_nonchiral": "V",
    "anomalous": "A",
    "survivor": "S",
}


def color_dimension(color: str) -> int:
    return 3 if color in ("3", "3bar") else 1


# ---------------------------------------------------------------------------
# Pinned upstream: the #314 matter-lift manifest and receipt
# ---------------------------------------------------------------------------


def load_pinned_matter_lift(base_dir: Path | None = None) -> dict[str, Any]:
    """Load, hash, and cross-check the pinned matter-lift packet.

    The menu charges are taken from the pinned manifest, and the receipt is
    checked for the rank-fifteen selection, the two-survivor scan with the
    derived integer block charges, the five multiplicity-free irreducible
    blocks with their dimensions, and the six emitted kernel residues.
    """

    base = base_dir or MODULE_DIR
    manifest_path = base / "manifests" / "super_tannakian_matter_reference.json"
    receipt_path = base / "receipts" / "super_tannakian_matter_reference.receipt.json"
    manifest = load_json(manifest_path)
    receipt = load_json(receipt_path)
    require(
        manifest.get("schema") == m314.SCHEMA,
        "PINNED_SCHEMA",
        f"the pinned matter-lift manifest does not carry schema {m314.SCHEMA}",
    )
    require(
        receipt.get("schema") == m314.RECEIPT_SCHEMA and receipt.get("issue") == 314,
        "PINNED_SCHEMA",
        "the pinned receipt is not a #314 matter-lift receipt",
    )
    manifest_sha256 = sha256_json(manifest)
    receipt_sha256 = sha256_json(receipt)
    require(
        receipt.get("manifest_sha256") == manifest_sha256,
        "PINNED_HASH_LINK",
        "the pinned receipt does not certify the pinned manifest",
    )
    conditional_gate = receipt.get("conditional_algebraic_gate", {})
    physical_gate = receipt.get("physical_source_gate", {})
    require(
        conditional_gate.get("passed") is True
        and physical_gate.get("passed") is False
        and physical_gate.get("matter_lift_source_bound") is False
        and physical_gate.get("upstream_current_representation_source_bound")
        is False,
        "PINNED_SCOPE",
        "the pinned #314 packet must remain a conditional matter fixture with "
        "its physical source gate closed",
    )

    charges = manifest["exterior_matter_contract"]["block_trace_charges"]
    y_color = m314.parse_rational(charges.get("color_block"), "PINNED_CHARGES")
    y_weak = m314.parse_rational(charges.get("weak_block"), "PINNED_CHARGES")
    require(
        3 * y_color + 2 * y_weak == 0,
        "PINNED_CHARGES",
        "the pinned charge pair is not trace balanced",
    )
    a = 6 * y_color
    b = 6 * y_weak
    require(
        a.denominator == 1 and b.denominator == 1 and (int(a), int(b)) in ((-2, 3), (2, -3)),
        "PINNED_CHARGES",
        "the pinned charge pair is not a primitive integral balanced pair",
    )
    a_int, b_int = int(a), int(b)

    scan = receipt.get("matter_selection_scan", {})
    require(
        scan.get("subsets_enumerated") == 1024
        and scan.get("survivor_count") == 2
        and scan.get("derived_block_charges")
        == {"a": a_int, "b": b_int, "normalization": 6},
        "PINNED_SCAN",
        "the pinned receipt scan does not match the derived block charges",
    )
    require(
        receipt.get("selection", {}).get("projector_rank") == 15
        and receipt.get("realized_package", {}).get("dimension") == 15,
        "PINNED_RANK",
        "the pinned receipt does not carry the rank-fifteen selection",
    )
    require(
        receipt.get("realized_package", {}).get("irreducible_block_commutants") == 5,
        "PINNED_COMMUTANTS",
        "the pinned receipt does not carry five scalar block commutants",
    )
    field_dimensions = {
        name: row.get("dimension")
        for name, row in receipt.get("realized_package", {}).get("fields", {}).items()
    }
    require(
        field_dimensions == PINNED_FIELD_DIMENSIONS,
        "PINNED_FIELDS",
        f"the pinned field dimensions drifted: {field_dimensions}",
    )
    invariant_sector = receipt.get("selection", {}).get("derived_invariant_sector", "")
    require(
        "vacuum" in invariant_sector and "top line" in invariant_sector,
        "PINNED_INVARIANT_SECTOR",
        "the pinned invariant sector is not the vacuum and top lines",
    )
    residues = receipt.get("kernel_emission", {}).get(
        "kernel_residues_modulo_deck_translations", []
    )
    require(
        len(residues) == 6,
        "PINNED_KERNEL",
        f"the pinned kernel emission does not carry six residues, got {len(residues)}",
    )
    kernel_residues = [
        (
            Fraction(row["u1_phase_turns"]),
            int(row["su3_center_power"]),
            int(row["su2_center_power"]),
        )
        for row in residues
    ]
    return {
        "manifest_path": PINNED_MANIFEST_RELPATH,
        "manifest_sha256": manifest_sha256,
        "receipt_path": PINNED_RECEIPT_RELPATH,
        "receipt_sha256": receipt_sha256,
        "conditional_algebraic_gate_passed": True,
        "physical_source_gate_passed": False,
        "y_color": y_color,
        "y_weak": y_weak,
        "a": a_int,
        "b": b_int,
        "kernel_residues": kernel_residues,
    }


# ---------------------------------------------------------------------------
# Part 1: the complete isotypic menu with exact projector completeness
# ---------------------------------------------------------------------------


def menu_ledger(
    pinned: Mapping[str, Any],
    *,
    drop_summand_index: int | None = None,
) -> dict[str, Any]:
    """Enumerate the complete summand menu of Lambda^* V and certify it.

    Each of the 32 Fock basis states is assigned to exactly one of twelve
    summands by its (color degree, weak degree) bidegree.  The diagonal
    summand projectors are orthogonal idempotents; their sum is verified to
    equal the identity on the 32-dimensional space.  Dropping a summand
    (control lane) breaks that identity and fails closed.
    """

    a = pinned["a"]
    b = pinned["b"]
    components = m314._scan_components(a, b)

    # Consistency of the component table with the bidegree grammar.
    for i, component in enumerate(components):
        k_c, k_w = COMPONENT_DEGREES[i]
        require(
            component["q"] == k_c * a + k_w * b,
            "MENU_CHARGE_MISMATCH",
            f"component {i} charge {component['q']} differs from the bidegree charge {k_c * a + k_w * b}",
        )
        require(
            component["parity"] == (k_c + k_w) % 2,
            "MENU_PARITY_MISMATCH",
            f"component {i} parity does not equal its exterior degree parity",
        )
        partner = next(
            j
            for j, other in enumerate(components)
            if other["pair"] == component["pair"] and other["side"] == 1 - component["side"]
        )
        pk_c, pk_w = COMPONENT_DEGREES[partner]
        require(
            (pk_c, pk_w) == (3 - k_c, 2 - k_w) and components[partner]["q"] == -component["q"],
            "MENU_CONJUGATION",
            f"component {i} and its conjugate {partner} are not complementary bidegrees with opposite charge",
        )

    fock = m314.FockSpace()
    summand_keys: list[Any] = ["vacuum"] + list(range(10)) + ["top"]
    assignment: dict[Any, list[int]] = {key: [] for key in summand_keys}
    for n, subset in enumerate(fock.subsets):
        k_c = sum(1 for i in subset if i < 3)
        k_w = sum(1 for i in subset if i >= 3)
        if (k_c, k_w) == (0, 0):
            key: Any = "vacuum"
        elif (k_c, k_w) == (3, 2):
            key = "top"
        else:
            key = DEGREES_TO_COMPONENT[(k_c, k_w)]
        assignment[key].append(n)
        expected_charge = 0 if key in ("vacuum", "top") else components[key]["q"]
        require(
            k_c * a + k_w * b == expected_charge,
            "MENU_CHARGE_MISMATCH",
            f"Fock basis state {n} carries charge {k_c * a + k_w * b}, summand expects {expected_charge}",
        )
    require(
        sum(len(indices) for indices in assignment.values()) == 32
        and len(assignment) == 12,
        "MENU_PARTITION",
        "the twelve summands do not partition the 32 Fock basis states",
    )

    rows: list[dict[str, Any]] = []
    for key in summand_keys:
        indices = assignment[key]
        if key == "vacuum":
            label, color, weak, q, degrees = "vacuum_line", "1", 1, 0, (0, 0)
        elif key == "top":
            label, color, weak, q, degrees = "top_line", "1", 1, 0, (3, 2)
        else:
            component = components[key]
            label = f"component_{key}"
            color = component["color"]
            weak = component["weak"]
            q = component["q"]
            degrees = COMPONENT_DEGREES[key]
        dimension = color_dimension(color) * weak
        require(
            len(indices) == dimension,
            "MENU_DIMENSION_MISMATCH",
            f"summand {label} collects {len(indices)} basis states, expected {dimension}",
        )
        rows.append(
            {
                "summand": label,
                "component_index": key if isinstance(key, int) else None,
                "field_label": FIELD_LABELS.get(key) if isinstance(key, int) else None,
                "exterior_degree": {
                    "color": degrees[0],
                    "weak": degrees[1],
                    "total": degrees[0] + degrees[1],
                },
                "color_factor": color,
                "weak_dimension": weak,
                "dimension": dimension,
                "charge_q6": q,
                "charge_y": frac_text(Fraction(q, 6)),
                "fermionic_parity": "odd" if (degrees[0] + degrees[1]) % 2 == 1 else "even",
                "triality": degrees[0] % 3,
                "duality": degrees[1] % 2,
                "gauss_vacuum_class": key in ("vacuum", "top"),
                "fock_basis_indices": indices,
            }
        )

    # Exact projector completeness on the 32-dimensional space.
    identity = m314.cidentity(fock.dim)
    projectors: dict[Any, list[list[Any]]] = {}
    for key in summand_keys:
        projector = m314.czeros(fock.dim)
        for n in assignment[key]:
            projector[n][n] = m314.IONE
        require(
            m314.c_is_zero(m314.csub(m314.cmul(projector, projector), projector))
            and m314.c_is_zero(m314.csub(m314.cdagger(projector), projector)),
            "MENU_PROJECTOR",
            f"the projector of summand {key} is not an orthogonal idempotent",
        )
        projectors[key] = projector
    for pos, key in enumerate(summand_keys):
        for other in summand_keys[pos + 1 :]:
            overlap = set(assignment[key]) & set(assignment[other])
            require(
                not overlap
                and all(
                    (projectors[key][n][n] * projectors[other][n][n]).is_zero()
                    for n in range(fock.dim)
                ),
                "MENU_PROJECTOR_ORTHOGONALITY",
                f"summand projectors {key} and {other} are not orthogonal",
            )
    total = m314.czeros(fock.dim)
    for key in summand_keys:
        if drop_summand_index is not None and key == drop_summand_index:
            continue
        total = m314.cadd(total, projectors[key])
    require(
        m314.c_is_zero(m314.csub(total, identity)),
        "PROJECTOR_COMPLETENESS",
        "the summand projectors do not sum to the identity on the 32-dimensional space",
    )

    # Grammar separation: charge alone separates the ten nontrivial
    # summands, and parity separates the two trivial-charge lines.
    nontrivial_charges = [component["q"] for component in components]
    require(
        len(set(nontrivial_charges)) == 10 and 0 not in nontrivial_charges,
        "MENU_SEPARATION",
        "the ten nontrivial summand charges are not distinct and nonzero",
    )
    signatures = {(row["charge_q6"], row["fermionic_parity"]) for row in rows}
    require(
        len(signatures) == 12,
        "MENU_SEPARATION",
        "the pair (charge, parity) does not separate the twelve summands",
    )

    return {
        "summands": rows,
        "summand_count": 12,
        "nontrivial_summand_count": 10,
        "total_dimension": 32,
        "projector_completeness": {
            "orthogonal_idempotents": 12,
            "pairwise_orthogonality_checks": 66,
            "sum_equals_identity_on_32": True,
            "statement": (
                "the twelve diagonal summand projectors are orthogonal "
                "idempotents whose sum is exactly the identity on Lambda^* V"
            ),
        },
        "grammar_separation": {
            "nontrivial_charges_distinct_and_nonzero": True,
            "charge_and_parity_separate_all_twelve": True,
            "statement": (
                "every nontrivial summand carries a distinct nonzero charge; "
                "the two trivial-charge lines are separated by fermionic "
                "parity, so the grammar signature determines the summand"
            ),
        },
        "isotypic_completeness": {
            "resolution_of_identity": True,
            "irreducibility_source": (
                "the pinned #314 receipt certifies the five realized blocks "
                "as irreducible with scalar commutants and the module as "
                "multiplicity free; the odd-sector components are their "
                "exact charge conjugates under the verified wedge pairing"
            ),
            "statement": (
                "every finite object of the declared exterior response "
                "algebra decomposes through the twelve summand projectors, "
                "and the decomposition is unique because the grammar "
                "signature separates the summands; factorization through "
                "the menu is exact inside the declared algebra"
            ),
            "verdict": "exact",
        },
    }


# ---------------------------------------------------------------------------
# Part 2: the exhaustive 1024-subset exclusion classification
# ---------------------------------------------------------------------------


def classify_subsets(
    components: list[dict[str, Any]],
    *,
    claimed_extra_survivor: int | None = None,
) -> dict[str, Any]:
    """Classify all 1024 subsets of the ten nontrivial summands.

    Bit i of the mask selects component i, in exactly the indexing of
    m314._scan_components and Lean/Screen/ExteriorSelection.lean.  Every
    subset carries one exclusion class: empty, vectorlike_nonchiral,
    anomalous (with the exact failing traces), or survivor.  The survivor
    set is verified against the pinned Lean parity masks, and the anomaly
    tally on both parity masks is verified to vanish exactly.
    """

    even_mask = sum(1 << i for i, c in enumerate(components) if c["parity"] == 0)
    odd_mask = sum(1 << i for i, c in enumerate(components) if c["parity"] == 1)
    require(
        even_mask == LEAN_EVEN_MASK and odd_mask == LEAN_ODD_MASK,
        "LEAN_MASK_PIN",
        f"the parity masks ({even_mask}, {odd_mask}) do not equal the pinned Lean masks (796, 227)",
    )
    conjugate_of = {
        i: next(
            j
            for j, other in enumerate(components)
            if other["pair"] == c["pair"] and other["side"] == 1 - c["side"]
        )
        for i, c in enumerate(components)
    }

    # The pinned anomaly tally: all four traces vanish on both parity masks.
    for mask in (even_mask, odd_mask):
        selection = [components[i] for i in range(10) if mask >> i & 1]
        anomalies = m314._scan_anomalies(selection)
        require(
            all(anomalies[name] == 0 for name in ANOMALY_NAMES),
            "PINNED_ANOMALY_TALLY",
            f"the anomaly tally on parity mask {mask} does not vanish: {anomalies}",
        )

    tally = {"empty": 0, "vectorlike_nonchiral": 0, "anomalous": 0, "survivor": 0}
    anomalous_breakdown: dict[str, int] = {}
    class_chars: list[str] = []
    detailed_table: list[dict[str, Any]] = []
    survivor_masks: set[int] = set()
    for mask in range(1024):
        selected = [i for i in range(10) if mask >> i & 1]
        selection = [components[i] for i in selected]
        anomalies = m314._scan_anomalies(selection)
        chiral = all(conjugate_of[i] not in selected for i in selected)
        failing = [name for name in ANOMALY_NAMES if anomalies[name] != 0]
        if mask == 0:
            exclusion_class = "empty"
        elif not chiral:
            exclusion_class = "vectorlike_nonchiral"
        elif failing:
            exclusion_class = "anomalous"
            signature = "+".join(failing)
            anomalous_breakdown[signature] = anomalous_breakdown.get(signature, 0) + 1
        else:
            exclusion_class = "survivor"
            survivor_masks.add(mask)
        tally[exclusion_class] += 1
        class_chars.append(CLASS_CODES[exclusion_class])
        detailed_table.append(
            {
                "mask": mask,
                "components": selected,
                "class": exclusion_class,
                "chiral": chiral,
                "failing_anomalies": failing,
                "anomaly_traces": {name: anomalies[name] for name in ANOMALY_NAMES},
            }
        )

    if claimed_extra_survivor is not None:
        survivor_masks = survivor_masks | {claimed_extra_survivor}
    require(
        len(survivor_masks) == 2,
        "SURVIVOR_COUNT",
        f"expected exactly two survivors, got {len(survivor_masks)}",
    )
    require(
        survivor_masks == {even_mask, odd_mask},
        "SURVIVOR_SET",
        f"the survivor masks {sorted(survivor_masks)} are not the parity masks",
    )

    survivor_rows = []
    for mask in sorted(survivor_masks):
        selected = [i for i in range(10) if mask >> i & 1]
        selection = [components[i] for i in selected]
        dimension = sum(color_dimension(c["color"]) * c["weak"] for c in selection)
        require(
            dimension == 15,
            "SURVIVOR_RANK",
            f"survivor mask {mask} has dimension {dimension}, expected the pinned rank fifteen",
        )
        doublet_slots = sum(
            color_dimension(c["color"]) for c in selection if c["weak"] == 2
        )
        require(
            doublet_slots == 4,
            "SURVIVOR_WITTEN",
            f"survivor mask {mask} carries {doublet_slots} weak doublet slots, expected four",
        )
        survivor_rows.append(
            {
                "mask": mask,
                "components": selected,
                "parity_sector": "even" if mask == even_mask else "odd",
                "dimension": dimension,
                "weak_doublet_slots": doublet_slots,
                "anomaly_traces": {name: 0 for name in ANOMALY_NAMES},
                "chiral": True,
            }
        )

    return {
        "subsets_enumerated": 1024,
        "indexing": (
            "bit i of the mask selects component i of the ten nontrivial "
            "summands, in exactly the order of the matter-lift scan table "
            "and the Fin 10 arrays of Lean/Screen/ExteriorSelection.lean"
        ),
        "exclusion_classes": {
            "empty": "the empty selection",
            "vectorlike_nonchiral": (
                "the selection contains a component together with its "
                "charge conjugate and fails the chirality clause"
            ),
            "anomalous": (
                "the selection is chiral and at least one of the four "
                "anomaly traces (grav, su3, su2, u1_cubed) is nonzero; the "
                "failing traces are recorded exactly"
            ),
            "survivor": "the selection is nonempty, chiral, and anomaly free",
        },
        "tally": tally,
        "anomalous_breakdown_by_failing_traces": dict(sorted(anomalous_breakdown.items())),
        "class_string_legend": {v: k for k, v in CLASS_CODES.items()},
        "class_string": "".join(class_chars),
        "detailed_table_sha256": sha256_json(detailed_table),
        "survivor_count": 2,
        "survivors": survivor_rows,
        "survivors_are_conjugate_pair": True,
        "lean_cross_reference": {
            "file": LEAN_RELPATH,
            "theorems": [
                "selection_unique",
                "parity_sectors_survive",
                "conj_exchanges_survivors",
                "witten_automatic",
            ],
            "even_mask": LEAN_EVEN_MASK,
            "odd_mask": LEAN_ODD_MASK,
            "agreement": True,
        },
    }


# ---------------------------------------------------------------------------
# Part 3: off-menu controls and the sterile countermodel
# ---------------------------------------------------------------------------


def sterile_countermodel(
    pinned: Mapping[str, Any],
    *,
    claim_source_visible: bool = False,
) -> dict[str, Any]:
    """The rank-one trivial-charge summand adjoined outside the algebra.

    The sterile line carries charge 0, triality 0, duality 0, color and
    weak dimension 1.  Every current observable of the declared grammar
    annihilates it exactly: the charge readback, the triality and duality
    classes, the weak Casimir indicator, every anomaly-form contribution,
    and every emitted kernel residue phase are zero.  A claim that the
    summand is source-visible fails closed against these exact zeros.
    """

    q, t, d = 0, 0, 0
    dimension = 1
    couplings: dict[str, Any] = {
        "u1_charge_readback": q,
        "triality_class": t,
        "duality_class": d,
        "weak_casimir_indicator": 0,
        "anomaly_form_contributions": {
            "grav": dimension * q,
            "su3": 0,
            "su2": 0,
            "u1_cubed": dimension * q**3,
            "witten_doublet_slots": 0,
        },
    }
    kernel_phases = []
    for r, a_c, b_c in pinned["kernel_residues"]:
        phase = r * q + Fraction(a_c * t, 3) + Fraction(b_c * d, 2)
        require(
            phase.denominator == 1,
            "STERILE_KERNEL",
            "an emitted kernel residue acts nontrivially on the sterile line",
        )
        kernel_phases.append(
            {
                "residue": {
                    "u1_phase_turns": frac_text(r),
                    "su3_center_power": a_c,
                    "su2_center_power": b_c,
                },
                "phase_integral": True,
            }
        )
    numeric_couplings = [
        couplings["u1_charge_readback"],
        couplings["triality_class"],
        couplings["duality_class"],
        couplings["weak_casimir_indicator"],
        *couplings["anomaly_form_contributions"].values(),
    ]
    all_zero = all(value == 0 for value in numeric_couplings)
    require(
        all_zero,
        "STERILE_CONSTRUCTION",
        "the sterile line construction carries a nonzero grammar coupling",
    )
    if claim_source_visible:
        require(
            any(value != 0 for value in numeric_couplings),
            "STERILE_SOURCE_INVISIBLE",
            "the sterile summand cannot be claimed source-visible: every "
            "declared readback and current-observable coupling is exactly zero",
        )
    return {
        "construction": (
            "one rank-one line with charge 0, triality 0, duality 0, "
            "adjoined as an external direct summand to the declared algebra"
        ),
        "couplings": couplings,
        "kernel_residue_phases": kernel_phases,
        "all_current_observable_couplings_zero": True,
        "menu_contrast": {
            "visible_nontrivial_summands": 10,
            "statement": (
                "every nontrivial menu summand carries a nonzero charge "
                "readback, so the grammar sees the full menu; the sterile "
                "line evades every declared readback"
            ),
        },
        "scalar_channel_arithmetic": {
            "admissible_pairing": ["L", "S", "sterile_line"],
            "charge_sum_q6": 0,
            "triality_sum_mod3": 0,
            "duality_sum_mod2": 0,
            "statement": (
                "the trilinear selection rules admit the pairing of the "
                "sterile line with the L block and the declared scalar; "
                "such a coupling lives in the scalar-sector interface owned "
                "by issues #616 and #609 and is not a current observable, "
                "so it does not make the sterile line source-visible here"
            ),
        },
        "verdict": "sterile_countermodel_source_invisible",
        "consequence": (
            "source data cannot exclude the sterile summand; the "
            "beyond-algebra completeness verdict is independence_limited "
            "with this witness"
        ),
    }


def off_menu_controls(
    pinned: Mapping[str, Any],
    menu: Mapping[str, Any],
    sterile: Mapping[str, Any],
) -> dict[str, Any]:
    """One row per off-menu control, each with its exact verdict."""

    components = m314._scan_components(pinned["a"], pinned["b"])

    # (a) A vectorlike pair inside the menu: component 0 with its conjugate 8.
    vectorlike_selection = [components[0], components[8]]
    vectorlike_anomalies = m314._scan_anomalies(vectorlike_selection)
    require(
        all(vectorlike_anomalies[name] == 0 for name in ANOMALY_NAMES),
        "OFFMENU_VECTORLIKE",
        f"the vectorlike control pair is not anomaly free: {vectorlike_anomalies}",
    )
    require(
        components[0]["pair"] == components[8]["pair"]
        and components[0]["side"] != components[8]["side"],
        "OFFMENU_VECTORLIKE",
        "the vectorlike control pair is not a conjugate pair",
    )

    # (b) Higher representations are absent from the declared algebra.
    color_factors = sorted({row["color_factor"] for row in menu["summands"]})
    weak_dimensions = sorted({row["weak_dimension"] for row in menu["summands"]})
    require(
        color_factors == ["1", "3", "3bar"] and weak_dimensions == [1, 2],
        "OFFMENU_HIGHER_REP",
        "the menu factor inventory drifted from {1, 3, 3bar} x {1, 2}",
    )
    single_factor_dimensions = sorted(
        {
            row["dimension"]
            for row in menu["summands"]
            if row["weak_dimension"] == 1 or row["color_factor"] == "1"
        }
    )
    require(
        set(single_factor_dimensions) <= {1, 2, 3},
        "OFFMENU_HIGHER_REP",
        "a single-factor summand exceeds the fundamental dimensions",
    )

    return {
        "vectorlike_pair": {
            "components": [0, 8],
            "mask": (1 << 0) | (1 << 8),
            "anomaly_traces": {name: 0 for name in ANOMALY_NAMES},
            "chiral": False,
            "verdict": "anomaly_free_but_fails_chirality_clause",
            "location": "inside the declared menu; excluded by the exact chirality clause alone",
        },
        "higher_representation": {
            "candidates": [
                "the six-dimensional symmetric square of the color triplet as a single color block",
                "a five-dimensional single-block current module (spin-2 style A5 irrep)",
            ],
            "menu_color_factors": color_factors,
            "menu_weak_dimensions": weak_dimensions,
            "single_factor_summand_dimensions": single_factor_dimensions,
            "verdict": "not_a_summand_of_the_declared_exterior_algebra",
            "location": "outside the declared menu; carried by the beyond-algebra boundary verdict",
        },
        "neutral_singlet_sterile": sterile,
        "scalar_duplicates_inert_doublets": {
            "verdict": "owned_by_issue_616",
            "analyzed_here": False,
            "note": (
                "scalar duplicates and inert doublets are scalar-response "
                "and multiplicity questions; issue #616 owns them"
            ),
        },
        "source_invisible_direct_sums": {
            "verdict": "covered_by_sterile_countermodel",
            "note": (
                "every declared readback is additive over direct summands, "
                "so a direct sum of trivial-charge trivial-isotype lines "
                "carries the same exact zero couplings as the rank-one "
                "sterile witness"
            ),
        },
    }


# ---------------------------------------------------------------------------
# Part 4: the declared light/heavy threshold
# ---------------------------------------------------------------------------


def light_heavy_threshold() -> dict[str, Any]:
    return {
        "status": "declared_threshold",
        "retained_light_sector": {
            "description": (
                "the unordered charge-conjugate pair of rank-fifteen chiral "
                "anomaly-free selections, exactly the two scan survivors"
            ),
            "masks": [LEAN_ODD_MASK, LEAN_EVEN_MASK],
            "rank": 15,
        },
        "in_algebra_exclusions": (
            "every other in-algebra configuration carries its exact "
            "exclusion clause from the 1024-subset classification: empty, "
            "vectorlike_nonchiral, or anomalous with its failing traces"
        ),
        "heavy_topological_side": {
            "lines": ["vacuum_line", "top_line"],
            "description": (
                "the Gauss-invariant trivial-charge lines of the algebra "
                "are the declared heavy/topological side of the threshold"
            ),
        },
        "physical_decoupling_interface": {
            "status": "separate_open_physical_interface",
            "note": (
                "a physical decoupling or laboratory pole identification "
                "justification for the threshold is not supplied by this "
                "certificate; the threshold is declared with its exact "
                "in-algebra clauses only"
            ),
        },
    }


# ---------------------------------------------------------------------------
# Part 6: fail-closed controls
# ---------------------------------------------------------------------------


def fail_closed_control_cases(
    pinned: Mapping[str, Any],
) -> list[tuple[str, Callable[[], Any], str, str]]:
    components = m314._scan_components(pinned["a"], pinned["b"])
    tampered = m314._scan_components(pinned["a"], pinned["b"])
    tampered[3] = dict(tampered[3], q=tampered[3]["q"] + 1)
    return [
        (
            "claimed_third_scan_survivor",
            lambda: classify_subsets(components, claimed_extra_survivor=273),
            "SURVIVOR_COUNT",
            "a claimed third survivor mask is rejected against the exact two-survivor scan",
        ),
        (
            "tampered_charge_assignment",
            lambda: classify_subsets(tampered),
            "PINNED_ANOMALY_TALLY",
            "shifting one component charge breaks the anomaly tally on the pinned parity masks",
        ),
        (
            "sterile_summand_claimed_source_visible",
            lambda: sterile_countermodel(pinned, claim_source_visible=True),
            "STERILE_SOURCE_INVISIBLE",
            "the sterile summand cannot be claimed source-visible: its couplings are exactly zero",
        ),
        (
            "dropped_menu_summand",
            lambda: menu_ledger(pinned, drop_summand_index=4),
            "PROJECTOR_COMPLETENESS",
            "dropping one summand breaks the projector resolution of the identity",
        ),
    ]


def fail_closed_control_rows(pinned: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name, runner, expected_code, statement in fail_closed_control_cases(pinned):
        actual_code = "ACCEPTED"
        try:
            runner()
        except CertificateError as exc:
            actual_code = exc.code
        require(
            actual_code == expected_code,
            "FAIL_CLOSED_CONTROL",
            f"{name}: expected {expected_code}, got {actual_code}",
        )
        rows.append(
            {
                "name": name,
                "expected_error": expected_code,
                "actual_error": actual_code,
                "statement": statement,
                "passed": True,
            }
        )
    return rows


# ---------------------------------------------------------------------------
# Certificate payload
# ---------------------------------------------------------------------------


def certificate_payload(base_dir: Path | None = None) -> dict[str, Any]:
    pinned = load_pinned_matter_lift(base_dir)
    components = m314._scan_components(pinned["a"], pinned["b"])
    menu = menu_ledger(pinned)
    classification = classify_subsets(components)
    python_scan = m314.exterior_selection_scan(pinned["y_color"], pinned["y_weak"])
    require(
        python_scan["survivor_count"] == 2
        and python_scan["survivors_equal_parity_sectors"] is True
        and python_scan["derived_block_charges"]
        == {"a": pinned["a"], "b": pinned["b"], "normalization": 6},
        "SCAN_AGREEMENT",
        "the matter-lift selection scan does not agree with the ledger classification",
    )
    sterile = sterile_countermodel(pinned)
    controls = off_menu_controls(pinned, menu, sterile)
    controls_rows = fail_closed_control_rows(pinned)

    return {
        "schema": SCHEMA,
        "issues": ISSUES,
        "description": (
            "Spectral ledger for matter-menu completeness inside the "
            "declared exterior response algebra Lambda^* V over the "
            "five-dimensional module with the derived primitive block "
            "charges. The twelve-summand isotypic menu is complete by an "
            "exact projector resolution of the identity, the 1024-subset "
            "exclusion classification is exhaustive with exactly the two "
            "kernel-checked survivors, the off-menu controls carry exact "
            "verdicts, the sterile trivial-charge countermodel bounds the "
            "beyond-algebra claim, and the light/heavy threshold is "
            "declared with its exact in-algebra clauses."
        ),
        "pinned_matter_lift": {
            "manifest_path": pinned["manifest_path"],
            "manifest_sha256": pinned["manifest_sha256"],
            "receipt_path": pinned["receipt_path"],
            "receipt_sha256": pinned["receipt_sha256"],
            "checks": {
                "conditional_fixture_only": True,
                "physical_source_gate_closed": True,
                "charge_pair_matches_derived_primitive_pair": True,
                "rank_fifteen_selection": True,
                "two_survivor_scan_with_derived_charges": True,
                "five_irreducible_blocks_multiplicity_free": True,
                "field_dimensions": PINNED_FIELD_DIMENSIONS,
                "kernel_residues": 6,
            },
        },
        "declared_algebra": {
            "carrier": (
                "the exterior algebra Lambda^* V over the five-dimensional "
                "matter carrier V = C (+) W of the pinned matter lift"
            ),
            "dimension": 32,
            "derived_block_charges": {
                "a": pinned["a"],
                "b": pinned["b"],
                "normalization": 6,
            },
            "charge_pair_y": {
                "color_block": frac_text(pinned["y_color"]),
                "weak_block": frac_text(pinned["y_weak"]),
            },
        },
        "menu": menu,
        "subset_classification": classification,
        "scan_agreement": {
            "matter_lift_scan_survivor_count": 2,
            "matter_lift_scan_matches_ledger": True,
            "lean_masks_match_ledger": True,
        },
        "off_menu_controls": controls,
        "light_heavy_threshold": light_heavy_threshold(),
        "verdicts": {
            "menu_completeness_inside_declared_algebra": "exact",
            "every_admissible_in_algebra_object_factors_through_menu": "exact",
            "beyond_declared_algebra": "independence_limited",
            "beyond_declared_algebra_witness": (
                "the sterile rank-one trivial-charge summand: every declared "
                "readback annihilates it, so source data cannot exclude it"
            ),
            "light_heavy_threshold": "declared_with_exact_in_algebra_clauses",
            "continuum_spin_statistics_and_laboratory_identification": (
                "separate_physical_interfaces"
            ),
        },
        "fail_closed_controls": controls_rows,
        "claim_boundary": {
            "proves": (
                "menu completeness inside the declared exterior response "
                "algebra by an exact projector resolution of the identity, "
                "the exhaustive exclusion classification of all 1024 "
                "in-algebra selections with exactly the two conjugate "
                "rank-fifteen survivors, exact off-menu control verdicts, "
                "and the declared light/heavy threshold with exact "
                "in-algebra clauses"
            ),
            "does_not_close": [
                "exclusion of sectors outside the declared algebra: the sterile countermodel is source-invisible, so the beyond-algebra verdict is independence_limited",
                "scalar multiplicity, scalar duplicates, and inert doublets (issue #616)",
                "response completeness of the declared algebra itself (issue #611)",
                "a physical decoupling or laboratory pole identification justification for the declared threshold",
                "continuum spin-statistics and laboratory identification of any matter observable",
            ],
        },
        "verifier_command": (
            "python3 code/a5_closure/matter_menu_spectral_ledger_certificate.py verify "
            "--manifest code/a5_closure/manifests/matter_menu_spectral_ledger_reference.json"
        ),
    }


def verify_manifest(stored: Mapping[str, Any], base_dir: Path | None = None) -> None:
    expected = certificate_payload(base_dir)
    require(
        stored == expected,
        "MANIFEST_MISMATCH",
        "the stored spectral-ledger manifest is stale, malformed, or tampered",
    )


def default_manifest_path() -> Path:
    return MODULE_DIR / "manifests" / "matter_menu_spectral_ledger_reference.json"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    emit = sub.add_parser("emit", help="derive and write the deterministic spectral-ledger manifest")
    emit.add_argument("--output", type=Path, default=default_manifest_path())
    verify = sub.add_parser("verify", help="recompute and compare a stored manifest")
    verify.add_argument("--manifest", type=Path, default=default_manifest_path())
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "emit":
        payload = certificate_payload()
        write_json(args.output, payload)
        print(
            json.dumps(
                {"status": "PASS", "manifest": str(args.output), "sha256": sha256_json(payload)},
                indent=2,
            )
        )
    else:
        stored = load_json(args.manifest)
        verify_manifest(stored)
        print(json.dumps({"status": "PASS", "manifest": str(args.manifest)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
