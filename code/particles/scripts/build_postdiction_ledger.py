#!/usr/bin/env python3
"""Aggregate the certified postdiction rows into one ledger.

The ledger is a deterministic aggregator.  Numeric values and measured
references are read live from their parent artifacts.  Structural rows are
derived from validated structured parents, with any direct algebraic
corollary identified as such.  A missing or inconsistent parent is a hard
failure, not a silently absent row.

Section one records the forced-structure layer: the machine-checked
icosahedral results that pin the gauge sector before any numeric lane runs.
The finite steps live in the Lean workspace under `Lean/Screen/`; the
receipt entries record the module paths and exact declaration names.  The
builder rejects a missing declaration and records the declared hypothesis
boundaries exactly as The Standard Model gauge paper states them.

The numeric sections carry the per-lane claim discipline of their parents:
interval rows report containment of the compare-only witness, conditional
rows carry their declared premises, chart coordinates keep their
NOT_EVALUABLE physical-comparison status, and the quark absolute-mass row
is an obstruction theorem rather than a number.

Run:
    python3 code/particles/scripts/build_postdiction_ledger.py
writes runs/status/postdiction_ledger.json and docs/POSTDICTION_LEDGER.md.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

SCRIPTS = Path(__file__).resolve().parent
PARTICLES = SCRIPTS.parent
CODE = PARTICLES.parent
REPO = CODE.parent
RUNS = PARTICLES / "runs"
RUNTIME = CODE / "P_derivation" / "runtime"
LEAN_SCREEN = REPO / "Lean" / "Screen"

PARENTS = {
    "mass_surface": RUNS / "status" / "source_only_mass_prediction_surface.json",
    "conditional_ew": RUNS / "calibration" / "conditional_ew_predictions_current.json",
    "endpoint": RUNTIME / "empirical_thomson_endpoint_current.json",
    "anchor_bridge": RUNTIME / "anchor_scheme_bridge_current.json",
    "kappa_rectangle": RUNS / "leptons" / "charged_kappa_interval_from_alpha_transport.json",
    "kappa_coherent": RUNS / "leptons" / "charged_kappa_interval_coherent_closure.json",
    "koide_balance": RUNS / "leptons" / "koide_balance_comparison.json",
    "clebsch_lane": RUNS / "flavor" / "down_type_register_clebsch_lane.json",
    "clebsch_selection": RUNS / "flavor" / "clebsch_register_pairing_selection.json",
    "fiber_obstruction": RUNS / "flavor" / "quark_spread_fiber_structure_transport_obstruction.json",
    "matter_receipt": CODE / "a5_closure" / "receipts" / "super_tannakian_matter_reference.receipt.json",
    "matter_menu": CODE / "a5_closure" / "manifests" / "matter_menu_spectral_ledger_reference.json",
    "port_current": CODE / "a5_closure" / "receipts" / "port_current_inner_reference.receipt.json",
    "axis_center_descent": CODE / "a5_closure" / "receipts" / "axis_center_descent_reference.receipt.json",
    "carrier_modes": RUNS / "status" / "carrier_mode_acceptance.json",
    "quantum_carrier_status": RUNS / "status" / "quantum_carrier_status.json",
    "alpha_hvp_verdict": PARTICLES / "alpha_hvp_audit" / "outputs" / "alpha_hvp_class_verdict.json",
    "hadron_payload": RUNS / "hadron" / "empirical_ee_hadronic_spectral_measure.json",
    "solver_standby": RUNS / "qcd" / "hadron_source_backend" / "qcd_ensemble" / "solver_on_standby.json",
}

LEAN_RECEIPTS = {
    "A2HolonomyBridge": LEAN_SCREEN / "A2HolonomyBridge.lean",
    "A5OPH": LEAN_SCREEN / "A5OPH.lean",
    "A5CharacterField": LEAN_SCREEN / "A5CharacterField.lean",
    "A5SixAxes": LEAN_SCREEN / "A5SixAxes.lean",
    "Z6Exact": LEAN_SCREEN / "Z6Exact.lean",
    "Z6Descent": LEAN_SCREEN / "Z6Descent.lean",
    "A5CouplingSymmetry": LEAN_SCREEN / "A5CouplingSymmetry.lean",
    "A5PortAction": LEAN_SCREEN / "A5PortAction.lean",
    "PortFrameGram": LEAN_SCREEN / "PortFrameGram.lean",
    "ExteriorSelection": LEAN_SCREEN / "ExteriorSelection.lean",
}

DEFAULT_OUT = RUNS / "status" / "postdiction_ledger.json"
DEFAULT_MD = REPO / "docs" / "POSTDICTION_LEDGER.md"


def _load(key: str, override: Path | None = None) -> dict[str, Any]:
    path = override or PARENTS[key]
    if not path.exists():
        raise SystemExit(f"postdiction ledger parent missing: {key} at {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _lean_receipt(
    *modules: str,
    declarations: dict[str, tuple[str, ...]] | None = None,
) -> list[str]:
    """Resolve Lean modules and fail closed on missing named declarations."""

    refs = []
    declarations = declarations or {}
    unknown = set(declarations) - set(modules)
    if unknown:
        raise SystemExit(
            "postdiction ledger declaration map names unbound modules: "
            + ", ".join(sorted(unknown))
        )
    for m in modules:
        path = LEAN_RECEIPTS[m]
        if not path.exists():
            raise SystemExit(f"postdiction ledger Lean receipt missing: {path}")
        source = path.read_text(encoding="utf-8")
        for declaration in declarations.get(m, ()):
            pattern = (
                rf"(?m)^\s*(?:theorem|lemma)\s+{re.escape(declaration)}"
                rf"(?:\s|:)"
            )
            if re.search(pattern, source) is None:
                raise SystemExit(
                    "postdiction ledger Lean declaration missing: "
                    f"{m}.{declaration} in {path}"
                )
        refs.append(path.relative_to(REPO).as_posix())
    return refs


def _rel(key: str) -> str:
    return PARENTS[key].relative_to(REPO).as_posix()


def _quantum_carrier_status_row(packet: dict[str, Any]) -> dict[str, Any]:
    core = {
        key: value for key, value in packet.items() if key != "receipt_sha256"
    }
    digest = "sha256:" + hashlib.sha256(
        json.dumps(
            core,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
    ).hexdigest()
    expected_verdicts = {
        "photon": "NOT_EVALUABLE_NO_SOURCE_SELECTED_MAXWELL_QUANTUM_SECTOR",
        "gluon": "NOT_EVALUABLE_NO_QCD",
        "graviton": "NOT_EVALUABLE_NO_INHABITED_EINSTEIN_QUANTUM_CARRIER",
    }
    rows = packet.get("rows")
    if not isinstance(rows, list):
        raise SystemExit("quantum-carrier status rows are absent")
    by_id = {
        row.get("carrier_id"): row for row in rows if isinstance(row, dict)
    }
    if (
        packet.get("schema") != "oph.quantum_carrier_status.v2"
        or packet.get("github_issue") != 552
        or packet.get("status") != "THREE_ROW_EXPLICIT_NOT_EVALUABLE"
        or packet.get("receipt_sha256") != digest
        or packet.get("comparison_values_consumed") is not False
        or packet.get("blind_prediction_eligible") is not False
        or packet.get("target_named_status_rows") is not True
        or packet.get("all_rows_at_allowed_exit") is not True
        or packet.get("continuum_spacetime_dimension") != 4
        or packet.get("classical_mode_vector_order")
        != ["photon", "gluon", "graviton"]
        or packet.get("classical_mode_vector") != [2, 16, 2]
        or set(by_id) != set(expected_verdicts)
    ):
        raise SystemExit("quantum-carrier status packet left its typed boundary")
    expected_multiplicities = {
        "photon": (1, "u1_lie_algebra_generator", 2),
        "gluon": (8, "su3_adjoint_generator", 16),
        "graviton": (1, "symmetric_metric_tensor_field", 2),
    }
    for carrier_id, verdict in expected_verdicts.items():
        row = by_id[carrier_id]
        capabilities = row.get("capabilities", {})
        baseline = row.get("classical_baseline", {})
        factor, role, total = expected_multiplicities[carrier_id]
        if (
            row.get("verdict") != verdict
            or row.get("verdict_class") != "EXPLICIT_NOT_EVALUABLE"
            or row.get("particle_promotion_allowed") is not False
            or not isinstance(row.get("blocking_frontier"), list)
            or not row["blocking_frontier"]
            or baseline.get("particle_claim") is not False
            or baseline.get("continuum_spacetime_dimension") != 4
            or baseline.get("multiplicity_factor") != factor
            or baseline.get("multiplicity_role") != role
            or baseline.get("exact_total_mode_count") != total
            or "gauge_algebra_dimension" in baseline
            or capabilities.get("state_space", {}).get(
                "physical_quantum_object_available"
            )
            is not False
            or capabilities.get("spectral_object", {}).get(
                "positive_physical_quantum_object_available"
            )
            is not False
            or capabilities.get("physical_current_residue", {}).get(
                "nonzero_positive_residue_available"
            )
            is not False
        ):
            raise SystemExit(
                f"quantum-carrier row {carrier_id} left its typed boundary"
            )
    return {
        "artifact_ref": _rel("quantum_carrier_status"),
        "classical_mode_vector": packet["classical_mode_vector"],
        "classical_mode_vector_order": packet["classical_mode_vector_order"],
        "receipt_sha256": packet["receipt_sha256"],
        "rows": [
            {
                "blocking_frontier": by_id[carrier_id]["blocking_frontier"],
                "carrier_id": carrier_id,
                "verdict": by_id[carrier_id]["verdict"],
            }
            for carrier_id in ("photon", "gluon", "graviton")
        ],
        "scope": (
            "The exact (2,16,2) vector contains differently typed conditional "
            "four-dimensional propagating-mode totals. Every quantum-particle "
            "row is explicitly not evaluable on the pinned declared corpus. "
            "The packet is target-named, comparison-value-free, and ineligible "
            "as a blind prediction."
        ),
    }


def _forced_structure(
    matter: dict[str, Any],
    matter_menu: dict[str, Any],
    port_current: dict[str, Any],
    axis_center_descent: dict[str, Any],
    carrier_modes: dict[str, Any],
) -> list[dict[str, Any]]:
    spectrum = matter["realized_package"]["charge_spectrum"]
    sm_spectrum = {"-1/2": 2, "-2/3": 3, "1": 1, "1/3": 3, "1/6": 6}
    scan = matter["matter_selection_scan"]
    classification = matter_menu["subset_classification"]
    scan_agreement = matter_menu["scan_agreement"]
    menu_verdicts = matter_menu["verdicts"]
    lean_cross_reference = classification["lean_cross_reference"]
    subsets_enumerated = classification["subsets_enumerated"]
    component_count = subsets_enumerated.bit_length() - 1
    survivors = classification["survivors"]
    survivor_dimensions = {row["dimension"] for row in survivors}
    if (
        2**component_count != subsets_enumerated
        or classification["survivor_count"] != len(survivors)
        or classification["survivor_count"] != scan["survivor_count"]
        or subsets_enumerated != scan["subsets_enumerated"]
        or not classification["survivors_are_conjugate_pair"]
        or not scan["survivors_are_conjugate_pair"]
        or survivor_dimensions != {matter["realized_package"]["dimension"]}
        or scan["derived_block_charges"]
        != matter_menu["declared_algebra"]["derived_block_charges"]
        or not scan_agreement["lean_masks_match_ledger"]
        or not scan_agreement["matter_lift_scan_matches_ledger"]
        or menu_verdicts["menu_completeness_inside_declared_algebra"] != "exact"
        or lean_cross_reference["agreement"] is not True
    ):
        raise SystemExit(
            "matter-menu, matter-lift, and Lean cross-reference parents disagree"
        )
    exterior_declarations = tuple(lean_cross_reference["theorems"])
    exterior_path = LEAN_RECEIPTS["ExteriorSelection"].relative_to(REPO).as_posix()
    if lean_cross_reference["file"] != exterior_path:
        raise SystemExit("matter-menu Lean source path does not match the ledger binding")
    realized_fields = matter["realized_package"]["fields"]
    field_order = ("Q", "u_c", "d_c", "L", "e_c")
    field_summary = ", ".join(
        f"{name}: {realized_fields[name]['charge']} x"
        f"{realized_fields[name]['dimension']}"
        for name in field_order
    )

    current_map = port_current["port_to_generator_map"]
    current_closure = port_current["closure"]
    derived_dimensions = current_closure["derived_block_dimensions"]
    color_adjoint_dimension = derived_dimensions["even_block_su3"]
    weak_adjoint_dimension = derived_dimensions["kernel_block_so3"]
    abelian_dimension = current_closure["center_dimension"]
    if (
        port_current["claim_boundary"]["status"]
        != "proved_conditional_on_declared_response_representation"
        or port_current["physical_source_gate"]["passed"] is not False
        or port_current["physical_source_gate"][
            "target_blind_impulse_readback_recomputed"
        ]
        is not True
        or port_current["source_definedness"][
            "response_model_declared_as_branch_premise"
        ]
        is not True
        or port_current["source_definedness"][
            "physical_response_source_bound"
        ]
        is not False
        or current_map["compact_lie_type"]
        != "u(3) (+) so(3) = u(1) (+) su(3) (+) su(2)"
        or current_map["block_dimensions_verified"] is not True
        or current_map["injective"] is not True
        or color_adjoint_dimension + weak_adjoint_dimension
        != current_closure["derived_dimension"]
        or (
            color_adjoint_dimension
            + weak_adjoint_dimension
            + abelian_dimension
            != current_map["image_real_dimension"]
        )
    ):
        raise SystemExit("port-current parent does not realize the declared product algebra")
    adjoint_branching = {
        "color_adjoint": [color_adjoint_dimension, 1, 0],
        "weak_adjoint": [1, weak_adjoint_dimension, 0],
        "abelian": [1, 1, 0],
        "mixed_xy_bifundamental_dimension": 0,
    }

    descent_gate = axis_center_descent["physical_global_form_gate"]
    declared_loop = axis_center_descent["carrier_deck_and_declared_loop_system"]
    tensor_kernel = axis_center_descent["kernel_on_realized_tensors"]
    effective_image = axis_center_descent["maximal_effective_image"]
    z6_bridge = axis_center_descent["two_z6_constructions"]
    if (
        axis_center_descent["claim_boundary"]["status"]
        != "conditional_exact_arithmetic_with_physical_global_form_open"
        or descent_gate["passed"] is not False
        or descent_gate["laboratory_global_form_attachment"] is not False
        or descent_gate["theta_periodicity_derived"] is not False
        or descent_gate["axis_relation_lattice_source_selected"] is not False
        or descent_gate["complete_character_category_source_derived"] is not False
        or descent_gate["same_source_loop_to_tensor_kernel_identification"]
        is not False
        or tensor_kernel["kernel_order"] != 6
        or tensor_kernel["matches_emitted_kernel_data"] is not True
        or effective_image["group"] != "(SU(3) x SU(2) x U(1)) / Z6"
        or declared_loop["six_axis_class_group_order"] != 6
        or declared_loop["declared_coefficient_system_menu"] != list(range(6))
        or declared_loop["axis_relation_lattice_source_selected"] is not False
        or z6_bridge["physical_loop_intertwiner_derived"] is not False
        or z6_bridge["conditional_algebraic_intertwiner_verified"] is not True
    ):
        raise SystemExit(
            "axis-centre descent parent has left its conditional arithmetic boundary"
        )

    carrier_by_id = {
        row["carrier_id"]: row for row in carrier_modes["carriers"]
    }
    if set(carrier_by_id) != {"photon", "gluon", "graviton"}:
        raise SystemExit("carrier-mode packet does not contain the expected rows")
    if not all(
        row["classical_carrier_gate"]["passed"]
        and not row["quantum_particle_gate"]["passed"]
        and row["hard_quadratic_mass_parameter_squared"] == 0
        for row in carrier_by_id.values()
    ):
        raise SystemExit("carrier-mode packet has left its declared boundary")

    rows = [
        {
            "id": "gauge_lie_algebra",
            "statement": (
                "The certified twelve-port carrier has module 1+3+3'+5 and "
                "one fixed line. Complete compact port response from A1 and "
                "endogenous proper-carrier transport from A2 force the abstract "
                "Lie type u(1)+su(2)+su(3). Target-blind impulse and readback "
                "separately determine R=-J. The charged-double-triplet matrices "
                "are an exact declared witness, while ordered source tomography "
                "and same-current holonomy remain open"
            ),
            "observed_counterpart": "Standard Model gauge Lie algebra su(3)+su(2)+u(1)",
            "match": "axiom-forced abstract Lie type; conditional matrix witness",
            "artifact_ref": _rel("port_current"),
            "machine_checked_steps": (
                "the carrier action and fixed-space dimension are exact; Lean "
                "checks the A2 holonomy-to-inner-action bridge, the centreless "
                "four-factor fixed-space exclusion, triviality of A5-actions "
                "on at most four objects, and unique "
                "partitions 11 = 8+3 and 12 = 3+3+3+3 over the compact-simple "
                "dimension list {3, 8, 10}; no compact semisimple algebra in "
                "dimensions 1, 2, 4, 5, 7; the characteristic-centre step; "
                "A5 not a subgroup of SU(2); the noncentrality witness "
                "[iS, iT] != 0; Galois stability of the two three-dimensional "
                "characters over Q(sqrt 5) with the centre-dimension list "
                "{0, 1, 5, 6, 7, 11, 12}; six-axis 2-transitivity and the "
                "dimension count of the dimension-six branch"
            ),
            "lean_declarations": {
                "A2HolonomyBridge": [
                    "internalImplementation_of_holonomy",
                    "four_factor_fixed_dimension_ne_one",
                    "compact_product_dimensions_of_fixed_space",
                ],
                "A5OPH": [
                    "sum_eq_eleven",
                    "sum_eq_twelve",
                    "action_trivial_of_card_le_four",
                    "sum_not_mem_excluded",
                    "quintet_noncentral",
                ],
                "A5CharacterField": [
                    "multiplicities_equal_of_galoisStable",
                    "centreDim_mem_trichotomy_list",
                ],
                "A5SixAxes": [
                    "two_transitive",
                    "V5_irreducible",
                    "no_three_plus_three_split",
                ],
            },
            "lean_receipts": _lean_receipt(
                "A2HolonomyBridge",
                "A5OPH",
                "A5CharacterField",
                "A5SixAxes",
                declarations={
                    "A2HolonomyBridge": (
                        "internalImplementation_of_holonomy",
                        "four_factor_fixed_dimension_ne_one",
                        "compact_product_dimensions_of_fixed_space",
                    ),
                    "A5OPH": (
                        "sum_eq_eleven",
                        "sum_eq_twelve",
                        "action_trivial_of_card_le_four",
                        "sum_not_mem_excluded",
                        "quintet_noncentral",
                    ),
                    "A5CharacterField": (
                        "multiplicities_equal_of_galoisStable",
                        "centreDim_mem_trichotomy_list",
                    ),
                    "A5SixAxes": (
                        "two_transitive",
                        "V5_irreducible",
                        "no_three_plus_three_split",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the abstract theorem uses the A1 faithful complete compact "
                "response and A2 endogenous holonomy clauses. The direct "
                "receipt verifies a declared charged-double-triplet witness "
                "and does not reconstruct it from ordered source histories or "
                "identify a laboratory current. Compact reductivity and the "
                "compact-simple classification are declared classical inputs"
            ),
            "paper_ref": "Standard Model gauge paper, Compact-Lie trichotomy section",
        },
        {
            "id": "global_form_z6",
            "statement": (
                "The common central kernel on every declared tensor is the "
                "order-six diagonal subgroup, so the maximal faithful image "
                "of that representation is (SU(3) x SU(2) x U(1))/Z6. "
                "The six-axis class has order six only after diagonal and "
                "zero-sum coefficient relations are declared. Source selection "
                "of those relations, a complete character category, and a "
                "same-source loop-to-kernel theorem remain open"
            ),
            "observed_counterpart": (
                "Standard Model global gauge-group form and its charge "
                "quantization pattern"
            ),
            "match": "exact conditional kernel and maximal faithful image",
            "artifact_refs": [
                _rel("matter_receipt"),
                _rel("axis_center_descent"),
            ],
            "lean_declarations": {
                "Z6Exact": [
                    "gauge_eq_kernel",
                    "residue_surjective",
                    "representative_formula",
                ],
                "Z6Descent": [
                    "kernel_on_realized_weights",
                    "four_admissible_global_forms",
                    "sixAxis_generator_maps_to_kernel_generator",
                    "sixAxisToKernel_intertwines_involutions",
                    "sixAxisToKernel_injective",
                    "sixAxisToKernel_range",
                ],
            },
            "lean_receipts": _lean_receipt(
                "Z6Exact",
                "Z6Descent",
                declarations={
                    "Z6Exact": (
                        "gauge_eq_kernel",
                        "residue_surjective",
                        "representative_formula",
                    ),
                    "Z6Descent": (
                        "kernel_on_realized_weights",
                        "four_admissible_global_forms",
                        "sixAxis_generator_maps_to_kernel_generator",
                        "sixAxisToKernel_intertwines_involutions",
                        "sixAxisToKernel_injective",
                        "sixAxisToKernel_range",
                    ),
                },
            ),
            "hypothesis_boundary": (
                "the result is exact for the declared matter table, central "
                "descent congruence, axis coefficient relations, and line "
                "lattice. The source current, physical matter action, relation-"
                "lattice selection, character completeness, loop-to-kernel "
                "identity, laboratory attachment, and continuum quantum field "
                "theory remain outside this result"
            ),
            "paper_ref": "Standard Model gauge paper, Z6 global-form section",
        },
        {
            "id": "hypercharge_spectrum",
            "statement": (
                f"Inside the declared {component_count}-component "
                "exterior-response algebra, an exhaustive scan of all "
                f"{subsets_enumerated} subsets selects exactly one "
                "unordered charge-conjugate pair of nonempty chiral "
                "anomaly-free rank-"
                f"{matter['realized_package']['dimension']} projectors. "
                "Primitive determinant balance "
                "fixes the block charges up to conjugation, and the selected "
                f"representative has multiset {{{field_summary}}}"
            ),
            "observed_counterpart": (
                "Standard Model one-generation hypercharge assignment"
            ),
            "realized_spectrum": spectrum,
            "match": "exact" if spectrum == sm_spectrum else "MISMATCH",
            "artifact_refs": [_rel("matter_receipt"), _rel("matter_menu")],
            "subset_count": subsets_enumerated,
            "survivor_count": classification["survivor_count"],
            "survivor_dimension": matter["realized_package"]["dimension"],
            "derived_block_charges": matter_menu["declared_algebra"][
                "derived_block_charges"
            ],
            "lean_declarations": {"ExteriorSelection": list(exterior_declarations)},
            "lean_receipts": _lean_receipt(
                "ExteriorSelection",
                declarations={"ExteriorSelection": exterior_declarations},
            ),
            "hypothesis_boundary": (
                "the selection is exhaustive inside the declared exterior "
                "algebra; completeness beyond that algebra, selection of one "
                "charge-conjugate representative, light-sector attachment, "
                "family multiplicity, scalar content, and laboratory "
                "identification remain separate"
            ),
            "paper_ref": "zoo paper, matter lift section",
        },
        {
            "id": "coupling_universality",
            "statement": (
                "A5-invariant readouts have port-independent group-averaged cap "
                "sums, so the per-cap ratio of any two averaged readouts is "
                "universal with zero spread"
            ),
            "observed_counterpart": (
                "universality clause of the Einstein-branch coupling law"
            ),
            "match": "structural",
            "lean_declarations": {
                "A5CouplingSymmetry": [
                    "groupAverage_port_independent",
                    "coupling_ratio_universal",
                ],
                "A5PortAction": ["transitive_on_ports"],
                "PortFrameGram": ["degree_five", "gram_sq"],
            },
            "lean_receipts": _lean_receipt(
                "A5CouplingSymmetry",
                "A5PortAction",
                "PortFrameGram",
                declarations={
                    "A5CouplingSymmetry": (
                        "groupAverage_port_independent",
                        "coupling_ratio_universal",
                    ),
                    "A5PortAction": ("transitive_on_ports",),
                    "PortFrameGram": ("degree_five", "gram_sq"),
                },
            ),
            "hypothesis_boundary": (
                "reduces the universality clause to A5-equivariance of the "
                "implemented source law; no coupling value is implied"
            ),
            "paper_ref": "Standard Model gauge paper, coupling symmetry section",
        },
    ]
    rows.extend(
        [
            {
                "id": "maxwell_classical_massless_kernel",
                "statement": (
                    "On the declared unbroken Maxwell action and deconfined "
                    "phase branch, the quadratic operator has zero hard mass "
                    "parameter and two transverse classical modes with "
                    "characteristic surface k^2=0"
                ),
                "observed_counterpart": (
                    "massless classical electromagnetic propagation"
                ),
                "match": "conditional structural",
                "artifact_ref": _rel("carrier_modes"),
                "hypothesis_boundary": (
                    "the Maxwell action, positive kinetic coefficient, field "
                    "content, and phase are supplied branch data; no photon "
                    "Hilbert space, positive-residue pole, or universal "
                    "zero-mass particle theorem is emitted"
                ),
                "paper_ref": "Observers paper, carrier-mode acceptance section",
            },
            {
                "id": "yang_mills_classical_massless_kernel",
                "statement": (
                    "On the declared pure Yang-Mills quadratic branch before "
                    "nonperturbative confinement, every color generator has "
                    "two transverse perturbative modes and zero hard "
                    "quadratic mass parameter"
                ),
                "observed_counterpart": (
                    "perturbative color-gauge kernel before confinement"
                ),
                "match": "conditional structural",
                "artifact_ref": _rel("carrier_modes"),
                "hypothesis_boundary": (
                    "this is not a free asymptotic-gluon claim and supplies "
                    "neither a continuum Yang-Mills gap nor a hadron mass"
                ),
                "paper_ref": "Observers paper, carrier-mode acceptance section",
            },
            {
                "id": "einstein_classical_massless_kernel",
                "statement": (
                    "On the declared pure Einstein-Hilbert linearization about "
                    "a suitable Ricci-flat background, the transverse-"
                    "traceless quadratic operator has zero hard mass parameter "
                    "and two classical modes with null characteristic"
                ),
                "observed_counterpart": (
                    "two massless classical gravitational-wave polarizations"
                ),
                "match": "conditional structural",
                "artifact_ref": _rel("carrier_modes"),
                "hypothesis_boundary": (
                    "the action and background are supplied branch data; no "
                    "graviton Hilbert space, quantum pole, or exclusion of "
                    "additional massive modes is emitted"
                ),
                "paper_ref": "Observers paper, carrier-mode acceptance section",
            },
            {
                "id": "simple_gut_xy_channel_absent",
                "statement": (
                    "The declared charged-double-triplet current fixture has "
                    "a direct-sum algebra with adjoint "
                    f"branch dimensions {color_adjoint_dimension}, "
                    f"{weak_adjoint_dimension}, and {abelian_dimension}. Its "
                    "adjoint therefore contains no mixed "
                    "(3,2,-5/6) (+) (bar3,2,+5/6) X/Y generator, so the "
                    "ordinary minimal simple-GUT X/Y exchange channel is absent"
                ),
                "observed_counterpart": (
                    "the Standard Model product adjoint contains no connected "
                    "simple-GUT X/Y generator"
                ),
                "match": "conditional algebraic channel exclusion",
                "artifact_ref": _rel("port_current"),
                "operator_census_ref": (
                    "code/a5_closure/receipts/"
                    "baryon_dimension_six_census.receipt.json"
                ),
                "derivation_kind": "direct_executable_algebraic_corollary",
                "adjoint_branching": adjoint_branching,
                "hypothesis_boundary": (
                    "the executable corollary applies to the declared "
                    "direct-sum matrix-current fixture. Its physical current "
                    "source gate is false, so the result is not a physical "
                    "current or proton-stability claim. General proton "
                    "stability does not follow. Conditional on the "
                    "declared one-generation matter table and baryon and "
                    "lepton labels, an exact dimension-six census admits "
                    "QQQL, QQUE, DUQL, and DUUE; the representatives remain "
                    "nonzero after the exterior-algebra relations. No "
                    "coefficient, physical decay amplitude, QCD matrix "
                    "element, or lifetime is supplied"
                ),
                "paper_ref": "Observers paper, gauge-channel boundary",
            },
        ]
    )
    return rows


def _alpha_rows(
    endpoint: dict[str, Any],
    bridge: dict[str, Any],
    alpha_hvp_verdict: dict[str, Any],
) -> list[dict[str, Any]]:
    ep = endpoint["endpoint"]
    co = endpoint["compare_only"]
    bridge_verdict = bridge["verdict"]
    if "reference_deficit_inside_certified_gap" not in bridge_verdict:
        raise SystemExit("anchor bridge artifact lacks the containment verdict field")
    cross_class = alpha_hvp_verdict["cross_class_agreement"]
    scope = alpha_hvp_verdict["scope"]
    if (
        alpha_hvp_verdict["verdict"]
        != "MULTI_CLASS_NOT_EVALUABLE__ONE_RECORDED_ACCOUNTING_REPLAY_COMPATIBLE"
        or cross_class["recorded_accounting_replay_count"] != 1
        or cross_class["independently_evaluated_class_count"] != 0
        or cross_class["verdict"] != "NOT_EVALUABLE_NO_INDEPENDENT_CLASS"
        or scope["comparison_timing"] != "retrospective"
        or scope["prospective_freeze"] is not False
        or scope["physical_alpha_prediction_emitted"] is not False
    ):
        raise SystemExit("alpha/HVP verdict has left its retrospective audit boundary")
    containment = bridge_verdict["reference_deficit_inside_certified_gap"]
    return [
        {
            "id": "alpha_inv_thomson_endpoint",
            "value_central": float(ep["alpha_inv_central"]),
            "value_interval": [float(v) for v in ep["alpha_inv_interval"]],
            "measured": float(co["codata_alpha_inv"]),
            "measured_source": "CODATA 2022 via the endpoint artifact, compare-only",
            "payload_release": endpoint["inputs"]["payload_release"],
            "row_class": endpoint["row_class"],
            "tier": alpha_hvp_verdict["row_class"],
            "anchor_gap_interval": [float(v) for v in co["same_scheme_anchor_gap_interval_inv_alpha"]],
            "reference_deficit_inside_recorded_accounting_interval": containment,
            "audit_verdict": alpha_hvp_verdict["verdict"],
            "cross_class_agreement": cross_class,
            "reading": (
                "one retrospective KNT19 accounting row is arithmetically "
                "compatible with the recorded same-scheme interval. The "
                "multi-class HVP test is not evaluable because no independent "
                "frozen class is present. Containment does not identify the "
                "physical source of the gap or close source-only transport"
            ),
            "artifact_refs": [
                _rel("endpoint"),
                _rel("anchor_bridge"),
                _rel("alpha_hvp_verdict"),
            ],
            "blocking_issues": [425, 545],
        }
    ]


def _lepton_rows(
    surface: dict[str, Any],
    rectangle: dict[str, Any],
    coherent: dict[str, Any],
    koide: dict[str, Any],
) -> list[dict[str, Any]]:
    witness_point = rectangle["compare_only"].get("witness_point")
    if witness_point is None:
        raise SystemExit(
            "rectangle artifact lacks the witness_point block; rebuild the "
            "rectangle lane first"
        )
    width_floor = coherent.get("width_floor_audit")
    if width_floor is None:
        raise SystemExit(
            "coherent artifact lacks the width_floor_audit block; rebuild "
            "the coherent lane first"
        )
    witnesses = rectangle["compare_only"]["witness_masses_gev"]
    particles = [r["particle"] for r in rectangle["conditional_mass_rows"]]
    family = next(f for f in surface["families"] if f["family"] == "charged leptons")
    mcpr = next(r for r in family["rows"] if r["lane"].startswith("MCPR"))
    rows: list[dict[str, Any]] = []
    rows.append(
        {
            "id": "charged_leptons_closure_target",
            "statement": (
                "The certified solve inverts exactly at the measured triple: "
                "one anchor-gap value closes the charged-lepton lane on the "
                "witness, that value lies inside the retrospective accounting "
                "interval, and its "
                "distance to the standard on-shell reference deficit is the "
                "live scheme term of the open anchor bridge. The lepton "
                "scale is localized only under the recorded accounting "
                "packet. A "
                "source-emitted bridge value is a sharp falsification "
                "target: landing on the closure value satisfies the "
                "conditional lane on the witness, while landing outside the "
                "recorded interval refutes "
                "the decomposition."
            ),
            "witness_point": witness_point,
            "width_floor": width_floor["floor_attribution"],
            "tier": "T1_empirical_closure",
            "artifact_refs": [_rel("kappa_rectangle"), _rel("kappa_coherent")],
            "blocking_issues": [545, 425],
        }
    )
    mcpr_masses = [float(m) / 1000.0 for m in mcpr["masses_MeV_display"]]
    rows.append(
        {
            "id": "charged_leptons_mcpr_conditional",
            "particles": particles,
            "values_gev": mcpr_masses,
            "measured_gev": witnesses,
            "measured_source": "PDG witness triple embedded in the kappa lane, compare-only",
            "relative_deltas": [v / w - 1.0 for v, w in zip(mcpr_masses, witnesses, strict=True)],
            "tier": mcpr["tier"],
            "row_class": mcpr["row_class"],
            "explanation": mcpr["explanation"],
            "artifact_ref": mcpr["artifact"],
            "blocking_objects": mcpr["blocking_objects"],
        }
    )
    for key, lane, artifact in (
        ("charged_leptons_kappa_rectangle", rectangle, "kappa_rectangle"),
        ("charged_leptons_kappa_coherent", coherent, "kappa_coherent"),
    ):
        mass_rows = lane["conditional_mass_rows"]
        containment = all(
            row["mass_interval"][0] < w < row["mass_interval"][1]
            for row, w in zip(mass_rows, witnesses, strict=True)
        )
        entry: dict[str, Any] = {
            "id": key,
            "particles": particles,
            "intervals_gev": [row["mass_interval"] for row in mass_rows],
            "centrals_gev": [row["mass_central"] for row in mass_rows],
            "measured_gev": witnesses,
            "measured_source": "PDG witness triple embedded in the lane, compare-only",
            "witness_inside_all_intervals": containment,
            "relative_half_widths": [
                (row["mass_interval"][1] - row["mass_interval"][0])
                / (2.0 * row["mass_central"])
                for row in mass_rows
            ],
            "logarithmic_half_width": (
                lane["kappa_interval"]["interval"][1]
                - lane["kappa_interval"]["interval"][0]
            )
            / 2.0,
            "one_sided_multiplicative_widths": {
                "lower": 1.0
                - math.exp(
                    -(
                        lane["kappa_interval"]["interval"][1]
                        - lane["kappa_interval"]["interval"][0]
                    )
                    / 2.0
                ),
                "upper": math.exp(
                    (
                        lane["kappa_interval"]["interval"][1]
                        - lane["kappa_interval"]["interval"][0]
                    )
                    / 2.0
                )
                - 1.0,
            },
            "tier": "T1_empirical_closure",
            "row_class": lane["row_class"],
            "epistemic_scope": lane["numerical_certificate"]["epistemic_scope"],
            "numerical_certificate": lane["numerical_certificate"],
            "artifact_ref": _rel(artifact),
            "blocking_issues": [425, 545],
        }
        if key.endswith("coherent"):
            entry["width_reduction_factor"] = lane["kappa_interval"]["width_reduction_factor"]
            entry["premise"] = "payload-coherent anchor-gap premise, declared"
        rows.append(entry)

    tau_row = koide["conditional_tau"]
    balance = koide["balance_comparison"]
    rows.append(
        {
            "id": "charged_leptons_koide_conditional_tau",
            "premises": tau_row["premises"],
            "inputs": tau_row["inputs"],
            "tau_enclosure_mev_outward": tau_row["tau_enclosure_mev_outward"],
            "tau_central_mev": tau_row["tau_central_mev"],
            "measured_tau_mev": tau_row["measured_tau_mev"],
            "distance_sigma": tau_row["distance_sigma"],
            "spurious_root_excluded_by": tau_row["spurious_root_excluded_by"],
            "balance_target_inside_enclosure": balance["target_inside_enclosure"],
            "balance_distance_half_widths": balance["distance_in_enclosure_half_widths"],
            "forward_test": (
                "the enclosure is three orders of magnitude narrower than "
                "the measurement uncertainty; an improving tau-mass average "
                "outside the window refutes the balanced-circulant premise"
            ),
            "tier": "T2_conditional",
            "row_class": tau_row["row_class"],
            "artifact_ref": _rel("koide_balance"),
        }
    )
    return rows


def _ew_rows(conditional: dict[str, Any]) -> list[dict[str, Any]]:
    cc = conditional["comparison_compare_only"]
    rows = []
    for key in ("mH_gev", "mt_pole_gev", "MW_chart_gev", "MZ_chart_gev"):
        block = cc[key]
        row: dict[str, Any] = {
            "id": f"ew_{key}",
            "value_central": block["conditional_central"],
            "value_envelope": block["conditional_envelope"],
            "physical_comparison_status": block["physical_comparison_status"],
            "tier": "T2_conditional",
            "row_class": conditional["row_class"],
            "artifact_ref": _rel("conditional_ew"),
        }
        if block["physical_comparison_status"] == "COMPARE_ONLY":
            row.update(
                {
                    "measured": block["measured"],
                    "measured_sigma": block["measured_sigma"],
                    "measured_source": block["measured_source"],
                    "delta": block["delta"],
                    "delta_over_sigma": block["delta_over_sigma"],
                    "envelope_overlaps_one_sigma_band": block["envelope_overlaps_one_sigma_band"],
                }
            )
        else:
            row["reason"] = block["reason"]
        rows.append(row)
    return rows


def _quark_rows(
    obstruction: dict[str, Any],
    clebsch: dict[str, Any],
    selection: dict[str, Any],
) -> list[dict[str, Any]]:
    if obstruction["fork"] != "ii_fiber_survives":
        raise SystemExit(
            "fiber obstruction artifact is not on the survives fork; "
            "rebuild the quark section before aggregating"
        )
    compare = clebsch["compare_only"]
    predictions = clebsch["predictions"]
    return [
        {
            "id": "quark_absolute_masses_obstruction",
            "statement": (
                "No absolute quark mass is emitted: the two-modulus spread "
                "fiber survives every certified structure transport, so the "
                "six absolute masses are non-identifiable from the corpus, by "
                "theorem rather than by omission"
            ),
            "fork": obstruction["fork"],
            "fiber_cut_detected": obstruction["fiber_cut_detected"],
            "tier": obstruction["claim_tier"],
            "artifact_ref": _rel("fiber_obstruction"),
            "blocking_issues": obstruction["github_issues"],
        },
        {
            "id": "quark_down_type_clebsch_route_rejected",
            "values": predictions,
            "measured_references": compare["references"],
            "relative_deltas": {
                k: v for k, v in compare.items() if k.endswith("_relative")
            },
            "flag_2024_compare_only": compare["flag_2024"],
            "retrospective_flag_rejection": clebsch[
                "retrospective_flag_rejection"
            ],
            "permutation_scan": clebsch["permutation_scan"],
            "promotion_allowed": clebsch["promotion_allowed"],
            "status": clebsch["status"],
            "tier": "T2_conditional_rejected_candidate",
            "row_class": clebsch["row_class"],
            "premise": (
                "a cross-sector register relation, independent Yukawa "
                "coefficient identification, and a physical generation order; "
                "the pairing receipt supplies channel compatibility only"
            ),
            "selection_artifact_ref": _rel("clebsch_selection"),
            "selection_status": selection["status"],
            "artifact_ref": _rel("clebsch_lane"),
            "reading": (
                "The target-free F1/F2 scan fixes only the unordered multiset. "
                "All six assignments fail the retrospective conservative FLAG "
                "gate. The displayed GST value is sqrt(md/ms) under an assumed "
                "texture, not a derived CKM angle; a simultaneous diagonal "
                "mass ansatz would instead give the identity CKM matrix."
            ),
        },
    ]


def _hadron_rows(payload: dict[str, Any], standby: dict[str, Any]) -> list[dict[str, Any]]:
    integral = payload["integral"]
    norm = integral["normalization"]
    return [
        {
            "id": "hadronic_correction_engine",
            "delta_alpha_had_5_MZ": integral["value"],
            "uncertainty_total": integral["uncertainty"],
            "source_compilation": payload["source_compilation"]["id"],
            "pin_factor": norm["pin_factor"],
            "policy": (
                "The published-compilation payload is the correction engine of "
                "the fine-structure lane; source-only hadron rows stay "
                "suppressed pending the source spectral measure (issue 425)."
            ),
            "artifact_ref": _rel("hadron_payload"),
        },
        {
            "id": "qcd_solver_on_standby",
            "status": standby["status"],
            "invocation_gate": standby["policy"]["invocation_gate"],
            "artifact_ref": _rel("solver_standby"),
        },
    ]


def _principal_results(sections: dict[str, Any]) -> list[dict[str, Any]]:
    """Digest the four strongest rows into the leading section, data-driven."""

    leptons = {r["id"]: r for r in sections["charged_leptons"]}
    target = leptons["charged_leptons_closure_target"]
    coherent = leptons["charged_leptons_kappa_coherent"]
    mcpr = leptons["charged_leptons_mcpr_conditional"]
    ew = {r["id"]: r for r in sections["electroweak"]}
    alpha = sections["alpha"][0]
    wp = target["witness_point"]
    glo, ghi = alpha["anchor_gap_interval"]
    log_hw = coherent["logarithmic_half_width"]
    one_sided = coherent["one_sided_multiplicative_widths"]
    ppm = abs(mcpr["relative_deltas"][0]) * 1.0e6
    mh, mt = ew["ew_mH_gev"], ew["ew_mt_pole_gev"]
    return [
        {
            "id": "lepton_closure_target",
            "statement": (
                "The anchor-gap value "
                f"{wp['required_anchor_gap_at_witness_inv_alpha']:.4f} closes "
                "the charged-lepton lane exactly on the measured triple, "
                f"inside the retrospective accounting interval "
                f"[{glo:.4f}, {ghi:.4f}]; the "
                f"distance {wp['scheme_term_difference_inv_alpha']:+.4f} to "
                "the standard on-shell reference deficit "
                f"{wp['reference_deficit_inv_alpha']:.4f} is the live scheme "
                "term of the open anchor bridge (issue 545). The lepton "
                "scale is localized only under that recorded accounting "
                "packet. A "
                "source-emitted bridge value is a falsification target: the "
                "closure value would satisfy the conditional lane, while a "
                "value outside the interval refutes the declared decomposition."
            ),
        },
        {
            "id": "lepton_certified_intervals",
            "statement": (
                "The target-anchored measured charged-lepton triple lies "
                "inside every outward-rounded diagnostic interval; the "
                "payload-coherent logarithmic half-width is "
                f"{log_hw * 100.0:.3f} percent, with one-sided multiplicative "
                f"widths -{one_sided['lower'] * 100.0:.2f} and "
                f"+{one_sided['upper'] * 100.0:.2f} percent. The conditional "
                f"eight-register triple sits {ppm:.0f} ppm from measurement "
                "with the architecture declared."
            ),
        },
        {
            "id": "koide_conditional_tau_window",
            "statement": (
                "Under the balanced-circulant and mass-ordering premises the "
                "measured electron and muon masses fix the tau mass inside "
                f"[{leptons['charged_leptons_koide_conditional_tau']['tau_enclosure_mev_outward'][0]}, "
                f"{leptons['charged_leptons_koide_conditional_tau']['tau_enclosure_mev_outward'][1]}] MeV, "
                f"{leptons['charged_leptons_koide_conditional_tau']['distance_sigma']} sigma from "
                "measurement; the window is three orders of magnitude "
                "narrower than the measurement uncertainty, so improving "
                "tau-mass averages test the premise directly. The premise "
                "ancestry is declared and the row stays conditional."
            ),
        },
        {
            "id": "higgs_top_envelopes",
            "statement": (
                f"The conditional Higgs envelope [{mh['value_envelope'][0]:.3f}, "
                f"{mh['value_envelope'][1]:.3f}] GeV sits "
                f"{mh['delta_over_sigma']:.2f} sigma from the measured "
                f"{mh['measured']} +- {mh['measured_sigma']} GeV, and the top "
                f"envelope [{mt['value_envelope'][0]:.2f}, "
                f"{mt['value_envelope'][1]:.2f}] GeV sits "
                f"{mt['delta_over_sigma']:.2f} sigma from "
                f"{mt['measured']} +- {mt['measured_sigma']} GeV, "
                "compare-only, conditional on the declared selection premises."
            ),
        },
        {
            "id": "forced_gauge_structure",
            "statement": (
                "Complete compact port response from A1 and endogenous overlap "
                "transport from A2 force the abstract Lie type "
                "su(3)+su(2)+u(1) on the twelve-port carrier. Target-blind "
                "readback independently derives R=-J. Inside the declared "
                "exterior-response algebra, an exhaustive scan selects the "
                "charge-conjugate rank-15 chiral anomaly-free pair and its "
                "one-generation hypercharge multiset; its common central "
                "kernel is Z6. Source reconstruction of the matrix current and "
                "matter action, physical global-form selection, laboratory "
                "attachment, and continuum quantum field theory remain separate."
            ),
        },
    ]


def build(
    out_path: Path = DEFAULT_OUT,
    md_path: Path | None = DEFAULT_MD,
    *,
    write: bool = True,
) -> dict[str, Any]:
    surface = _load("mass_surface")
    conditional = _load("conditional_ew")
    endpoint = _load("endpoint")
    bridge = _load("anchor_bridge")
    rectangle = _load("kappa_rectangle")
    coherent = _load("kappa_coherent")
    koide = _load("koide_balance")
    clebsch = _load("clebsch_lane")
    selection = _load("clebsch_selection")
    obstruction = _load("fiber_obstruction")
    matter = _load("matter_receipt")
    matter_menu = _load("matter_menu")
    port_current = _load("port_current")
    axis_center_descent = _load("axis_center_descent")
    carrier_modes = _load("carrier_modes")
    quantum_carrier_status = _load("quantum_carrier_status")
    alpha_hvp_verdict = _load("alpha_hvp_verdict")
    payload = _load("hadron_payload")
    standby = _load("solver_standby")

    sections = {
        "forced_structure": _forced_structure(
            matter,
            matter_menu,
            port_current,
            axis_center_descent,
            carrier_modes,
        ),
        "quantum_carrier_status": _quantum_carrier_status_row(
            quantum_carrier_status
        ),
        "alpha": _alpha_rows(endpoint, bridge, alpha_hvp_verdict),
        "charged_leptons": _lepton_rows(surface, rectangle, coherent, koide),
        "electroweak": _ew_rows(conditional),
        "quarks": _quark_rows(obstruction, clebsch, selection),
        "hadrons": _hadron_rows(payload, standby),
        "neutrinos": [
            {
                "id": "neutrino_dimensionless_pointer",
                "statement": (
                    "dimensionless PMNS and mass-splitting-ratio "
                    "comparisons live on the results status surface; the "
                    "absolute attachment stays compare-only"
                ),
                "artifact_ref": "code/particles/RESULTS_STATUS.md",
            }
        ],
    }
    result = {
        "artifact": "oph_postdiction_ledger",
        "generator": "code/particles/scripts/build_postdiction_ledger.py",
        "schema_version": 2,
        "row_class": "compare_only_postdiction_ledger",
        "guards": {
            "compare_only": True,
            "public_promotion_allowed": False,
            "changes_any_solve_path": False,
            "new_axiom_introduced": False,
            "hand_typed_measured_values": False,
        },
        "aggregation_policy": (
            "numeric values and measured references are read live from cited "
            "parents; structural rows are derived from validated structured "
            "parents and identify direct algebraic corollaries explicitly; a "
            "missing or inconsistent parent aborts the build"
        ),
        "principal_results": _principal_results(sections),
        "sections": sections,
    }
    if write:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        if md_path is not None:
            md_path.write_text(_render_md(result), encoding="utf-8")
    return result


def _fmt(x: float, digits: int = 6) -> str:
    return f"{x:.{digits}g}"


def _render_md(ledger: dict[str, Any]) -> str:
    s = ledger["sections"]
    lines: list[str] = []
    add = lines.append
    add("# Postdiction Ledger")
    add("")
    add("Generated deterministically by `scripts/build_postdiction_ledger.py`; "
        "the JSON artifact is `runs/status/postdiction_ledger.json`.")
    add("")
    add("Numeric values and measured references on this page are read live from "
        "the cited parent artifacts. Structural rows are derived from validated "
        "structured parents, and direct algebraic corollaries are identified. "
        "The ledger promotes nothing and changes no solve path. Interval rows "
        "report containment of the compare-only witness; conditional rows carry "
        "their declared premises; chart coordinates keep their NOT_EVALUABLE "
        "physical-comparison status.")
    add("")
    add("## Principal results")
    add("")
    for entry in ledger["principal_results"]:
        add(f"- {entry['statement']}")
    add("")
    add("## Forced structure")
    add("")
    add("The icosahedral screen results pin the gauge sector before any numeric "
        "lane runs. The finite steps are machine checked in the Lean workspace; "
        "the recorded hypothesis boundaries are the exact classical inputs and "
        "open premises of The Standard Model gauge paper.")
    add("")
    add("| Result | Observed counterpart | Match | Receipts |")
    add("| --- | --- | --- | --- |")
    for row in s["forced_structure"]:
        receipts = list(row.get("lean_receipts", []))
        if row.get("artifact_ref"):
            receipts.append(row["artifact_ref"])
        receipts.extend(row.get("artifact_refs", []))
        receipt_txt = ", ".join(f"`{r}`" for r in receipts)
        add(f"| {row['statement']} | {row['observed_counterpart']} | "
            f"`{row['match']}` | {receipt_txt} |")
    add("")
    add("Lean declaration bindings:")
    add("")
    for row in s["forced_structure"]:
        bindings = row.get("lean_declarations")
        if not bindings:
            continue
        rendered = "; ".join(
            f"`{module}`: " + ", ".join(f"`{name}`" for name in names)
            for module, names in bindings.items()
        )
        add(f"- `{row['id']}`: {rendered}")
    add("")
    add("Hypothesis boundaries:")
    add("")
    for row in s["forced_structure"]:
        add(f"- `{row['id']}`: {row['hypothesis_boundary']}")
    add("")
    add("## Quantum carrier gate")
    add("")
    carrier_status = s["quantum_carrier_status"]
    add(
        "The exact conditional four-dimensional propagating-mode vector is "
        "`(2, 16, 2)`: two Maxwell modes from one U(1) generator, sixteen "
        "perturbative color modes from eight SU(3) adjoint generators, and two "
        "Einstein transverse-traceless modes from one metric tensor field. "
        "The entries are neither particle counts nor one uniform gauge-algebra "
        "dimension vector."
    )
    add("")
    add("| Carrier | Quantum verdict | Blocking frontier |")
    add("| --- | --- | --- |")
    for row in carrier_status["rows"]:
        blockers = ", ".join(f"`{item}`" for item in row["blocking_frontier"])
        add(
            f"| `{row['carrier_id']}` | `{row['verdict']}` | {blockers} |"
        )
    add("")
    add(
        "The target-named status packet consumes no laboratory comparison value, "
        "permits no particle promotion, and is ineligible as a blind prediction. "
        "Its receipt is "
        f"`{carrier_status['artifact_ref']}`."
    )
    add("")
    add("## Fine-structure lane")
    add("")
    for row in s["alpha"]:
        lo, hi = row["value_interval"]
        glo, ghi = row["anchor_gap_interval"]
        add(f"- `alpha_em^-1` Thomson endpoint: `{_fmt(row['value_central'], 10)}` "
            f"in `[{_fmt(lo, 10)}, {_fmt(hi, 10)}]` against CODATA "
            f"`{_fmt(row['measured'], 10)}` (compare-only). Payload release "
            f"`{row['payload_release']}`.")
        inside = (
            "inside"
            if row["reference_deficit_inside_recorded_accounting_interval"]
            else "outside"
        )
        add(f"- Recorded retrospective same-scheme accounting interval "
            f"`[{_fmt(glo, 4)}, {_fmt(ghi, 4)}]` inverse-alpha units; the "
            f"standard reference deficit sits {inside} that interval.")
        add(f"- Independent-class verdict: `{row['audit_verdict']}`; evaluated "
            f"independent classes: "
            f"`{row['cross_class_agreement']['independently_evaluated_class_count']}`.")
        add(f"- Reading: {row['reading']}")
        add(f"- Blocking issues: {', '.join(f'#{i}' for i in row['blocking_issues'])}")
    add("")
    add("## Charged leptons")
    add("")
    for row in s["charged_leptons"]:
        if row["id"].endswith("closure_target"):
            wp = row["witness_point"]
            add(f"- Closure target ({row['tier']}): the anchor-gap value "
                f"`{wp['required_anchor_gap_at_witness_inv_alpha']:.4f}` closes the "
                "lane exactly on the measured triple (inversion machine-checked); "
                f"the distance `{wp['scheme_term_difference_inv_alpha']:+.4f}` to the "
                f"on-shell reference deficit `{wp['reference_deficit_inv_alpha']:.4f}` "
                "is the live scheme term of the bridge. The certified width floor "
                "is the scheme-band ambiguity; no budget is shrunk without the "
                "source bridge.")
            continue
        if row["id"].endswith("koide_conditional_tau"):
            lo, hi = row["tau_enclosure_mev_outward"]
            measured, sigma = row["measured_tau_mev"]
            add(f"- Koide conditional tau ({row['tier']}): under the "
                "balanced-circulant and mass-ordering premises the measured "
                "electron and muon masses fix the tau mass inside "
                f"`[{lo}, {hi}]` MeV, `{row['distance_sigma']}` sigma from "
                f"the measured `{measured} +- {sigma}` MeV; the premise "
                "ancestry is declared and improving tau-mass averages test "
                "the premise directly.")
            continue
        if row["id"].endswith("mcpr_conditional"):
            deltas = ", ".join(
                f"{p} `{_fmt(d * 1e6, 3)} ppm`"
                for p, d in zip(row["particles"], row["relative_deltas"], strict=True)
            )
            add(f"- MCPR conditional triple ({row['tier']}): {deltas} against the "
                "PDG witness triple; the eight-register architecture is a "
                "declared model input.")
        else:
            kind = "coherent closure" if row["id"].endswith("coherent") else "rectangle"
            contained = "inside" if row["witness_inside_all_intervals"] else "OUTSIDE"
            one_sided = row["one_sided_multiplicative_widths"]
            add(f"- Kappa interval, {kind} ({row['tier']}): outward-rounded "
                "target-anchored diagnostic intervals with logarithmic half-width "
                f"`{_fmt(row['logarithmic_half_width'] * 100, 4)}%` and one-sided "
                f"multiplicative widths `-{_fmt(one_sided['lower'] * 100, 3)}%` / "
                f"`+{_fmt(one_sided['upper'] * 100, 3)}%`; the witness triple lies "
                f"{contained} every interval.")
            if "width_reduction_factor" in row:
                add(f"  - Width reduction over the rectangle: "
                    f"`{_fmt(row['width_reduction_factor'], 3)}x`; premise: "
                    f"{row['premise']}.")
    add("")
    add("## Electroweak sector")
    add("")
    add("| Quantity | Conditional central | Envelope | Measured | Delta/sigma | Status |")
    add("| --- | ---: | --- | --- | ---: | --- |")
    for row in s["electroweak"]:
        env = f"[{_fmt(row['value_envelope'][0], 8)}, {_fmt(row['value_envelope'][1], 8)}]"
        if row["physical_comparison_status"] == "COMPARE_ONLY":
            add(f"| `{row['id'][3:]}` | `{_fmt(row['value_central'], 8)}` | `{env}` | "
                f"`{row['measured']} +- {row['measured_sigma']}` ({row['measured_source']}) | "
                f"`{_fmt(row['delta_over_sigma'], 3)}` | compare-only |")
        else:
            add(f"| `{row['id'][3:]}` | `{_fmt(row['value_central'], 8)}` | `{env}` | "
                "chart coordinate | n/a | NOT_EVALUABLE |")
    add("")
    add("W/Z rows are running/tree chart coordinates. The strict one-loop "
        "consumer has a separate external fixture: interval receipts exclude "
        "scalar zeros in the declared principal-sheet boxes and isolate, for "
        "each of W and Z, one simple scalar zero with derivative and scalar-residue "
        "balls in its declared lower-half pole box on a channel-specific algebraic "
        "chart. They identify neither chart with the physical resonance sheet and "
        "prove no unique continuation, sign bridge, full-matrix Laurent "
        "residue, physical-current amplitude, or independent numerical replay. "
        "The fixture is not composed with the OPH chart, so no physical W/Z pole "
        "or mass comparison is defined. The Higgs and top rows are conditional "
        "on the declared selection axioms.")
    add("")
    add("## Quarks")
    add("")
    for row in s["quarks"]:
        if row["id"].endswith("obstruction"):
            add(f"- Absolute masses ({row['tier']}): {row['statement']} "
                f"(issues {', '.join(f'#{i}' for i in row['blocking_issues'])}).")
        else:
            vals = row["values"]
            refs = row["measured_references"]
            flag = row["flag_2024_compare_only"]
            flag_refs = ", ".join(
                f"Nf={entry['nf']}: {_fmt(entry['reference_ms_over_md'], 4)}"
                for entry in flag
            )
            add(
                f"- Down-type register-Clebsch route, rejected "
                f"({row['tier']}): `ms/md = {_fmt(vals['ms_over_md'], 4)}` "
                f"against FLAG 2024 ({flag_refs}); all six generation "
                f"assignments are rejected by the retrospective conservative "
                f"gate. The diagnostic `sqrt(md/ms) = "
                f"{_fmt(vals['cabibbo_gst_sqrt_md_over_ms'], 4)}` is not a "
                f"derived Cabibbo angle. Premise: {row['premise']}. "
                f"{row['reading']}"
            )
    add("")
    add("## Hadrons")
    add("")
    for row in s["hadrons"]:
        if row["id"] == "hadronic_correction_engine":
            add(f"- Correction engine payload: `Delta alpha_had^(5)(M_Z^2) = "
                f"{row['delta_alpha_had_5_MZ']} +- {row['uncertainty_total']}` "
                f"from `{row['source_compilation']}` "
                f"(pin factor `{_fmt(row['pin_factor'], 7)}`). {row['policy']}")
        else:
            add(f"- QCD solver: `{row['status']}`; invocation is gated on the "
                "source-side parameter emissions recorded in the standby receipt.")
    add("")
    add("## Neutrinos")
    add("")
    for row in s["neutrinos"]:
        add(f"- {row['statement']} (`{row['artifact_ref']}`).")
    add("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--md", type=Path, default=DEFAULT_MD)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if the committed JSON or Markdown differs from a live rebuild",
    )
    args = parser.parse_args()
    result = build(args.out, args.md, write=not args.check)
    if args.check:
        expected_json = json.dumps(result, indent=2, sort_keys=True) + "\n"
        expected_md = _render_md(result)
        if not args.out.is_file() or args.out.read_text(encoding="utf-8") != expected_json:
            raise SystemExit(f"postdiction ledger JSON drift: {args.out}")
        if not args.md.is_file() or args.md.read_text(encoding="utf-8") != expected_md:
            raise SystemExit(f"postdiction ledger Markdown drift: {args.md}")
        print("postdiction ledger parity OK")
        return
    for name, rows in result["sections"].items():
        print(f"{name}: {len(rows)} rows")
    print(f"wrote {args.out}")
    print(f"wrote {args.md}")


if __name__ == "__main__":
    main()
