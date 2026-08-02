#!/usr/bin/env python3
"""Build the target-clean seam-current edge-ray candidate packet.

The finite source machine supplies thirty incidence seams.  The exact Lean
source proves that their antipodal-odd load readback has image ``D6`` and
that its response-selected metric completion is the same Euclidean
three-carrier obtained from the port Gram form.  This producer independently
replays the finite carrier geometry over ``Q(sqrt(5))``.  The thirty seam
differences form the edge-axis orbit and have the exact even moments

    sum (w.n)^2 = 10,
    sum (w.n)^4 = 6,
    sum (w.n)^6 = 30/7 - (2/7) I6(n).

If several additional physical premises turn those finite seam readbacks
into the sole homogeneous scalar propagation stencil, continuum
normalization gives the candidate symbol

    Λ_a(k,n) = (1/(5 a^2)) sum_j [1 - cos(a k w_j.n)]

and hence

    C4 = -a^2/20,  B0 = a^4/840,  B6 = -a^4/12600.

The exact source ray is not a closed physical producer.  In particular, the
finite incidence theorem does not prove that a seam is a spatial
translation, that the action is homogeneous, that this orbit is the sole
kinetic support, or that a laboratory field realizes the scalar symbol.
Every such premise is explicit below.  This program reads no target,
comparison, or public measurement data and cannot arm a comparison.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Any

import a5_multipole_fixed_point_certificate as base

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
RECEIPT_PATH = RUNTIME / "seam_current_edge_prediction_receipt.json"
SEAM_PROOF_PATH = REPO_ROOT / "Lean" / "Screen" / "SeamCurrentCarrierQuotient.lean"
ORBIT_PROOF_PATH = REPO_ROOT / "Lean" / "Screen" / "A5OrbitRaySeparation.lean"
SEAM_MOMENT_PROOF_PATH = (
    REPO_ROOT / "Lean" / "Screen" / "SeamCurrentEdge30Moment.lean"
)

SCHEMA = "oph.seam_current_edge_prediction_candidate.v1"
STATUS = (
    "EXACT_SOURCE_NATIVE_EDGE_RAY__PROSPECTIVE_PHYSICAL_BRANCH_UNARMED__"
    "PHYSICAL_PRODUCER_OPEN"
)

Vec3 = tuple[base.Q5, base.Q5, base.Q5]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise base.FingerprintError(message)


def load_typed_receipt(
    path: Path, *, schema: str, status: str
) -> tuple[bytes, dict[str, Any]]:
    """Load a canonical self-hashed receipt and refuse every typed drift."""

    raw = path.read_bytes()
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError(f"invalid parent JSON: {path.name}") from error
    require(raw == base.canonical_json_bytes(payload), f"noncanonical {path.name}")
    require(payload.get("schema") == schema, f"schema drift in {path.name}")
    require(payload.get("status") == status, f"status drift in {path.name}")
    body = {key: value for key, value in payload.items() if key != "receipt_sha256"}
    require(
        payload.get("receipt_sha256")
        == base.tagged_sha256(base.canonical_json_bytes(body)),
        f"self-digest drift in {path.name}",
    )
    return raw, payload


def load_exact_source(path: Path, required_fragments: tuple[str, ...]) -> bytes:
    """Pin a proof source and require its claim-bearing statements verbatim.

    Whitespace is normalized before checking theorem statements.  The full
    source bytes are pinned in the output, so any other edit is detected by
    :func:`verify_committed_receipt` as well.
    """

    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise base.FingerprintError(f"non-UTF-8 proof source: {path.name}") from error
    normalized = re.sub(r"\s+", " ", text)
    for fragment in required_fragments:
        require(fragment in normalized, f"claim-bearing source drift in {path.name}")
    require("sorry" not in text and "admit" not in text, f"placeholder in {path.name}")
    return raw


def q5_dot(left: Vec3, right: Vec3) -> base.Q5:
    total = base.ZERO
    for x, y in zip(left, right):
        total = base.q5_add(total, base.q5_mul(x, y))
    return total


def q5_vec_add(left: Vec3, right: Vec3) -> Vec3:
    return tuple(base.q5_add(x, y) for x, y in zip(left, right))  # type: ignore[return-value]


def q5_vec_sub(left: Vec3, right: Vec3) -> Vec3:
    return tuple(base.q5_sub(x, y) for x, y in zip(left, right))  # type: ignore[return-value]


def q5_vec_neg(vector: Vec3) -> Vec3:
    return tuple(base.q5_neg(x) for x in vector)  # type: ignore[return-value]


def sphere_equal(left: base.Poly, right: base.Poly) -> bool:
    difference = base.p_add(left, base.p_scale(right, -1))
    return base.p_is_zero(base.p_reduce_sphere(difference))


def exact_seam_orbit() -> dict[str, Any]:
    """Reconstruct the seam directions and their moments in exact arithmetic."""

    vertices = base.cartesian_vertices()
    inverse_norm = base.q5_div(base.ONE, base.NORM_SQ)

    def unit_dot(i: int, j: int) -> base.Q5:
        return base.q5_mul(q5_dot(vertices[i], vertices[j]), inverse_norm)

    seams = [
        (i, j)
        for i in range(12)
        for j in range(i + 1, 12)
        if unit_dot(i, j) == base.INV_SQRT5
    ]
    require(len(seams) == 30, "source seam census drift")
    require(
        all(sum(vertex in seam for seam in seams) == 5 for vertex in range(12)),
        "source seam degree drift",
    )

    differences = [q5_vec_sub(vertices[j], vertices[i]) for i, j in seams]
    midpoints = [q5_vec_add(vertices[i], vertices[j]) for i, j in seams]
    difference_norms = {q5_dot(vector, vector) for vector in differences}
    midpoint_norms = {q5_dot(vector, vector) for vector in midpoints}
    require(difference_norms == {base.q5(4)}, "seam-difference norm drift")
    require(
        midpoint_norms == {base.q5(6, 2)}, "edge-midpoint norm drift"
    )
    difference_norm = next(iter(difference_norms))
    midpoint_norm = next(iter(midpoint_norms))

    # Every oriented seam difference lies on one edge-midpoint axis.  Since
    # the midpoint list contains both signs, there are two matches per row.
    parallel_counts = []
    for difference in differences:
        matches = 0
        for midpoint in midpoints:
            dot_squared = base.q5_pow(q5_dot(difference, midpoint), 2)
            if dot_squared == base.q5_mul(difference_norm, midpoint_norm):
                matches += 1
        parallel_counts.append(matches)
    require(set(parallel_counts) == {2}, "seam-to-edge-axis binding drift")

    directed = differences + [q5_vec_neg(vector) for vector in differences]
    directed_multiplicity = Counter(directed)
    require(
        len(directed_multiplicity) == 30
        and set(directed_multiplicity.values()) == {2},
        "directed seam-to-edge multiplicity drift",
    )

    midpoint_m2 = base.normalized_moment(midpoints, 2, midpoint_norm)
    midpoint_m4 = base.normalized_moment(midpoints, 4, midpoint_norm)
    midpoint_m6 = base.normalized_moment(midpoints, 6, midpoint_norm)
    difference_m2 = base.normalized_moment(differences, 2, difference_norm)
    difference_m4 = base.normalized_moment(differences, 4, difference_norm)
    difference_m6 = base.normalized_moment(differences, 6, difference_norm)
    require(
        difference_m2 == midpoint_m2
        and difference_m4 == midpoint_m4
        and difference_m6 == midpoint_m6,
        "seam-difference and edge-midpoint moments diverged",
    )

    cartesian = base.build_cartesian_frame()
    i6 = cartesian["_i6_poly_object"]
    require(
        midpoint_m2 == base.p_scale(base.radial_power(1), 10),
        "edge second moment drift",
    )
    require(
        midpoint_m4 == base.p_scale(base.radial_power(2), 6),
        "edge fourth moment drift",
    )
    sixth_target = base.p_add(
        base.p_scale(base.radial_power(3), Fraction(30, 7)),
        base.p_scale(i6, Fraction(-2, 7)),
    )
    require(
        sphere_equal(midpoint_m6, sixth_target),
        "edge sixth-moment decomposition drift",
    )
    require(
        cartesian["edge_value"] == "-5/16 on all thirty edge midpoints",
        "normalized I6 edge value drift",
    )

    return {
        "source_ports": 12,
        "source_seams": len(seams),
        "port_degree": 5,
        "unoriented_axes": 15,
        "signed_edge_directions": len(directed_multiplicity),
        "directed_seam_labels": len(directed),
        "directed_labels_per_signed_direction": 2,
        "seam_difference_norm_squared": base.q5_str(difference_norm),
        "edge_midpoint_norm_squared": base.q5_str(midpoint_norm),
        "exact_axis_binding": (
            "every oriented incidence difference is parallel to exactly the "
            "two signed representatives of one edge-midpoint axis"
        ),
        "even_moments_on_unit_sphere": {
            "sum_w_dot_n_squared": "10",
            "sum_w_dot_n_fourth": "6",
            "sum_w_dot_n_sixth": "30/7 - (2/7) I6(n)",
        },
        "orientation_boundary": (
            "the finite source fixes incidence orientation only as bookkeeping; "
            "the cosine candidate is orientation independent"
        ),
    }


def build_receipt() -> dict[str, Any]:
    base_bytes, base_receipt = load_typed_receipt(
        base.RECEIPT_PATH, schema=base.SCHEMA, status=base.STATUS
    )
    require(
        base_receipt["decision_rules_and_ledger"]["comparison_boundary"][
            "public_measurement_read"
        ]
        is False,
        "target-free geometry parent became exposed",
    )
    require(
        base_receipt["cartesian_frame"]["edge_value"]
        == "-5/16 on all thirty edge midpoints",
        "edge-orbit parent drift",
    )

    seam_bytes = load_exact_source(
        SEAM_PROOF_PATH,
        (
            "theorem seam_table_complete :",
            "theorem seamAxisCurrent_eq_table (c : SeamCurrent) : seamAxisCurrent c = seamAxisTable c := by",
            "theorem exists_seamCurrent_iff_even (z : Axis → ℤ) : (\u2203 c : SeamCurrent, seamAxisCurrent c = z) ↔ EvenAxisTotal z := by",
            "noncomputable def d6CompletionEquivEuclidean3 : UniformSpace.Completion D6Point ≃ᵤ EuclideanVec3 :=",
            "it does not rename seams as axis translations",
        ),
    )
    orbit_bytes = load_exact_source(
        ORBIT_PROOF_PATH,
        (
            "| .edge30 => -1 / 12600",
            "theorem common_normalization : C4 = -1 / 20 ∧ B0 = 1 / 840 ∧ B0 / C4 ^ 2 = 10 / 21 := by",
            "theorem edge30_b6_over_c4_squared : b6OverC4Squared .edge30 = -2 / 63 := by",
            "theorem edge30_b6_over_b0 : b6OverB0 .edge30 = -1 / 15 := by",
            "No theorem here derives a row from OPH repair data",
        ),
    )
    seam_moment_bytes = load_exact_source(
        SEAM_MOMENT_PROOF_PATH,
        (
            "theorem carrierSeamDifference_eq_endpointDifference (e : Fin 30) : carrierSeamDifference e = endpointDifference e := by",
            "theorem projective_class_multisets_equal : (Finset.univ.val.map seamAxisClass : Multiset (Fin 15)) = Finset.univ.val.map edgeAxisClass := by",
            "theorem seamMoment2_eq (k : Vec3) : seamMoment2 k = 10 * radiusSquared k := by",
            "theorem seamMoment4_eq (k : Vec3) : seamMoment4 k = 6 * radiusSquared k ^ 2 := by",
            "theorem seamMoment6_eq (k : Vec3) : seamMoment6 k = (30 / 7 : ℝ) * radiusSquared k ^ 3 - (2 / 7 : ℝ) * I6 k := by",
            "theorem source_seam_edge30_control_certificate :",
            "No seam is declared to be a physical hop",
        ),
    )

    source_orbit = exact_seam_orbit()
    common_weight = Fraction(1, 5)
    c4 = -common_weight * Fraction(6, 24)
    b0 = common_weight * Fraction(30, 7 * 720)
    b6 = common_weight * Fraction(-2, 7 * 720)
    require(c4 == Fraction(-1, 20), "C4 derivation drift")
    require(b0 == Fraction(1, 840), "B0 derivation drift")
    require(b6 == Fraction(-1, 12600), "B6 derivation drift")

    b0_over_c4_squared = b0 / (c4 * c4)
    b6_over_c4_squared = b6 / (c4 * c4)
    b6_over_b0 = b6 / b0
    require(b0_over_c4_squared == Fraction(10, 21), "B0/C4^2 drift")
    require(b6_over_c4_squared == Fraction(-2, 63), "B6/C4^2 drift")
    require(b6_over_b0 == Fraction(-1, 15), "B6/B0 drift")

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 666,
        "producer_scope": {
            "name": "source-native seam-current edge-ray candidate",
            "type": "target-clean prospective conditional physical-branch candidate",
            "physical_producer_closed": False,
            "frozen_prediction_registered": False,
            "statement": (
                "the source incidence/readback chain fixes the finite edge ray; "
                "a nature-facing prediction exists only if every physical "
                "promotion gate in this receipt is discharged"
            ),
        },
        "parent_pins": [
            {
                "path": "code/a5_fingerprint/runtime/a5_multipole_fixed_point_receipt.json",
                "bytes": len(base_bytes),
                "sha256": base.tagged_sha256(base_bytes),
                "role": "exact target-free carrier geometry and normalized I6",
            },
            {
                "path": "Lean/Screen/SeamCurrentCarrierQuotient.lean",
                "bytes": len(seam_bytes),
                "sha256": base.tagged_sha256(seam_bytes),
                "role": "exact thirty-seam source image D6 and carrier completion",
            },
            {
                "path": "Lean/Screen/A5OrbitRaySeparation.lean",
                "bytes": len(orbit_bytes),
                "sha256": base.tagged_sha256(orbit_bytes),
                "role": "kernel-checked edge-ray arithmetic and orbit separation",
            },
            {
                "path": "Lean/Screen/SeamCurrentEdge30Moment.lean",
                "bytes": len(seam_moment_bytes),
                "sha256": base.tagged_sha256(seam_moment_bytes),
                "role": (
                    "kernel-checked source-seam to edge-orbit binding, exact "
                    "moments, and coefficient derivation"
                ),
            },
        ],
        "exact_source_result": {
            "seam_current_image": (
                "D6 = {z in Z^6 : the coordinate sum is even}; the residual "
                "cokernel is one parity bit"
            ),
            "metric_completion": (
                "the D6 readback is dense in, and completes to, the same "
                "response-selected Euclidean three-carrier"
            ),
            "finite_geometry_replay": source_orbit,
            "scope_boundary": (
                "these are exact finite incidence, quotient, orbit, and metric "
                "facts; they do not identify seam currents with physical motion"
            ),
        },
        "conditional_physical_candidate": {
            "symbol_name": "Λ (a spatial kinetic eigenvalue, not automatically a frequency squared)",
            "domain": "k >= 0 and n is a unit direction on S^2",
            "operator": (
                "Λ_a(k,n) = (1/(5 a^2)) sum_{j=1}^{30} "
                "[1 - cos(a k w_j.n)]"
            ),
            "template": (
                "I6 is the unique normalized A5-invariant rank-six harmonic "
                "with I6=1 on the twelve primitive port directions"
            ),
            "expansion": (
                "Λ_a = k^2 - (a^2/20) k^4 + (a^4/840) k^6 - "
                "(a^4/12600) k^6 I6(n) + O(a^6 k^8)"
            ),
            "coefficients": {
                "C4_over_a2": str(c4),
                "B0_over_a4": str(b0),
                "B6_over_a4": str(b6),
            },
            "scale_free_relations": {
                "B0_over_C4_squared": str(b0_over_c4_squared),
                "B6_over_C4_squared": str(b6_over_c4_squared),
                "B6_over_B0": str(b6_over_b0),
            },
            "signs": {"C4": "negative", "B0": "positive", "B6": "negative"},
            "harmonic_nulls": (
                "intrinsic anisotropic ranks one through five vanish on the "
                "declared complete equal-weight edge orbit"
            ),
            "fit_freedom_after_C4": (
                "one shared orientation in SO(3)/A5; no scale, amplitude, or "
                "independent rank-six shape coefficient remains"
            ),
        },
        "physical_premises": {
            "seam_as_displacement": (
                "the algebraic boundary/readback of each source seam is a "
                "physical spatial displacement in the reconstructed carrier"
            ),
            "homogeneous_translation_action": (
                "the finite displacement extends at every carrier point as one "
                "translation action with the same local coefficients"
            ),
            "complete_support": (
                "the complete thirty-direction edge orbit is the sole intrinsic "
                "directional support through the displayed derivative order"
            ),
            "equal_weights": (
                "proper-carrier covariance and transitivity select one common "
                "coefficient, with no independent invariant orbit channel"
            ),
            "continuum_normalization": (
                "the quadratic term is k^2, fixing the common coefficient to "
                "1/(5 a^2)"
            ),
            "finite_scale": "a is finite, strictly positive, and common to the tested domain",
            "positive_scale_lower_bound": (
                "a source-derived physical lower bound a >= a_min > 0 is "
                "available in the tested sector and units; mathematical a>0 "
                "alone supplies no experimental power against a null"
            ),
            "physical_sector": (
                "the tested field realizes this scalar Λ symbol; a photon "
                "reading also requires identical action on both transverse polarizations"
            ),
            "cofinal_gluing": (
                "the local finite actions glue consistently across scale and "
                "observer charts without changing the edge ray"
            ),
            "frame_and_boost": (
                "one carrier frame and its SO(3)/A5 orientation have a declared "
                "transport and boost map into the comparison frame"
            ),
            "readout_and_nuisance": (
                "source, medium, gravity, instrument, polarization, boost, and "
                "orientation effects are identified or profiled under a frozen model"
            ),
            "exclusivity": (
                "the fitted intrinsic coefficients isolate this carrier term "
                "from every other allowed contribution at the displayed orders"
            ),
        },
        "promotion_gates": {
            "all_discharged": False,
            "physical_producer_closed": False,
            "comparison_eligible": False,
            "gates": [
                {"gate": "seam-as-displacement identification", "status": "OPEN", "owner": "#666"},
                {"gate": "homogeneous translation action", "status": "OPEN", "owner": "#663/#666"},
                {"gate": "sole complete edge support and equal weights", "status": "OPEN", "owner": "#655/#666"},
                {"gate": "physical scalar or polarization-independent sector", "status": "OPEN", "owner": "#655/#666"},
                {"gate": "cofinal scale and observer gluing", "status": "OPEN", "owner": "#663"},
                {"gate": "finite physical carrier scale", "status": "OPEN", "owner": "#664"},
                {"gate": "source-derived positive physical scale lower bound", "status": "OPEN", "owner": "#664"},
                {"gate": "frame, boost, readout, nuisance, and exclusivity contract", "status": "OPEN", "owner": "#666"},
                {"gate": "dataset-specific post-custody preregistration", "status": "OPEN", "owner": "#639"},
            ],
            "fail_closed_rule": (
                "no comparison may be armed and no OPH-wide falsification claim "
                "may be made while any gate is open"
            ),
        },
        "baseline_contrast": {
            "baseline": (
                "minimal locally Lorentz-invariant Standard Model plus General "
                "Relativity in local vacuum"
            ),
            "baseline_prediction": "C4 = B0 = B6 = 0 for intrinsic vacuum propagation",
            "nonuniqueness": (
                "nonminimal effective operators or another icosahedral medium "
                "can imitate the edge relation; support would distinguish this "
                "branch from the baseline without identifying OPH uniquely"
            ),
        },
        "prospective_decision_rule": {
            "eligibility": (
                "all physical promotion gates are discharged before exposure, "
                "and a later data release supplies a joint likelihood or full "
                "covariance for same-sector C4, isotropic B0, and the complete "
                "rank-six vector"
            ),
            "trigger": (
                "C4 is measured negative at at least five standard deviations "
                "with calibrated sensitivity to both linked sixth-order terms"
            ),
            "fail": (
                "after profiling the predeclared SO(3)/A5 orientation, the linked "
                "B0 and negative-B6 edge relation is excluded at at least five "
                "standard deviations with calibrated joint coverage"
            ),
            "support": (
                "the zero-coefficient baseline is excluded at at least five "
                "standard deviations, the edge manifold agrees within two "
                "standard deviations, named systematics are rejected, and an "
                "independent release replicates the result"
            ),
            "inconclusive": (
                "C4 is not detected, the sixth-order terms are below sensitivity, "
                "the covariance is incomplete, the carrier contribution cannot "
                "be isolated, or issue #664 supplies no source-derived positive "
                "physical lower bound on a"
            ),
            "no_null_verdict": (
                "a null result cannot reject the branch unless issue #664 first "
                "provides a source-derived lower bound that places a nonzero "
                "signal inside the calibrated sensitivity domain; the premise "
                "a>0 by itself is insufficient"
            ),
            "scope_of_failure": (
                "failure rejects this seam-current edge propagation branch; it "
                "rejects OPH as a whole only after a separate theorem proves the "
                "branch forced and exclusive"
            ),
        },
        "fz11_separation": {
            "fz11_register_id": "FZ-11",
            "fz11_prediction_receipt_read": False,
            "fz11_bytes_modified": False,
            "supersedes_fz11": False,
            "relationship": (
                "FZ-11 freezes the conditional primitive-vertex ray; this packet "
                "records the distinct source-native edge ray and cannot amend, "
                "reinterpret, or score FZ-11"
            ),
            "vertex_B6_over_C4_squared": "32/315",
            "edge_B6_over_C4_squared": str(b6_over_c4_squared),
            "opposite_rank_six_sign": True,
        },
        "exposure_and_custody_boundary": {
            "target_values_read": False,
            "comparison_data_read": False,
            "public_measurement_read": False,
            "comparison_permitted": False,
            "comparison_state": "INELIGIBLE_UNARMED_PHYSICAL_PRODUCER_OPEN",
            "comparison_inputs": [],
            "source_inputs_only": [
                "exact local carrier receipt",
                "exact local Lean seam-current source",
                "exact local Lean orbit-ray source",
                "exact local Lean source-seam moment source",
            ],
            "excluded_data_class": (
                "the 2026-07-17 WMAP campaign, its CMB likelihood class, every "
                "data product inspected for FZ-11, and all data seen before a "
                "separate edge-branch custody freeze are ineligible"
            ),
            "candidate_registration": (
                "this runtime receipt is not the frozen prediction register or "
                "a dataset-specific preregistration"
            ),
        },
    }
    receipt["receipt_sha256"] = base.tagged_sha256(base.canonical_json_bytes(receipt))
    return receipt


def verify_committed_receipt() -> dict[str, Any]:
    """Refuse parent, source, payload, or self-digest drift."""

    raw = RECEIPT_PATH.read_bytes()
    try:
        committed = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError("invalid committed edge receipt") from error
    require(raw == base.canonical_json_bytes(committed), "noncanonical edge receipt")
    require(committed.get("schema") == SCHEMA, "edge receipt schema drift")
    require(committed.get("status") == STATUS, "edge receipt status drift")
    body = {key: value for key, value in committed.items() if key != "receipt_sha256"}
    require(
        committed.get("receipt_sha256")
        == base.tagged_sha256(base.canonical_json_bytes(body)),
        "edge receipt self-digest drift",
    )
    rebuilt = build_receipt()
    require(raw == base.canonical_json_bytes(rebuilt), "edge receipt parent/source drift")
    return committed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)
    if args.write and args.verify:
        parser.error("choose --write or --verify")
    if args.verify:
        receipt = verify_committed_receipt()
    else:
        receipt = build_receipt()
    if args.write:
        RUNTIME.mkdir(exist_ok=True)
        RECEIPT_PATH.write_bytes(base.canonical_json_bytes(receipt))
        print(RECEIPT_PATH)
        return 0
    print(json.dumps(receipt["conditional_physical_candidate"], indent=1))
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
