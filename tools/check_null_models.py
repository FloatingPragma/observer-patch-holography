#!/usr/bin/env python3
"""Generate and check the OPH null-model scorecard.

The two checks are deliberately negative controls:

* substitute the constants in the declared gauge-width P-closure map and
  interval-certify every resulting fixed point;
* expose the declared deltahedral/port menu and distinguish "available" from
  "uniquely selected".

This tool never promotes a physical claim.  Missing alternative producers and
undeclared selector menus are reported as unavailable, not silently counted as
failures of the alternatives.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCORECARD = Path("code/P_derivation/null_model_scorecard.md")
TARGET_ALPHA_INV = Decimal("137.035999177")
THRESHOLDS = (
    Decimal("1e-4"),
    Decimal("1e-5"),
    Decimal("2.5e-6"),
)
MODE = "thomson_structured_running_plus_gauge_width"
POINT_DPS = 40
INTERVAL_DPS = 40
SU2_CUTOFF = 40
SU3_CUTOFF = 30
ALPHA_INTERVAL_HALF_WIDTH = "1e-13"
NONTRIVIAL_MIN_COUNT = 2
NONTRIVIAL_MIN_FRACTION = Decimal("0.05")
class NullModelError(RuntimeError):
    """Raised when a null-model source or gate fails closed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise NullModelError(message)


@dataclass(frozen=True)
class ConstantChoice:
    key: str
    expression: str
    equivalence_key: str


C1_CHOICES = (
    ConstantChoice("phi", "φ", "phi"),
    ConstantChoice("1.60", "1.60", "8/5"),
    ConstantChoice("1.62", "1.62", "81/50"),
    ConstantChoice("1.65", "1.65", "33/20"),
    ConstantChoice("e_over_1_7", "e/1.7", "10e/17"),
    ConstantChoice("eight_fifths", "8/5", "8/5"),
    ConstantChoice("sqrt_e", "√e", "sqrt(e)"),
    ConstantChoice("three_minus_sqrt_two", "3−√2", "3-sqrt(2)"),
)
C2_CHOICES = (
    ConstantChoice("pi", "π", "pi"),
    ConstantChoice("one", "1", "1"),
    ConstantChoice("two", "2", "2"),
    ConstantChoice("e", "e", "e"),
    ConstantChoice("three", "3", "3"),
    ConstantChoice("pi_over_two", "π/2", "pi/2"),
)

DELTAHEDRAL_MENU = (
    ("tetrahedron", 4),
    ("triangular_bipyramid", 5),
    ("octahedron", 6),
    ("pentagonal_bipyramid", 7),
    ("snub_disphenoid", 8),
    ("triaugmented_triangular_prism", 9),
    ("gyroelongated_square_bipyramid", 10),
    ("icosahedron", 12),
)

INPUT_PATHS = (
    Path("code/P_derivation/interval_contraction_certificate.py"),
    Path("code/a5_closure/echosahedral_selector_certificate.py"),
    Path("code/a5_closure/manifests/echosahedral_federation_reference.json"),
    Path("code/a5_closure/port_current_inner_certificate.py"),
    Path("code/a5_closure/manifests/port_current_response_reference.json"),
    Path("code/a5_closure/a5_screen_sm_closure.json"),
    Path("code/a5_closure/receipts/echosahedral_federation_reference.receipt.json"),
    Path("code/a5_closure/receipts/port_current_inner_reference.receipt.json"),
    Path("code/a5_closure/receipts/super_tannakian_matter_reference.receipt.json"),
    Path("code/a5_closure/claim_boundary_certificates.py"),
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise NullModelError(f"cannot load JSON input {path}: {exc}") from exc
    require(isinstance(value, dict), f"JSON input must be an object: {path}")
    return value


def _resolve_root(root: Path) -> Path:
    resolved = root.resolve()
    require(resolved.is_dir(), f"repository root does not exist: {resolved}")
    missing = [str(path) for path in INPUT_PATHS if not (resolved / path).is_file()]
    require(not missing, f"null-model inputs are missing under {resolved}: {missing}")
    return resolved


def declared_constant_grid() -> tuple[tuple[ConstantChoice, ConstantChoice], ...]:
    """Return the ordered, explicitly declared substitution grid."""

    return tuple((c1, c2) for c1 in C1_CHOICES for c2 in C2_CHOICES)


def _load_module(path: Path, name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _constant_value(backend: Any, key: str) -> Any:
    """Evaluate one declared expression over either point or interval backend."""

    if key == "phi":
        return (backend.one + backend.sqrt(backend.num(5))) / backend.two
    if key == "1.60" or key == "eight_fifths":
        return backend.num(8) / backend.num(5)
    if key == "1.62":
        return backend.num(81) / backend.num(50)
    if key == "1.65":
        return backend.num(33) / backend.num(20)
    if key == "e_over_1_7":
        return backend.exp(backend.one) * backend.num(10) / backend.num(17)
    if key == "sqrt_e":
        return backend.sqrt(backend.exp(backend.one))
    if key == "three_minus_sqrt_two":
        return backend.num(3) - backend.sqrt(backend.num(2))
    if key == "pi":
        return backend.pi
    if key == "one":
        return backend.one
    if key == "two":
        return backend.two
    if key == "e":
        return backend.exp(backend.one)
    if key == "three":
        return backend.num(3)
    if key == "pi_over_two":
        return backend.pi / backend.two
    raise NullModelError(f"unknown declared constant key: {key}")


def _decimal(text: str) -> Decimal:
    try:
        return Decimal(text)
    except Exception as exc:  # pragma: no cover - defensive parse boundary
        raise NullModelError(f"invalid decimal emitted by interval backend: {text}") from exc


def _distance_bounds(lo: Decimal, hi: Decimal) -> tuple[Decimal, Decimal]:
    require(lo <= hi, "root enclosure is reversed")
    if TARGET_ALPHA_INV < lo:
        return lo - TARGET_ALPHA_INV, hi - TARGET_ALPHA_INV
    if TARGET_ALPHA_INV > hi:
        return TARGET_ALPHA_INV - hi, TARGET_ALPHA_INV - lo
    return Decimal(0), max(TARGET_ALPHA_INV - lo, hi - TARGET_ALPHA_INV)


def _certify_constant_pair(
    interval_module: Any,
    c1: ConstantChoice,
    c2: ConstantChoice,
) -> dict[str, Any]:
    """Certify one substituted fixed point with the existing interval chain."""

    icc = interval_module
    mpb = icc.MpBackend(POINT_DPS)
    mpb.phi = _constant_value(mpb, c1.key)
    mpb.sqrt_pi = mpb.sqrt(_constant_value(mpb, c2.key))
    point = icc.PointSolver(mpb, SU2_CUTOFF, SU3_CUTOFF)
    alpha_star = point.fixed_point(icc.MODE_GAUGE_WIDTH)

    ivb = icc.IvBackend(INTERVAL_DPS)
    ivb.phi = _constant_value(ivb, c1.key)
    ivb.sqrt_pi = ivb.sqrt(_constant_value(ivb, c2.key))
    chain = icc.IntervalChain(ivb, point, SU2_CUTOFF, SU3_CUTOFF)
    # DualBackend initializes the canonical constants independently.  Replace
    # those two constant duals with the declared substitution.
    chain.dual.phi = icc.Dual(ivb.phi, ivb.zero)
    chain.dual.sqrt_pi = icc.Dual(ivb.sqrt_pi, ivb.zero)

    # These are bracket-inflation hints, not assumed derivative bounds.  Every
    # inner root is subsequently accepted only after interval endpoint signs
    # and a sign-definite residual derivative are verified by IntervalChain.
    floor = mpb.num("1e-34")
    chain.scales = {
        "dudp_abs": mpb.num(1),
        "dmz_rel_dp": mpb.num(10),
        "dmz_rel_du": mpb.num(100),
        "u_floor": floor,
        "mz_floor": floor,
    }

    radius = mpb.num(ALPHA_INTERVAL_HALF_WIDTH)
    interval = ivb.hull(
        ivb.thin(alpha_star - radius).a,
        ivb.thin(alpha_star + radius).b,
    )
    g_box, diagnostics = chain.g_dual(interval, icc.MODE_GAUGE_WIDTH)
    lipschitz = ivb.sup_abs(g_box.d)
    require(lipschitz < ivb.one.a, f"{c1.expression}, {c2.expression}: not a contraction")

    midpoint = ivb.thin(alpha_star)
    g_mid, midpoint_diagnostics = chain.g_dual(midpoint, icc.MODE_GAUGE_WIDTH)
    g_enclosure = g_mid.x + g_box.d * (interval - midpoint)
    self_map = interval.a < g_enclosure.a and g_enclosure.b < interval.b
    require(self_map, f"{c1.expression}, {c2.expression}: centered form is not a self-map")

    inverse_enclosure = ivb.one / g_enclosure
    inverse_pair = icc._iv_pair(inverse_enclosure)
    root_lo = _decimal(inverse_pair["lo"])
    root_hi = _decimal(inverse_pair["hi"])
    distance_lo, distance_hi = _distance_bounds(root_lo, root_hi)
    with localcontext() as context:
        context.prec = 50
        relative_lo = +(distance_lo / TARGET_ALPHA_INV)
        relative_hi = +(distance_hi / TARGET_ALPHA_INV)
        point_alpha_inv = +(Decimal(mpb.ctx.nstr(mpb.one / alpha_star, 38)))

    threshold_results = {
        format(threshold, "E").replace("E+", "e").replace("E-", "e-"): {
            "certified_inside": relative_hi <= threshold,
            "certified_outside": relative_lo > threshold,
        }
        for threshold in THRESHOLDS
    }
    require(
        all(row["certified_inside"] != row["certified_outside"] for row in threshold_results.values()),
        f"{c1.expression}, {c2.expression}: a threshold intersects the root enclosure",
    )
    alpha_u_checks = diagnostics["alpha_u_verification"]
    mz_checks = diagnostics["mz_verification"]
    require(
        alpha_u_checks["endpoint_signs_verified"]
        and alpha_u_checks["R_u_sign_definite"]
        and mz_checks["endpoint_signs_verified"]
        and mz_checks["h_m_sign_definite"],
        f"{c1.expression}, {c2.expression}: an inner-root certificate failed",
    )
    midpoint_alpha_u = midpoint_diagnostics["alpha_u_verification"]
    midpoint_mz = midpoint_diagnostics["mz_verification"]
    require(
        midpoint_alpha_u["endpoint_signs_verified"]
        and midpoint_alpha_u["R_u_sign_definite"]
        and midpoint_mz["endpoint_signs_verified"]
        and midpoint_mz["h_m_sign_definite"],
        f"{c1.expression}, {c2.expression}: midpoint inner-root certificate failed",
    )

    return {
        "c1_key": c1.key,
        "c1": c1.expression,
        "c1_equivalence_key": c1.equivalence_key,
        "c2_key": c2.key,
        "c2": c2.expression,
        "c2_equivalence_key": c2.equivalence_key,
        "canonical_pair": c1.key == "phi" and c2.key == "pi",
        "numeric_pair_key": f"{c1.equivalence_key}|{c2.equivalence_key}",
        "alpha_inv_point": format(point_alpha_inv, "f"),
        "alpha_inv_enclosure": {
            "lo": format(root_lo, "f"),
            "hi": format(root_hi, "f"),
        },
        "absolute_distance_enclosure": {
            "lo": format(distance_lo, "f"),
            "hi": format(distance_hi, "f"),
        },
        "relative_distance_enclosure": {
            "lo": format(relative_lo, "E").replace("E+", "e").replace("E-", "e-"),
            "hi": format(relative_hi, "E").replace("E+", "e").replace("E-", "e-"),
        },
        "threshold_results": threshold_results,
        "certificate": {
            "status": "INTERVAL_BANACH_ROOT_CERTIFIED",
            "mode": MODE,
            "self_map": True,
            "contraction": True,
            "lipschitz_upper": icc._iv_upper_str(lipschitz),
            "alpha_interval_half_width": ALPHA_INTERVAL_HALF_WIDTH,
            "pixel_inner_root_endpoint_signs_verified": True,
            "pixel_inner_root_derivative_sign_definite": True,
            "mz_inner_root_endpoint_signs_verified": True,
            "mz_inner_root_derivative_sign_definite": True,
            "edge_sum_tails_included": True,
        },
    }


def _deduplicated_alternatives(rows: Sequence[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    seen: set[str] = set()
    output: list[Mapping[str, Any]] = []
    for row in rows:
        key = str(row["numeric_pair_key"])
        if key in seen:
            continue
        seen.add(key)
        if not row["canonical_pair"]:
            output.append(row)
    return output


def summarize_constant_scan(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    require(len(rows) == len(C1_CHOICES) * len(C2_CHOICES), "constant grid is incomplete")
    canonical = [row for row in rows if row["canonical_pair"]]
    require(len(canonical) == 1, "constant grid must contain exactly one canonical pair")
    alternatives = _deduplicated_alternatives(rows)
    unique_keys = {str(row["numeric_pair_key"]) for row in rows}

    thresholds: dict[str, dict[str, Any]] = {}
    for threshold in THRESHOLDS:
        key = format(threshold, "E").replace("E+", "e").replace("E-", "e-")
        inside = sum(bool(row["threshold_results"][key]["certified_inside"]) for row in alternatives)
        with localcontext() as context:
            context.prec = 30
            fraction = Decimal(inside) / Decimal(len(alternatives))
        thresholds[key] = {
            "certified_alternative_hits": inside,
            "unique_alternative_pairs": len(alternatives),
            "fraction": format(fraction, "f"),
        }

    decisive_key = "2.5e-6"
    decisive = thresholds[decisive_key]
    decisive_count = int(decisive["certified_alternative_hits"])
    decisive_fraction = Decimal(str(decisive["fraction"]))
    trigger = (
        decisive_count >= NONTRIVIAL_MIN_COUNT
        and decisive_fraction >= NONTRIVIAL_MIN_FRACTION
    )
    return {
        "declared_pair_count": len(rows),
        "numerically_distinct_pair_count": len(unique_keys),
        "unique_alternative_pair_count": len(alternatives),
        "duplicate_alias_pairs": len(rows) - len(unique_keys),
        "threshold_basis": "relative absolute distance |root-target|/target",
        "thresholds": thresholds,
        "interpretation_rule": {
            "threshold": decisive_key,
            "nontrivial_minimum_count": NONTRIVIAL_MIN_COUNT,
            "nontrivial_minimum_fraction": format(NONTRIVIAL_MIN_FRACTION, "f"),
            "triggered": trigger,
            "if_triggered": (
                "The phi/sqrt(pi) detuning story carries no evidential weight "
                "because a non-trivial fraction of declared substitutes reaches "
                "the OPH-scale residual."
            ),
            "current_interpretation": (
                "NO_EVIDENTIAL_WEIGHT: non-trivial null substitutions reach the "
                "OPH-scale residual."
                if trigger
                else
                "NO_POSITIVE_WEIGHT_FROM_THIS_SCAN: the trigger is not met, but "
                "the audit-declared grid is neither a source-derived nor an "
                "exhaustive/pre-registered alternative menu, so the converse "
                "does not establish the phi/sqrt(pi) story."
            ),
        },
        "canonical": canonical[0],
    }


def build_constant_scan(root: Path) -> dict[str, Any]:
    path = root / "code/P_derivation/interval_contraction_certificate.py"
    module = _load_module(path, f"_oph_interval_null_models_{hash(path)}")
    rows = [
        _certify_constant_pair(module, c1, c2)
        for c1, c2 in declared_constant_grid()
    ]
    return {
        "status": "PASS_INTERVAL_CERTIFICATES_COMPLETE",
        "claim_boundary": (
            "This is a sensitivity/null-model audit of the declared incomplete "
            "gauge-width numerical map. It is not a physical Thomson-endpoint "
            "derivation and does not promote any alpha claim."
        ),
        "target_alpha_inv_compare_only": format(TARGET_ALPHA_INV, "f"),
        "grid_status": (
            "audit_declared_probe_menu_not_source_derived_not_exhaustive_"
            "not_a_pre_registered_hypothesis_space"
        ),
        "c1_choices": [
            {
                "key": row.key,
                "expression": row.expression,
                "equivalence_key": row.equivalence_key,
            }
            for row in C1_CHOICES
        ],
        "c2_choices": [
            {
                "key": row.key,
                "expression": row.expression,
                "equivalence_key": row.equivalence_key,
            }
            for row in C2_CHOICES
        ],
        "solver": {
            "source": "code/P_derivation/interval_contraction_certificate.py",
            "mode": MODE,
            "point_dps": POINT_DPS,
            "interval_dps": INTERVAL_DPS,
            "su2_cutoff": SU2_CUTOFF,
            "su3_cutoff": SU3_CUTOFF,
            "edge_sum_tail_bounds_extend_to_infinite_sums": True,
        },
        "rows": rows,
        "summary": summarize_constant_scan(rows),
    }


@dataclass(frozen=True)
class _SourceReduct:
    port_count: int
    total_coordination_charge: int


def _execute_small_nonuniqueness_producer(path: Path) -> dict[str, Any]:
    """Execute only the dependency-free small no-go function from its source."""

    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    matches = [
        node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == "source_completion_nonuniqueness_certificate"
    ]
    require(len(matches) == 1, "source completion non-uniqueness producer is missing or duplicated")
    function = matches[0]
    module = ast.Module(body=[function], type_ignores=[])

    def local_require(condition: bool, message: str) -> None:
        if not condition:
            raise NullModelError(message)

    namespace: dict[str, Any] = {
        "Any": Any,
        "SourceReduct": _SourceReduct,
        "require": local_require,
    }
    exec(compile(module, str(path), "exec"), namespace)
    result = namespace["source_completion_nonuniqueness_certificate"]()
    require(isinstance(result, dict), "source completion no-go producer returned no object")
    return result


def _run_exact_verifier(
    root: Path,
    script: Path,
    manifest: Path,
    receipt: Path,
) -> None:
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    command = [
        sys.executable,
        str(root / script),
        "verify",
        "--manifest",
        str(root / manifest),
        "--receipt",
        str(root / receipt),
    ]
    completed = subprocess.run(
        command,
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    require(
        completed.returncode == 0,
        (
            f"exact selector/current verifier failed for {script}: "
            f"{completed.stderr.strip() or completed.stdout.strip()}"
        ),
    )


def build_selector_ablation(root: Path) -> dict[str, Any]:
    screen_path = root / "code/a5_closure/a5_screen_sm_closure.json"
    selector_path = (
        root
        / "code/a5_closure/receipts/echosahedral_federation_reference.receipt.json"
    )
    current_path = (
        root / "code/a5_closure/receipts/port_current_inner_reference.receipt.json"
    )
    matter_path = (
        root / "code/a5_closure/receipts/super_tannakian_matter_reference.receipt.json"
    )
    survival_path = root / "code/a5_closure/claim_boundary_certificates.py"
    screen = _load_json(screen_path)
    selector = _load_json(selector_path)
    current = _load_json(current_path)
    matter = _load_json(matter_path)
    no_go = _execute_small_nonuniqueness_producer(survival_path)
    _run_exact_verifier(
        root,
        Path("code/a5_closure/echosahedral_selector_certificate.py"),
        Path("code/a5_closure/manifests/echosahedral_federation_reference.json"),
        Path("code/a5_closure/receipts/echosahedral_federation_reference.receipt.json"),
    )
    _run_exact_verifier(
        root,
        Path("code/a5_closure/port_current_inner_certificate.py"),
        Path("code/a5_closure/manifests/port_current_response_reference.json"),
        Path("code/a5_closure/receipts/port_current_inner_reference.receipt.json"),
    )

    vertex = screen.get("vertex_module", {})
    require(vertex.get("exact_match") is True, "12-port A5/SM adjoint character match drifted")
    require(
        vertex.get("decomposition") == {"1": 1, "3": 1, "3prime": 1, "5": 1},
        "12-port A5 vertex decomposition drifted",
    )
    require(
        selector.get("schema") == "oph.echosahedral_selector_receipt.v1",
        "echosahedral selector receipt schema drifted",
    )
    require(
        selector.get("carrier", {}).get("ports") == 12,
        "echosahedral selector no longer emits twelve ports",
    )
    selector_result = selector.get("icosahedral_selector", {})
    require(
        selector_result.get("orientation_preserving_order") == 60
        and selector_result.get("port_action_faithful") is True
        and selector_result.get("port_action_transitive") is True,
        "echosahedral A5 selector certificate drifted",
    )
    trichotomy = screen.get("compact_lie_trichotomy")
    require(isinstance(trichotomy, list), "compact Lie trichotomy is missing")
    lie_types = [row.get("lie_type") for row in trichotomy if row.get("viable") is True]
    expected_lie_types = [
        "su(3)+su(2)+u(1)",
        "su(2)+su(2)+u(1)^6",
        "u(1)^12",
    ]
    require(lie_types == expected_lie_types, "12-port compact Lie alternative menu drifted")
    current_gate = current.get("issue_closure_condition", {})
    require(
        current_gate.get("conditional_algebraic_gate_passed") is True,
        "12-port conditional current gate no longer passes",
    )
    require(
        current_gate.get("physical_source_realization_gate_passed") is False
        and current_gate.get("met_locally") is False,
        "12-port current receipt must fail closed at source realization",
    )
    require(
        isinstance(current.get("semantic_response_binding"), dict)
        and current["semantic_response_binding"].get("sector_structure_recomputed")
        is True,
        "12-port response constraints are not recomputed",
    )
    require(
        current["semantic_response_binding"].get("current_lift_source_selected")
        is False
        and current["semantic_response_binding"].get(
            "overlap_holonomy_internality_certified"
        )
        is False,
        "12-port response artifact crossed the current-source boundary",
    )
    matter_selection = matter.get("selection", {})
    require(matter_selection.get("projector_rank") == 15, "matter projector rank drifted")
    require(
        matter.get("candidate_matter_class", {}).get("uniqueness_promoted")
        is False,
        "matter receipt unexpectedly promotes candidate-class uniqueness",
    )
    does_not_close = matter.get("claim_boundary", {}).get("does_not_close", [])
    require(
        any("exclusion of other anomaly-free light sectors" in str(item) for item in does_not_close),
        "matter receipt no longer exposes the no-extra-sector gate",
    )
    require(
        no_go.get("status") == "EXACT_FINITE_NONIDENTIFIABILITY_THEOREM",
        "source completion non-identifiability status drifted",
    )
    require(
        no_go.get("inequivalent_current_completions")
        == ["abelian_12", "compact_sm_lie_type"],
        "current completion countermodels drifted",
    )
    require(
        no_go.get("inequivalent_matter_completions")
        == ["rank15_exterior_packet", "rank15_plus_sterile_singlet"],
        "matter completion countermodels drifted",
    )
    require(
        no_go.get("source_only_reconstruction_of_every_completion") is False,
        "source reduct unexpectedly became completion-unique",
    )

    names = [name for name, _ in DELTAHEDRAL_MENU]
    counts = [count for _, count in DELTAHEDRAL_MENU]
    require(len(names) == len(set(names)) == 8, "deltahedral menu names are not unique")
    require(counts == [4, 5, 6, 7, 8, 9, 10, 12], "deltahedral port-count menu drifted")

    rows: list[dict[str, Any]] = []
    for name, count in DELTAHEDRAL_MENU:
        if count != 12:
            rows.append(
                {
                    "configuration": name,
                    "port_count": count,
                    "source_selector_producer": "UNDECLARED",
                    "port_response_producer": "UNDECLARED",
                    "compatible_compact_lie_types": None,
                    "sm_lie_type_available": "UNKNOWN_FAIL_CLOSED",
                    "sm_lie_type_uniquely_selected": False,
                    "status": "UNDECLARED_ALTERNATIVE_PRODUCER",
                }
            )
            continue
        rows.append(
            {
                "configuration": name,
                "port_count": count,
                "source_selector_producer": (
                    "code/a5_closure/receipts/"
                    "echosahedral_federation_reference.receipt.json"
                ),
                "port_response_producer": (
                    "code/a5_closure/receipts/"
                    "port_current_inner_reference.receipt.json"
                ),
                "compatible_compact_lie_types": lie_types,
                "sm_lie_type_available": True,
                "sm_lie_type_uniquely_selected": False,
                "conditional_response_gate_passed": True,
                "physical_source_binding": False,
                "status": "AXIOM_DOMAIN_FORCED__COUNTERFACTUAL_MENU_UNTESTED",
            }
        )

    computed = [row for row in rows if row["compatible_compact_lie_types"] is not None]
    undeclared = [row for row in rows if row["compatible_compact_lie_types"] is None]
    unique_selected = [row for row in rows if row["sm_lie_type_uniquely_selected"] is True]
    available = [row for row in rows if row["sm_lie_type_available"] is True]
    return {
        "status": "FAIL_CLOSED_SELECTOR_MENU_INCOMPLETE",
        "exact_reference_verifiers_recomputed": [
            "echosahedral_selector_certificate.py verify",
            "port_current_inner_certificate.py verify",
        ],
        "menu_status": (
            "audit_declared_eight_convex_deltahedra_probe_menu; no source theorem "
            "declares this to be the admissible physical carrier menu"
        ),
        "rows": rows,
        "summary": {
            "declared_probe_menu_size": len(rows),
            "configurations_with_executable_port_response_producer": len(computed),
            "configurations_with_undeclared_port_response_producer": len(undeclared),
            "sm_lie_type_available_count": len(available),
            "sm_lie_type_uniquely_selected_count": len(unique_selected),
            "verdict": (
                "Within A1's declared twelve-port icosahedral boundary and the "
                "A2 same-response internal-holonomy premise, the abstract Lie "
                "type u(1)+su(2)+su(3) is forced. This is an axiom-domain "
                "theorem, not a selection among eight deltahedral carriers. "
                "The seven counterfactual carriers lie outside A1 and have no "
                "compatible producer; they are not counted as exclusions. "
                "Source reconstruction of the matrix current remains open."
            ),
        },
        "twelve_port_completion_nonuniqueness": no_go,
        "matter_nonuniqueness": {
            "rank15_projector_verified": True,
            "inequivalent_completions": no_go["inequivalent_matter_completions"],
            "sterile_singlet_excluded": False,
            "completion_unique": False,
        },
    }


def build_report(root: Path = ROOT) -> dict[str, Any]:
    root = _resolve_root(root)
    hashes = {
        path.as_posix(): _sha256(root / path)
        for path in INPUT_PATHS
    }
    report = {
        "schema": "oph.null_model_scorecard.v1",
        "status": "FAIL_CLOSED_NULL_MODEL_AUDIT",
        "physical_promotion": False,
        "input_sha256": hashes,
        "constant_substitution": build_constant_scan(root),
        "selector_ablation": build_selector_ablation(root),
    }
    validate_report(report)
    return report


def validate_report(report: Mapping[str, Any]) -> None:
    require(report.get("schema") == "oph.null_model_scorecard.v1", "null-model schema mismatch")
    require(report.get("physical_promotion") is False, "null-model audit cannot promote physics")
    constant = report.get("constant_substitution", {})
    require(
        constant.get("status") == "PASS_INTERVAL_CERTIFICATES_COMPLETE",
        "constant-substitution certificates are incomplete",
    )
    rows = constant.get("rows")
    require(isinstance(rows, list), "constant-substitution rows are missing")
    require(
        all(row.get("certificate", {}).get("status") == "INTERVAL_BANACH_ROOT_CERTIFIED" for row in rows),
        "at least one constant-substitution root is uncertified",
    )
    selector = report.get("selector_ablation", {})
    require(
        selector.get("summary", {}).get("sm_lie_type_uniquely_selected_count") == 0,
        "selector report overstates unique Standard Model selection",
    )
    require(
        selector.get("matter_nonuniqueness", {}).get("completion_unique") is False,
        "selector report hides matter completion non-uniqueness",
    )


def render_scorecard(report: Mapping[str, Any]) -> str:
    validate_report(report)
    constant = report["constant_substitution"]
    summary = constant["summary"]
    selector = report["selector_ablation"]

    lines = [
        "# OPH null-model scorecard",
        "",
        "Generated by `python3 tools/check_null_models.py`; do not edit by hand.",
        "This page records negative controls for declared numerical and selector",
        "constructions. It promotes no physical claim.",
        "",
    ]
    lines += [
        "## W3a: constant substitution in the pixel closure",
        "",
        "Map: `P = c₁ + sqrt(c₂)/A_T(P)` using the existing",
        "`thomson_structured_running_plus_gauge_width` chain. Every row below",
        "has an interval-verified inner pixel root, interval-verified `m_Z` root,",
        "edge-sum tail bounds, a centered self-map, and a Banach contraction.",
        "The comparison target is `137.035999177`; distances are relative",
        "`|root−target|/target`.",
        "",
        (
            f"Declared pairs: **{summary['declared_pair_count']}**; numerically "
            f"distinct pairs: **{summary['numerically_distinct_pair_count']}**; "
            f"distinct alternatives after removing the canonical `φ,π` pair: "
            f"**{summary['unique_alternative_pair_count']}**."
        ),
        "",
        "| Threshold | Certified alternative hits | Fraction |",
        "|---:|---:|---:|",
    ]
    for threshold in ("1e-4", "1e-5", "2.5e-6"):
        row = summary["thresholds"][threshold]
        lines.append(
            f"| `{threshold}` | {row['certified_alternative_hits']} / "
            f"{row['unique_alternative_pairs']} | {row['fraction']} |"
        )
    rule = summary["interpretation_rule"]
    lines += [
        "",
        "Interpretation rule: at least "
        f"`{rule['nontrivial_minimum_count']}` distinct alternatives and at least "
        f"`{rule['nontrivial_minimum_fraction']}` of the distinct-alternative menu "
        f"must reach relative distance `{rule['threshold']}` to strip the "
        "`φ`/`√π` detuning story of evidential weight. "
        f"Triggered: **{str(rule['triggered']).lower()}**.",
        "",
        f"Gate interpretation: **{rule['current_interpretation']}**",
        "",
        "The declared grid is an audit probe, not a source-derived, exhaustive,",
        "or preregistered hypothesis menu. Failure of the trigger therefore grants",
        "no positive evidential weight.",
        "",
        f"Threshold provenance: `{rule['threshold']}` is the published OPH-scale",
        "residual of the gauge-width map, the canonical pair's certified relative",
        "distance rounded up to two significant figures. It was set together with",
        "this scan, after that distance was known, so it is not preregistered, and",
        "the canonical pair itself sits about one percent inside it.",
        "",
        "| c₁ | c₂ | Certified inverse-α root (point display) | Absolute-distance enclosure | Relative-distance enclosure | Alias |",
        "|---|---|---:|---:|---:|---|",
    ]
    seen: set[str] = set()
    for row in constant["rows"]:
        numeric_key = row["numeric_pair_key"]
        alias = "duplicate numeric pair" if numeric_key in seen else ""
        seen.add(numeric_key)
        absolute = row["absolute_distance_enclosure"]
        relative = row["relative_distance_enclosure"]
        lines.append(
            f"| `{row['c1']}` | `{row['c2']}` | `{row['alpha_inv_point']}` | "
            f"`[{absolute['lo']}, {absolute['hi']}]` | "
            f"`[{relative['lo']}, {relative['hi']}]` | {alias} |"
        )

    lines += [
        "",
        "## W3b: selector and port-menu ablation",
        "",
        selector["summary"]["verdict"],
        "",
        "| Vertex configuration | Ports | Producer state | SM Lie type |",
        "|---|---:|---|---|",
    ]
    for row in selector["rows"]:
        if row["sm_lie_type_available"] is True:
            sm = "**forced in the A1/A2 domain**; physical matrix current open"
        else:
            sm = "outside A1; no producer"
        lines.append(
            f"| `{row['configuration']}` | {row['port_count']} | "
            f"`{row['status']}` | {sm} |"
        )
    no_go = selector["twelve_port_completion_nonuniqueness"]
    lines += [
        "",
        "For the A1 twelve-port icosahedral carrier, the coefficient-module",
        "classification by itself leaves three viable compact Lie types:",
        "`su(3)+su(2)+u(1)`, `su(2)+su(2)+u(1)^6`, and `u(1)^12`.",
        "The complete faithful A1 response and A2 same-response internal",
        "holonomy supply the stronger fixed-space theorem that forces the first",
        "abstract Lie type. The charged-double-triplet matrices are an exact",
        "declared witness. The response artifact does not reconstruct those",
        "generators, their bracket, or same-current holonomy, so physical source",
        "binding and laboratory gauge-current identification remain open.",
        "",
        "The small completion no-go uses a bare source reduct that omits the",
        "complete-response and internal-holonomy premises. That weaker reduct admits",
        "both",
        f"`{'` and `'.join(no_go['inequivalent_current_completions'])}`.",
        "Its nonuniqueness does not contradict the A1/A2 abstract theorem.",
        "",
        "Matter completion is likewise non-unique: the same reduct admits",
        f"`{'` and `'.join(no_go['inequivalent_matter_completions'])}`. The",
        "rank-15 projector is verified, but the sterile-singlet completion is not",
        "excluded. Seven alternative deltahedra lie outside A1 and lack a compatible",
        "producer. They are counterfactual ablations, not competing OPH branches.",
        "",
        "## Reproducibility",
        "",
        "| Input | SHA-256 |",
        "|---|---|",
    ]
    for path, digest in sorted(report["input_sha256"].items()):
        lines.append(f"| `{path}` | `{digest}` |")
    lines.append("")
    return "\n".join(lines)


def check_scorecard(path: Path, expected: str, report: Mapping[str, Any]) -> None:
    require(path.is_file(), f"null-model scorecard is missing: {path}")
    current = path.read_text(encoding="utf-8")
    require(
        current == expected,
        (
            f"{path} has drifted from the null-model inputs; regenerate with: "
            "python3 tools/check_null_models.py"
        ),
    )


def _output_path(root: Path, requested: Path | None) -> Path:
    if requested is None:
        return root / DEFAULT_SCORECARD
    return requested if requested.is_absolute() else root / requested


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=ROOT,
        help="Repository root to audit (supports isolated mutation tests).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output path, absolute or relative to --root.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Recompute and fail unless the tracked scorecard is byte-exact.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        root = _resolve_root(args.root)
        report = build_report(root)
        rendered = render_scorecard(report)
        output = _output_path(root, args.output)
        if args.check:
            check_scorecard(output, rendered, report)
            print(
                "null-model scorecard OK: "
                f"{len(report['constant_substitution']['rows'])} interval roots, "
                f"{report['selector_ablation']['summary']['declared_probe_menu_size']} "
                "selector probes"
            )
            return 0
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8", newline="\n")
        print(f"wrote {output.relative_to(root)}")
        return 0
    except (NullModelError, ImportError, OSError, ValueError) as exc:
        print(f"null-model check failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
