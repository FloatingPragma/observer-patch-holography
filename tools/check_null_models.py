#!/usr/bin/env python3
"""Generate and check the OPH null-model scorecard.

The three checks are deliberately negative controls:

* substitute the constants in the declared gauge-width P-closure map and
  interval-certify every resulting fixed point;
* expose the declared deltahedral/port menu and distinguish "available" from
  "uniquely selected";
* recompute the RSCC lower-order ablation and require its better fit to appear
  at the top of the generated scorecard.
* calibrate a fully declared shallow expression grammar against 2,000
  deterministic pseudorandom correction targets before any quark replacement
  constant can receive evidential weight.

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
import math
import os
import subprocess
import sys
from dataclasses import dataclass
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCORECARD = Path("tracking/null_model_scorecard.md")
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
RSCC_DISCLOSURE_MARKER = "**RSCC negative-control result: the lower-order ablation fits better.**"
EXPRESSION_NULL_SEED = "oph.w-f6.expression-grammar.v1"
EXPRESSION_NULL_TARGET_COUNT = 2000
EXPRESSION_NULL_LOWER = Decimal("0.6")
EXPRESSION_NULL_UPPER = Decimal("0.8")
EXPRESSION_NULL_TOLERANCES = (
    ("0.20%", Decimal("0.0020")),
    ("0.10%", Decimal("0.0010")),
    ("0.05%", Decimal("0.0005")),
)
EXPRESSION_SMALL_INTEGERS = (1, 2, 3, 4, 5)
QUARK_CORRECTION_TARGETS = (
    ("m_b multiplicative correction", Decimal("0.693432")),
    ("m_s multiplicative correction", Decimal("0.665516")),
    ("m_d multiplicative correction", Decimal("0.768575")),
    ("m_s/m_d bridge factor", Decimal("1.1631")),
)


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


@dataclass(frozen=True)
class ExpressionValue:
    """One syntactic member of the bounded W-F6 expression grammar."""

    expression: str
    value: Decimal
    binary_operation_count: int


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
    Path("code/particles/flavor/quark_rscc_completion_candidate.py"),
    Path("code/particles/flavor/audit_quark_rscc_completion_candidate.py"),
    Path("code/particles/runs/flavor/quark_rscc_completion_candidate_audit.json"),
    Path("code/particles/flavor/derive_down_type_register_clebsch_lane.py"),
    Path("code/particles/runs/flavor/down_type_register_clebsch_lane.json"),
    Path("code/particles/flavor/derive_clebsch_register_pairing_selection.py"),
    Path("code/particles/runs/flavor/clebsch_register_pairing_selection.json"),
    Path("code/particles/data/flag_2024_light_quark_ratio_fixture.json"),
    Path(
        "code/particles/scripts/"
        "generate_flag_2024_light_quark_ratio_fixture.py"
    ),
    Path("code/particles/leptons/derive_charged_mcpr_completion_conditional.py"),
    Path("code/particles/runs/leptons/charged_mcpr_completion_conditional.json"),
    Path("code/particles/leptons/derive_charged_koide_orientation_isometry.py"),
    Path("code/particles/runs/leptons/charged_koide_orientation_isometry.json"),
    Path("code/particles/runs/leptons/charged_stage5_frozen_candidate_audit.json"),
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


def expression_grammar_atoms(
    root: Path = ROOT,
) -> tuple[ExpressionValue, ...]:
    """Return the fully declared terminal alphabet for the W-F6 audit.

    ``P`` is read from the already tracked MCPR source-input receipt.  The
    remaining transcendental constants are fixed decimal literals evaluated
    at substantially more precision than is relevant to the loosest 0.20%
    audit tolerance.
    """

    root = root.resolve()
    mcpr = _load_json(
        root
        / "code/particles/runs/leptons/"
        "charged_mcpr_completion_conditional.json"
    )
    try:
        p_value = Decimal(str(mcpr["source_input"]["P"]))
    except (KeyError, TypeError, ValueError) as exc:
        raise NullModelError(
            "charged MCPR receipt does not expose a decimal source_input.P"
        ) from exc

    with localcontext() as context:
        context.prec = 70
        pi = Decimal(
            "3.141592653589793238462643383279502884197169399375105820974944"
        )
        e = Decimal(
            "2.718281828459045235360287471352662497757247093699959574966968"
        )
        ln_two = Decimal(
            "0.693147180559945309417232121458176568075500134360255254120680"
        )
        sqrt_five = Decimal(5).sqrt()
        phi = (Decimal(1) + sqrt_five) / Decimal(2)
        atoms = [
            ExpressionValue("P", +p_value, 0),
            ExpressionValue("phi", +phi, 0),
            ExpressionValue("pi", +pi, 0),
            ExpressionValue("e", +e, 0),
            ExpressionValue("sqrt(P)", +p_value.sqrt(), 0),
            ExpressionValue("sqrt(pi)", +pi.sqrt(), 0),
            ExpressionValue("ln(2)", +ln_two, 0),
            ExpressionValue("P-phi", +(p_value - phi), 0),
            ExpressionValue("sqrt(5)", +sqrt_five, 0),
        ]
        atoms.extend(
            ExpressionValue(str(integer), Decimal(integer), 0)
            for integer in EXPRESSION_SMALL_INTEGERS
        )
    return tuple(atoms)


def declared_expression_grammar(
    root: Path = ROOT,
) -> tuple[ExpressionValue, ...]:
    """Enumerate the exact, depth-one expression grammar used by W-F6.

    The source audit's phrase "280-expression grammar" does not specify its
    integer range, tree depth, operand ordering, duplicate policy, or value
    window, so reproducing that count would be spurious.  This replacement is
    deliberately simple and exhaustive:

    * the fourteen declared atoms are included;
    * ``+`` and ``*`` use unordered atom pairs with repetition;
    * ``-`` uses ordered distinct atom pairs;
    * ``/`` uses all ordered atom pairs (all atoms are nonzero);
    * no result is fed into another operation.

    This gives 602 syntactic expressions.  Alias-equivalent expressions remain
    present, but aliases cannot change the yes/no question "does any expression
    hit this target?" used by the random-target calibration.
    """

    atoms = expression_grammar_atoms(root)
    expressions: list[ExpressionValue] = list(atoms)
    with localcontext() as context:
        context.prec = 70
        for left_index, left in enumerate(atoms):
            for right_index, right in enumerate(atoms):
                if left_index <= right_index:
                    expressions.extend(
                        (
                            ExpressionValue(
                                f"({left.expression}+{right.expression})",
                                +(left.value + right.value),
                                1,
                            ),
                            ExpressionValue(
                                f"({left.expression}*{right.expression})",
                                +(left.value * right.value),
                                1,
                            ),
                        )
                    )
                if left_index != right_index:
                    expressions.append(
                        ExpressionValue(
                            f"({left.expression}-{right.expression})",
                            +(left.value - right.value),
                            1,
                        )
                    )
                expressions.append(
                    ExpressionValue(
                        f"({left.expression}/{right.expression})",
                        +(left.value / right.value),
                        1,
                    )
                )
    require(
        len(atoms) == 14 and len(expressions) == 602,
        "declared W-F6 expression grammar cardinality drifted",
    )
    return tuple(expressions)


def deterministic_expression_null_targets(
    count: int = EXPRESSION_NULL_TARGET_COUNT,
    seed: str = EXPRESSION_NULL_SEED,
) -> tuple[Decimal, ...]:
    """Generate platform-stable pseudorandom targets in the open interval.

    Each target is a SHA-256 counter-mode draw.  Taking the midpoint of one of
    the ``2**64`` equal bins excludes both interval endpoints and avoids a
    dependency on the implementation details of Python's ``random`` module.
    """

    require(count > 0, "expression-null target count must be positive")
    denominator = Decimal(2**64)
    span = EXPRESSION_NULL_UPPER - EXPRESSION_NULL_LOWER
    targets: list[Decimal] = []
    with localcontext() as context:
        context.prec = 70
        for index in range(count):
            digest = hashlib.sha256(
                f"{seed}:{index}".encode("utf-8")
            ).digest()
            integer = int.from_bytes(digest[:8], "big")
            unit = (Decimal(integer) + Decimal("0.5")) / denominator
            target = EXPRESSION_NULL_LOWER + span * unit
            require(
                EXPRESSION_NULL_LOWER
                < target
                < EXPRESSION_NULL_UPPER,
                "expression-null target escaped the open interval",
            )
            targets.append(+target)
    return tuple(targets)


def _relative_expression_distance(
    value: Decimal,
    target: Decimal,
) -> Decimal:
    require(target != 0, "expression-null target cannot be zero")
    with localcontext() as context:
        context.prec = 70
        return +(abs(value - target) / abs(target))


def build_expression_grammar_null_model(
    root: Path = ROOT,
) -> dict[str, Any]:
    """Calibrate the declared correction-expression search against null targets."""

    expressions = declared_expression_grammar(root)
    atoms = expression_grammar_atoms(root)
    targets = deterministic_expression_null_targets()

    thresholds: dict[str, dict[str, Any]] = {}
    for display, tolerance in EXPRESSION_NULL_TOLERANCES:
        hit_count = sum(
            any(
                _relative_expression_distance(expression.value, target)
                <= tolerance
                for expression in expressions
            )
            for target in targets
        )
        with localcontext() as context:
            context.prec = 50
            fraction = Decimal(hit_count) / Decimal(len(targets))
        thresholds[display] = {
            "relative_tolerance": format(tolerance, "f"),
            "random_target_hit_count": hit_count,
            "random_target_count": len(targets),
            "random_target_hit_fraction": format(fraction, "f"),
            "random_target_hit_percent": format(
                fraction * Decimal(100), "f"
            ),
        }

    correction_rows: list[dict[str, Any]] = []
    for label, target in QUARK_CORRECTION_TARGETS:
        ranked = sorted(
            (
                (
                    _relative_expression_distance(expression.value, target),
                    expression.binary_operation_count,
                    expression.expression,
                    expression,
                )
                for expression in expressions
            ),
            key=lambda row: row[:3],
        )
        distance, _, _, nearest = ranked[0]
        correction_rows.append(
            {
                "label": label,
                "audit_supplied_rounded_target": format(target, "f"),
                "inside_random_calibration_interval": (
                    EXPRESSION_NULL_LOWER
                    < target
                    < EXPRESSION_NULL_UPPER
                ),
                "nearest_expression": nearest.expression,
                "nearest_value": format(nearest.value, "f"),
                "nearest_relative_distance": format(distance, "f"),
                "nearest_relative_distance_percent": format(
                    distance * Decimal(100), "f"
                ),
                "threshold_hits": {
                    display: distance <= tolerance
                    for display, tolerance in EXPRESSION_NULL_TOLERANCES
                },
                "target_informed_diagnostic": True,
                "promotion_allowed": False,
            }
        )

    positive_finite_count = sum(expression.value > 0 for expression in expressions)
    interval_expression_count = sum(
        EXPRESSION_NULL_LOWER
        < expression.value
        < EXPRESSION_NULL_UPPER
        for expression in expressions
    )
    return {
        "status": "BOUNDED_DECLARED_GRAMMAR_NULL_COMPLETE_NO_EVIDENTIAL_WEIGHT",
        "promotion_allowed": False,
        "evidential_weight_granted": False,
        "source_claimed_280_expression_rates_reproduced": False,
        "source_grammar_reconstruction_status": (
            "NOT_REPRODUCIBLE_AS_STATED: the cited 280-expression menu does "
            "not declare its small-integer range, tree depth, operand ordering, "
            "duplicate equivalence, or value filter. The rates below are fresh "
            "results for the fully declared grammar in this receipt."
        ),
        "grammar": {
            "schema": "oph.w_f6.expression_grammar.v1",
            "atoms": [atom.expression for atom in atoms],
            "small_integer_set": list(EXPRESSION_SMALL_INTEGERS),
            "binary_operators": ["+", "-", "*", "/"],
            "maximum_binary_operations": 1,
            "commutative_pair_policy": (
                "unordered_with_repetition_for_plus_and_multiply"
            ),
            "noncommutative_pair_policy": (
                "ordered_distinct_for_subtract; ordered_with_repetition_for_divide"
            ),
            "recursive_composition": False,
            "syntactic_expression_count": len(expressions),
            "positive_expression_count": positive_finite_count,
            "expressions_inside_random_target_interval": interval_expression_count,
            "numeric_alias_policy": (
                "syntactic aliases retained; aliases cannot change any-hit rates"
            ),
            "P_source": (
                "code/particles/runs/leptons/"
                "charged_mcpr_completion_conditional.json#/source_input/P"
            ),
            "P_value": next(
                format(atom.value, "f")
                for atom in atoms
                if atom.expression == "P"
            ),
        },
        "random_target_sampler": {
            "algorithm": (
                "SHA-256(seed + ':' + decimal_counter), first 64 digest bits, "
                "open-bin midpoint mapped affinely to (0.6,0.8)"
            ),
            "seed": EXPRESSION_NULL_SEED,
            "target_count": len(targets),
            "open_interval": [
                format(EXPRESSION_NULL_LOWER, "f"),
                format(EXPRESSION_NULL_UPPER, "f"),
            ],
        },
        "thresholds": thresholds,
        "candidate_correction_diagnostics": correction_rows,
        "claim_entry_gate": {
            "required_for_any_future_proposed_constant": [
                "exact finite grammar",
                "relative tolerance",
                "deterministic null-target hit rate",
            ],
            "current_candidate_corrections_promotable": False,
            "verdict": (
                "NO_EVIDENTIAL_WEIGHT: proximity to a short expression is a "
                "target-informed diagnostic until a source derivation and a "
                "preregistered grammar/tolerance test both exist."
            ),
        },
    }


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
        current_gate.get("physical_source_realization_gate_passed") is True,
        "12-port current receipt lost its source-bound response artifact gate",
    )
    require(
        isinstance(current.get("semantic_response_binding"), dict)
        and current["semantic_response_binding"].get("sector_structure_recomputed")
        is True,
        "12-port current gate passes without a recomputed semantic binding",
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
                "physical_source_binding": True,
                "status": "AVAILABLE_NOT_UNIQUELY_SELECTED",
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
                "u(1)+su(2)+su(3) is available on the source-bound twelve-port "
                "response branch, not uniquely selected across the carrier "
                "menu. Seven alternative deltahedra have no compatible "
                "producer and are unknown; they are not counted as exclusions."
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


def _load_rscc_audit_module(root: Path) -> Any:
    module_dir = root / "code/particles/flavor"
    audit_path = module_dir / "audit_quark_rscc_completion_candidate.py"
    old_candidate = sys.modules.pop("quark_rscc_completion_candidate", None)
    sys.path.insert(0, str(module_dir))
    try:
        module = _load_module(audit_path, f"_oph_rscc_null_models_{hash(audit_path)}")
    finally:
        sys.path.remove(str(module_dir))
        sys.modules.pop("quark_rscc_completion_candidate", None)
        if old_candidate is not None:
            sys.modules["quark_rscc_completion_candidate"] = old_candidate
    return module


def build_rscc_ablation(root: Path) -> dict[str, Any]:
    module = _load_rscc_audit_module(root)
    recomputed = module.build_audit()
    stored_path = (
        root / "code/particles/runs/flavor/quark_rscc_completion_candidate_audit.json"
    )
    stored = _load_json(stored_path)
    full = recomputed["descriptive_mixed_chart_comparison"]
    ablation = recomputed["negative_control"]
    stored_full = stored.get("descriptive_mixed_chart_comparison", {})
    stored_ablation = stored.get("negative_control", {})

    metric_names = ("max_abs_relative_error_percent", "raw_diagonal_residual_sum")
    for metric in metric_names:
        require(full.get(metric) == stored_full.get(metric), f"stored RSCC full {metric} drifted")
        require(
            ablation.get(metric) == stored_ablation.get(metric),
            f"stored RSCC ablation {metric} drifted",
        )
    beats_max = (
        float(ablation["max_abs_relative_error_percent"])
        < float(full["max_abs_relative_error_percent"])
    )
    beats_sum = (
        float(ablation["raw_diagonal_residual_sum"])
        < float(full["raw_diagonal_residual_sum"])
    )
    require(
        ablation.get("beats_full_rscc_maximum_error") is beats_max,
        "RSCC maximum-error ablation verdict is inconsistent with its metrics",
    )
    require(
        ablation.get("beats_full_rscc_raw_residual_sum") is beats_sum,
        "RSCC residual-sum ablation verdict is inconsistent with its metrics",
    )
    require(
        stored_ablation.get("beats_full_rscc_maximum_error") is beats_max,
        "stored RSCC maximum-error ablation verdict drifted",
    )
    require(
        stored_ablation.get("beats_full_rscc_raw_residual_sum") is beats_sum,
        "stored RSCC residual-sum ablation verdict drifted",
    )
    return {
        "status": "NEGATIVE_CONTROL_BEATS_FULL_RSCC" if beats_max or beats_sum else "FULL_RSCC_NOT_BEATEN",
        "full_rscc": {
            "max_abs_relative_error_percent": full["max_abs_relative_error_percent"],
            "raw_diagonal_residual_sum": full["raw_diagonal_residual_sum"],
        },
        "zero_w2_zero_delta_g_ablation": {
            "max_abs_relative_error_percent": ablation[
                "max_abs_relative_error_percent"
            ],
            "raw_diagonal_residual_sum": ablation["raw_diagonal_residual_sum"],
        },
        "ablation_beats_full_maximum_error": beats_max,
        "ablation_beats_full_residual_sum": beats_sum,
        "ablation_beats_full_model": beats_max or beats_sum,
        "interpretation": (
            "The detailed RSCC covariance ledger is not selected by its own "
            "retrospective comparison metric. These mixed-chart residuals are "
            "diagnostics, not likelihoods."
        ),
    }


def _load_down_type_clebsch_module(root: Path) -> Any:
    """Load the lane from an arbitrary repository root without stale imports."""

    flavor_dir = root / "code/particles/flavor"
    particles_dir = root / "code/particles"
    calibration_name = "calibration.derive_d11_criticality_boundary_scan"
    fallback_name = "derive_d11_criticality_boundary_scan"
    saved_modules = {
        name: sys.modules.pop(name, None)
        for name in (calibration_name, fallback_name)
    }
    old_path = sys.path[:]
    sys.path.insert(0, str(particles_dir))
    sys.path.insert(0, str(flavor_dir))
    path = flavor_dir / "derive_down_type_register_clebsch_lane.py"
    try:
        module = _load_module(path, f"_oph_clebsch_null_models_{hash(path)}")
    finally:
        sys.path[:] = old_path
        for name in (calibration_name, fallback_name):
            sys.modules.pop(name, None)
            if saved_modules[name] is not None:
                sys.modules[name] = saved_modules[name]
    return module


def _relative_decimal_gap(left: Decimal, right: Decimal) -> Decimal:
    require(right != 0, "relative-gap reference is zero")
    return abs(left / right - Decimal(1))


def _koide_provenance_audit(
    mcpr: Mapping[str, Any],
    finite_gns: Mapping[str, Any],
    stage5_audit: Mapping[str, Any],
) -> dict[str, Any]:
    """Classify the square-root-mass invariant without promoting the bridge."""

    with localcontext() as context:
        context.prec = 90
        chi = Decimal(str(mcpr["response_map"]["fixed_point"]["chi"]))
        direct = Decimal(str(mcpr["regular_C3_shape"]["sqrt_mass_invariant"]))
        formula = (Decimal(1) + (-Decimal(2) * chi).exp()) / Decimal(3)
        two_thirds = Decimal(2) / Decimal(3)
        formula_gap = abs(direct - formula)
        two_thirds_gap = abs(direct - two_thirds)

    interpretation = str(
        stage5_audit.get("koide_invariant", {}).get("interpretation", "")
    )
    require(
        "imposed balanced carrier" in interpretation,
        "Stage-5 audit no longer exposes the imposed balanced-carrier boundary",
    )
    require(
        mcpr.get("charged_reference_data_consumed") is False
        and mcpr.get("public_charged_mass_promotion_allowed") is False,
        "MCPR artifact changed its no-target/promotion boundary",
    )
    require(
        mcpr.get("checks", {}).get("sqrt_mass_invariant_identity") is True,
        "MCPR square-root-mass identity check no longer passes",
    )
    require(
        formula_gap < Decimal("1e-58"),
        "MCPR square-root-mass invariant no longer equals its declared formula",
    )
    require(
        finite_gns.get("public_koide_promotion_allowed") is False,
        "finite GNS artifact unexpectedly promotes physical Koide",
    )
    require(
        finite_gns.get("checks_pass") is True,
        "finite GNS conditional theorem no longer passes",
    )
    return {
        "mcpr_sqrt_mass_invariant": format(direct, "f"),
        "mcpr_formula": "(1 + exp(-2*chi))/3",
        "mcpr_chi": format(chi, "f"),
        "formula_identity_absolute_gap": format(formula_gap, "E"),
        "distance_from_two_thirds": format(two_thirds_gap, "E"),
        "mcpr_is_exact_two_thirds": direct == two_thirds,
        "runtime_charged_reference_consumed": False,
        "finite_gns_status": finite_gns.get("status"),
        "finite_gns_public_koide_promotion_allowed": False,
        "stage5_audit_interpretation": interpretation,
        "classification": (
            "CONDITIONAL_ALGEBRAIC_CONSEQUENCE_OF_DECLARED_ARCHITECTURE: "
            "the direct invariant is derived from the MCPR roots and equals "
            "(1+exp(-2 chi))/3 without a charged target input, but the response "
            "architecture, phase, and balanced amplitude are stipulated and "
            "historically target-informed rather than blinded. The companion "
            "finite tracial-GNS theorem derives balance inside its conditional "
            "event model; the physical chiral mass-response attachment remains "
            "open. Therefore physical Koide is neither a runtime target imposed "
            "on the MCPR solve nor a blind/source-derived physical prediction."
        ),
    }


def _flag_2024_compare_only(
    fixture: Mapping[str, Any],
    oph_ms_over_md: Decimal,
) -> dict[str, Any]:
    """Derive FLAG ms/md and apply the declared conservative comparison gate."""

    require(
        fixture.get("schema")
        == "oph.flag_2024_light_quark_ratio_fixture.v1",
        "FLAG light-quark fixture schema drifted",
    )
    require(
        fixture.get("status") == "COMPARE_ONLY_HAND_TRANSCRIBED_REFERENCE",
        "FLAG light-quark fixture is not compare-only",
    )
    boundary = fixture.get("claim_boundary", {})
    require(
        boundary.get("comparison_only") is True
        and boundary.get("oph_fit_or_selection_input") is False
        and boundary.get("oph_theory_uncertainty_supplied") is False
        and boundary.get("prediction_preexisted_audit") is True
        and boundary.get("significance_gate_preregistered") is False
        and boundary.get("comparison_is_retrospective") is True,
        "FLAG fixture crossed its compare-only/no-theory-uncertainty boundary",
    )
    derived = fixture.get("derived_quantity", {})
    require(
        derived.get("identity")
        == "ms/md = (ms/mud) * (1 + mu/md) / 2",
        "FLAG derived-ratio identity drifted",
    )
    require(
        derived.get("input_covariance_available") is False,
        "FLAG fixture unexpectedly claims an input covariance",
    )
    policy = derived.get("uncertainty_policy", {})
    threshold = Decimal(
        str(policy.get("conservative_rejection_threshold_sigma"))
    )
    require(
        threshold > 0,
        "FLAG comparison rejection threshold must be positive",
    )

    source_rows = fixture.get("averages")
    require(
        isinstance(source_rows, list) and len(source_rows) == 2,
        "FLAG fixture must contain the two declared Nf averages",
    )
    rows: list[dict[str, Any]] = []
    for source in source_rows:
        with localcontext() as context:
            context.prec = 60
            ms_mud = Decimal(str(source["ms_over_mud"]["value"]))
            sigma_ms_mud = Decimal(
                str(source["ms_over_mud"]["standard_uncertainty"])
            )
            mu_md = Decimal(str(source["mu_over_md"]["value"]))
            sigma_mu_md = Decimal(
                str(source["mu_over_md"]["standard_uncertainty"])
            )
            central = ms_mud * (Decimal(1) + mu_md) / Decimal(2)
            from_ms_mud = (
                (Decimal(1) + mu_md) / Decimal(2) * sigma_ms_mud
            )
            from_mu_md = ms_mud / Decimal(2) * sigma_mu_md
            independent_sigma = (
                from_ms_mud * from_ms_mud
                + from_mu_md * from_mu_md
            ).sqrt()
            positive_sigma = from_ms_mud + from_mu_md
            absolute_gap = abs(oph_ms_over_md - central)
            independent_gap_sigma = absolute_gap / independent_sigma
            conservative_gap_sigma = absolute_gap / positive_sigma
        stored_derived = source.get("derived_ms_over_md", {})
        require(
            stored_derived
            == {
                "value": format(central, "f"),
                "independent_standard_uncertainty": format(
                    independent_sigma, "f"
                ),
                "rho_plus_one_standard_uncertainty": format(
                    positive_sigma, "f"
                ),
            },
            f"FLAG derived ms/md fields drifted for Nf={source.get('nf')}",
        )
        rows.append(
            {
                "nf": source["nf"],
                "transcribed_inputs": {
                    "ms_over_mud": source["ms_over_mud"],
                    "mu_over_md": source["mu_over_md"],
                },
                "derived_ms_over_md": format(central, "f"),
                "uncertainty_contributions": {
                    "from_ms_over_mud": format(from_ms_mud, "f"),
                    "from_mu_over_md": format(from_mu_md, "f"),
                },
                "independent_propagation": {
                    "standard_uncertainty": format(
                        independent_sigma, "f"
                    ),
                    "gap_sigma": format(independent_gap_sigma, "f"),
                },
                "maximally_positive_correlation_propagation": {
                    "rho": "1",
                    "standard_uncertainty": format(positive_sigma, "f"),
                    "gap_sigma": format(conservative_gap_sigma, "f"),
                },
                "conservative_rejection_gate": {
                    "threshold_sigma": format(threshold, "f"),
                    "gap_sigma": format(conservative_gap_sigma, "f"),
                    "rejected": conservative_gap_sigma >= threshold,
                },
            }
        )

    require(
        {row["nf"] for row in rows} == {"2+1+1", "2+1"},
        "FLAG fixture Nf menu drifted",
    )
    require(
        all(
            row["conservative_rejection_gate"]["rejected"]
            for row in rows
        ),
        "conditional Clebsch lane is not rejected by every FLAG comparison",
    )
    return {
        "status": "CONDITIONAL_LANE_OUTPUT_REJECTED_BY_FLAG_COMPARE_ONLY_GATE",
        "source": fixture["source"],
        "oph_conditional_ms_over_md": format(oph_ms_over_md, "f"),
        "input_covariance_available": False,
        "oph_theory_uncertainty_supplied": False,
        "prediction_solve_or_physical_selection_input": False,
        "prediction_preexisted_audit": True,
        "significance_gate_preregistered": False,
        "comparison_is_retrospective": True,
        "uncertainty_policy": {
            "independent_case_reported": True,
            "rho_plus_one_case_reported": True,
            "conservative_rejection_gate_uses_rho_plus_one": True,
            "threshold_sigma": format(threshold, "f"),
            "qualification": (
                "FLAG does not supply the covariance between the two "
                "transcribed marginal ratios. Independent propagation and "
                "rho=+1 propagation are therefore separate declared cases, "
                "not bounds inferred from a hidden likelihood."
            ),
        },
        "rows": rows,
        "rejection_gate_triggered_for_all_declared_nf_rows": True,
        "interpretation": (
            "The current conditional Clebsch lane output is falsified by this "
            "retrospective >=5-sigma compare-only rejection gate even under "
            "the larger rho=+1 propagated uncertainty. The prediction "
            "preexisted this audit, but the significance rule was not "
            "preregistered. No OPH theory uncertainty is fabricated."
        ),
    }


def build_quark_clebsch_audit(root: Path) -> dict[str, Any]:
    """Enumerate all six factor assignments and audit the protected ratio."""

    lane = _load_down_type_clebsch_module(root)
    mcpr_path = (
        root / "code/particles/runs/leptons/charged_mcpr_completion_conditional.json"
    )
    stored_path = (
        root / "code/particles/runs/flavor/down_type_register_clebsch_lane.json"
    )
    selection_path = (
        root / "code/particles/runs/flavor/clebsch_register_pairing_selection.json"
    )
    finite_gns_path = (
        root / "code/particles/runs/leptons/charged_koide_orientation_isometry.json"
    )
    stage5_audit_path = (
        root / "code/particles/runs/leptons/charged_stage5_frozen_candidate_audit.json"
    )
    flag_path = (
        root / "code/particles/data/flag_2024_light_quark_ratio_fixture.json"
    )
    matter_receipt_path = (
        root
        / "code/a5_closure/receipts/"
        "super_tannakian_matter_reference.receipt.json"
    )
    mcpr = _load_json(mcpr_path)
    stored = _load_json(stored_path)
    selection = _load_json(selection_path)
    finite_gns = _load_json(finite_gns_path)
    stage5_audit = _load_json(stage5_audit_path)
    flag_fixture = _load_json(flag_path)
    matter_receipt = _load_json(matter_receipt_path)

    require(
        selection.get("weight_set_scan", {}).get("surviving_weight_set")
        == ["1/3", "1", "3"],
        "Clebsch selection artifact no longer emits the declared unordered set",
    )
    require(
        selection.get("weight_set_scan", {})
        .get("order_assignment", {})
        .get("status")
        == "open",
        "generation-order assignment unexpectedly closed",
    )
    require(
        stored.get("promotion_allowed") is False,
        "stored Clebsch lane unexpectedly allows promotion",
    )
    yukawa_sector = matter_receipt.get("yukawa_sector", {})
    channel_dimensions = {
        "/".join(channel["channel"]): channel["invariant_dimension"]
        for channel in yukawa_sector.get("channels", [])
    }
    forbidden = yukawa_sector.get("forbidden_channel_control", {})
    require(
        channel_dimensions.get("Q/Sbar/d_c") == 1
        and channel_dimensions.get("L/Sbar/e_c") == 1
        and forbidden.get("channel") == ["Q", "S", "d_c"]
        and forbidden.get("invariant_dimension") == 0,
        "matter receipt Yukawa-line facts drifted",
    )

    factor_names = ("b_over_tau", "s_over_mu", "d_over_e")
    artifact = lane.build_artifact(
        mcpr,
        selection,
        flag_fixture,
        input_hashes=stored["dependency_audit"]["input_sha256"],
    )
    require(
        artifact["predictions"] == stored["predictions"]
        and artifact["permutation_scan"] == stored["permutation_scan"],
        "stored down-type Clebsch artifact drifted from its producer",
    )
    rows: list[dict[str, Any]] = []
    for source_row in artifact["permutation_scan"]["rows"]:
        boundary = source_row["boundary_values_at_mu_U"]
        predictions = source_row["predictions"]
        transport = source_row["common_transport"]
        with localcontext() as context:
            context.prec = 50
            output_ratio = Decimal(str(predictions["ms_over_md"]))
            protected_rhs = Decimal(
                str(transport["register_scale_ratio_identity_rhs"])
            )
            protected_relative_gap = _relative_decimal_gap(
                output_ratio, protected_rhs
            )
        rows.append(
            {
                "assignment": source_row["assignment_labels"],
                "adopted_assignment": source_row[
                    "is_current_assumed_order"
                ],
                "retrospective_metric": source_row[
                    "retrospective_metric"
                ],
                "retrospective_unique_least_discrepant": source_row[
                    "retrospective_unique_least_discrepant"
                ],
                "conservative_flag_rejected_for_all_nf_rows": source_row[
                    "conservative_flag_rejected_for_all_nf_rows"
                ],
                "boundary_values_at_mu_U": {
                    key: boundary[key]
                    for key in (
                        "y_tau",
                        "y_mu",
                        "y_e",
                        "y_b",
                        "y_s",
                        "y_d",
                    )
                },
                "conditional_outputs": {
                    key: predictions[key]
                    for key in (
                        "mb_mb_gev",
                        "ms_2gev_gev",
                        "md_2gev_gev",
                        "ms_over_md",
                        "cabibbo_gst_sqrt_md_over_ms",
                    )
                },
                "arithmetic_checks_pass": all(
                    source_row["arithmetic_checks"].values()
                ),
                "rg_ratio_identity": {
                    "factor_ratio_s_over_d": transport[
                        "factor_ratio_s_over_d"
                    ],
                    "output_ms_over_md": format(output_ratio, "f"),
                    "boundary_rhs": format(protected_rhs, "f"),
                    "relative_gap": format(
                        protected_relative_gap, "E"
                    ),
                    "passes": protected_relative_gap < Decimal("1e-13"),
                },
            }
        )

    require(len(rows) == 6, "Clebsch permutation menu is incomplete")
    require(
        all(row["rg_ratio_identity"]["passes"] for row in rows),
        "common down-lane RG factor does not protect every permutation ratio",
    )
    require(
        all(row["arithmetic_checks_pass"] for row in rows),
        "at least one Clebsch permutation failed arithmetic checks",
    )
    winners = [
        row
        for row in rows
        if row["retrospective_unique_least_discrepant"]
    ]
    require(
        len(winners) == 1,
        "retrospective discrepancy metric is not uniquely minimized",
    )
    winner = winners[0]
    require(
        winner["adopted_assignment"] is True,
        "current assumed order is not retrospectively least discrepant",
    )
    require(
        all(
            row["conservative_flag_rejected_for_all_nf_rows"]
            for row in rows
        ),
        "at least one declared Clebsch order evades the FLAG rejection gate",
    )

    adopted = next(row for row in rows if row["adopted_assignment"])

    input_masses = [
        Decimal(str(value))
        for value in mcpr["optional_scale_display"]["masses_MeV"]
    ]
    with localcontext() as context:
        context.prec = 60
        input_low_scale_rhs = input_masses[1] / input_masses[0] / Decimal(9)
        output_ratio = Decimal(
            adopted["rg_ratio_identity"]["output_ms_over_md"]
        )
        low_scale_relative_gap = _relative_decimal_gap(
            output_ratio, input_low_scale_rhs
        )
    scale_qualified_identity = {
        "exact_algebraic_statement": (
            "m_s(2 GeV)/m_d(2 GeV) = "
            "[y_mu(mu_U)/y_e(mu_U)]/9"
        ),
        "factor_origin": "(1/3)/3 = 1/9",
        "common_down_lane_rg_and_qcd_factor_cancels": True,
        "passes_for_all_six_permutations": True,
        "adopted_output_ms_over_md": adopted["rg_ratio_identity"][
            "output_ms_over_md"
        ],
        "adopted_register_scale_rhs": adopted["rg_ratio_identity"][
            "boundary_rhs"
        ],
        "adopted_register_scale_relative_gap": adopted[
            "rg_ratio_identity"
        ]["relative_gap"],
        "literal_low_scale_mcpr_input_rhs": format(input_low_scale_rhs, "f"),
        "literal_low_scale_relative_gap": format(
            low_scale_relative_gap, "E"
        ),
        "literal_low_scale_identity_exact": low_scale_relative_gap
        < Decimal("1e-13"),
        "qualification": (
            "The identity is RG-protected only after the charged-lepton ratio "
            "is evaluated at the Clebsch boundary mu_U. The current lepton "
            "runner contains flavor-dependent 1.5*y_i^2 self terms, so replacing "
            "that boundary ratio with the low-scale MCPR display m_mu/m_e is not "
            "an exact identity in the implementation. Common down-lane running "
            "and same-scale QCD transport cancel; a generation-dependent "
            "threshold or matching packet would break that cancellation and "
            "therefore define a modified lane, not a protected consequence of "
            "the current one."
        ),
    }
    require(
        scale_qualified_identity["literal_low_scale_identity_exact"] is False,
        "low-scale and register-scale lepton ratios unexpectedly became identical",
    )

    return {
        "status": "BOUNDED_SIX_PERMUTATION_AUDIT_COMPLETE",
        "promotion_allowed": False,
        "menu": {
            "factor_set": ["1", "1/3", "3"],
            "assignment_slots": list(factor_names),
            "permutation_count": 6,
            "exhaustive": True,
        },
        "retrospective_target_informed_metric": {
            "order": artifact["permutation_scan"][
                "retrospective_metric"
            ]["order"],
            "target_informed": True,
            "preregistered": False,
            "current_assumed_order_uniquely_least_discrepant": True,
            "physical_generation_order_selected": False,
            "remaining_open_premises": [
                "REGISTER_RELATION_EXISTENCE",
                "INDEPENDENT_YUKAWA_COEFFICIENT_IDENTIFICATION",
                "GENERATION_REGISTER_ORDER",
                "THIRD_GENERATION_REGISTER_FACTOR",
            ],
            "all_six_permutations_rejected": True,
        },
        "pairing_and_weight_set_provenance": {
            "matter_receipt_exact_conditional_facts": {
                "Q/Sbar/d_c_invariant_dimension": 1,
                "L/Sbar/e_c_invariant_dimension": 1,
                "Q/S/d_c_invariant_dimension": 0,
            },
            "constrains_candidate_channel_pairing": True,
            "equates_independent_yukawa_coefficients": False,
            "pairing_status": (
                "CONDITIONAL_CHANNEL_COMPATIBILITY_ONLY: the two allowed "
                "Sbar invariant lines and the forbidden S control constrain "
                "which candidate channels can be paired. They do not identify "
                "the independent down-type and charged-lepton Yukawa "
                "coefficients."
            ),
            "F1_F2_status": (
                "CONDITIONAL_DECLARED_ALGEBRA: F1 measure balance and F2 "
                "register faithfulness select the unordered {1/3,1,3} set "
                "inside the declared alphabet. They remain frozen constraints, "
                "not source-derived physical laws."
            ),
        },
        "rows": rows,
        "scale_qualified_rg_identity": scale_qualified_identity,
        "flag_2024_compare_only": _flag_2024_compare_only(
            flag_fixture,
            Decimal(adopted["rg_ratio_identity"]["output_ms_over_md"]),
        ),
        "koide_and_sqrt_mass_invariant": _koide_provenance_audit(
            mcpr, finite_gns, stage5_audit
        ),
        "claim_boundary": (
            "This is an exhaustive audit of one declared six-element assignment "
            "menu. The current assumed order is uniquely least discrepant under "
            "a retrospective, target-informed metric; this does not select a "
            "physical family order, and all six assignments fail the conservative "
            "FLAG gate. The rejection applies to the current declared-model "
            "route, not every possible generation-dependent threshold lane. The "
            "prediction preexisted this audit, but the significance gate was not "
            "preregistered. No external or OPH uncertainty is manufactured."
        ),
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
        "rscc_ablation": build_rscc_ablation(root),
        "quark_clebsch_audit": build_quark_clebsch_audit(root),
        "expression_grammar_null_model": (
            build_expression_grammar_null_model(root)
        ),
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
    rscc = report.get("rscc_ablation", {})
    if rscc.get("ablation_beats_full_model"):
        require(
            rscc.get("status") == "NEGATIVE_CONTROL_BEATS_FULL_RSCC",
            "RSCC negative-control win is not reflected in status",
        )
    clebsch = report.get("quark_clebsch_audit", {})
    require(
        clebsch.get("status") == "BOUNDED_SIX_PERMUTATION_AUDIT_COMPLETE",
        "bounded Clebsch permutation audit is incomplete",
    )
    require(
        clebsch.get("menu", {}).get("permutation_count") == 6
        and len(clebsch.get("rows", [])) == 6,
        "Clebsch assignment menu does not contain six permutations",
    )
    retrospective = clebsch.get(
        "retrospective_target_informed_metric", {}
    )
    require(
        sum(
            bool(row.get("retrospective_unique_least_discrepant"))
            for row in clebsch["rows"]
        )
        == 1
        and retrospective.get(
            "current_assumed_order_uniquely_least_discrepant"
        )
        is True
        and retrospective.get("target_informed") is True
        and retrospective.get("preregistered") is False
        and retrospective.get("physical_generation_order_selected") is False
        and retrospective.get("all_six_permutations_rejected") is True,
        "Clebsch retrospective metric boundary or all-miss result drifted",
    )
    require(
        clebsch["scale_qualified_rg_identity"][
            "passes_for_all_six_permutations"
        ]
        is True,
        "Clebsch RG-ratio identity is not protected over the full menu",
    )
    flag = clebsch.get("flag_2024_compare_only", {})
    require(
        flag.get("status")
        == "CONDITIONAL_LANE_OUTPUT_REJECTED_BY_FLAG_COMPARE_ONLY_GATE"
        and flag.get("rejection_gate_triggered_for_all_declared_nf_rows")
        is True,
        "FLAG compare-only falsification gate is not closed",
    )
    require(
        flag.get("input_covariance_available") is False
        and flag.get("oph_theory_uncertainty_supplied") is False
        and flag.get("prediction_solve_or_physical_selection_input") is False
        and flag.get("prediction_preexisted_audit") is True
        and flag.get("significance_gate_preregistered") is False
        and flag.get("comparison_is_retrospective") is True,
        "FLAG comparison crossed its external-input claim boundary",
    )
    expression_null = report.get("expression_grammar_null_model", {})
    require(
        expression_null.get("status")
        == "BOUNDED_DECLARED_GRAMMAR_NULL_COMPLETE_NO_EVIDENTIAL_WEIGHT",
        "W-F6 expression-grammar null audit is incomplete",
    )
    require(
        expression_null.get("promotion_allowed") is False
        and expression_null.get("evidential_weight_granted") is False
        and expression_null.get(
            "source_claimed_280_expression_rates_reproduced"
        )
        is False,
        "W-F6 expression audit crossed its fail-closed claim boundary",
    )
    grammar = expression_null.get("grammar", {})
    require(
        grammar.get("schema") == "oph.w_f6.expression_grammar.v1"
        and grammar.get("small_integer_set") == [1, 2, 3, 4, 5]
        and grammar.get("binary_operators") == ["+", "-", "*", "/"]
        and grammar.get("maximum_binary_operations") == 1
        and grammar.get("recursive_composition") is False
        and grammar.get("syntactic_expression_count") == 602,
        "W-F6 finite grammar declaration drifted",
    )
    sampler = expression_null.get("random_target_sampler", {})
    require(
        sampler.get("seed") == EXPRESSION_NULL_SEED
        and sampler.get("target_count") == EXPRESSION_NULL_TARGET_COUNT
        and sampler.get("open_interval") == ["0.6", "0.8"],
        "W-F6 deterministic random-target sampler drifted",
    )
    expression_thresholds = expression_null.get("thresholds", {})
    hit_counts: list[int] = []
    for display, tolerance in EXPRESSION_NULL_TOLERANCES:
        row = expression_thresholds.get(display, {})
        require(
            row.get("relative_tolerance") == format(tolerance, "f")
            and row.get("random_target_count")
            == EXPRESSION_NULL_TARGET_COUNT
            and isinstance(row.get("random_target_hit_count"), int),
            f"W-F6 threshold row is incomplete: {display}",
        )
        hit_counts.append(row["random_target_hit_count"])
    require(
        hit_counts == sorted(hit_counts, reverse=True),
        "W-F6 hit counts are not monotone in tolerance",
    )
    corrections = expression_null.get(
        "candidate_correction_diagnostics", []
    )
    require(
        len(corrections) == len(QUARK_CORRECTION_TARGETS)
        and all(
            row.get("target_informed_diagnostic") is True
            and row.get("promotion_allowed") is False
            for row in corrections
        )
        and expression_null.get("claim_entry_gate", {}).get(
            "current_candidate_corrections_promotable"
        )
        is False,
        "W-F6 correction diagnostics did not fail closed",
    )


def _fmt_decimal(value: Any, digits: int = 12) -> str:
    return f"{float(value):.{digits}g}"


def render_scorecard(report: Mapping[str, Any]) -> str:
    validate_report(report)
    constant = report["constant_substitution"]
    summary = constant["summary"]
    selector = report["selector_ablation"]
    rscc = report["rscc_ablation"]
    clebsch = report["quark_clebsch_audit"]
    expression_null = report["expression_grammar_null_model"]

    lines = [
        "# OPH null-model scorecard",
        "",
        "Generated by `python3 tools/check_null_models.py`; do not edit by hand.",
        "This page records negative controls for declared numerical and selector",
        "constructions. It promotes no physical claim.",
        "",
    ]
    if rscc["ablation_beats_full_model"]:
        lines += [
            RSCC_DISCLOSURE_MARKER,
            (
                "Removing every RSCC `w²` term and `δ_g` lowers the maximum "
                f"mixed-chart residual from "
                f"`{_fmt_decimal(rscc['full_rscc']['max_abs_relative_error_percent'], 8)}%` "
                f"to `{_fmt_decimal(rscc['zero_w2_zero_delta_g_ablation']['max_abs_relative_error_percent'], 8)}%` "
                "and lowers the raw diagonal residual sum from "
                f"`{_fmt_decimal(rscc['full_rscc']['raw_diagonal_residual_sum'], 8)}` "
                f"to `{_fmt_decimal(rscc['zero_w2_zero_delta_g_ablation']['raw_diagonal_residual_sum'], 8)}`."
            ),
            "Therefore the detailed RSCC covariance ledger is not selected by its",
            "own retrospective comparison metric. These numbers are diagnostics,",
            "not likelihoods.",
            "",
        ]

    lines += [
        "## W3a — constant substitution in the pixel closure",
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
        "## W3b — selector and port-menu ablation",
        "",
        selector["summary"]["verdict"],
        "",
        "| Vertex configuration | Ports | Producer state | SM Lie type |",
        "|---|---:|---|---|",
    ]
    for row in selector["rows"]:
        if row["sm_lie_type_available"] is True:
            sm = "available; **not uniquely selected**"
        else:
            sm = "unknown (fail closed)"
        lines.append(
            f"| `{row['configuration']}` | {row['port_count']} | "
            f"`{row['status']}` | {sm} |"
        )
    no_go = selector["twelve_port_completion_nonuniqueness"]
    lines += [
        "",
        "For the only executable port-response entry (the 12-port icosahedron),",
        "the exact coefficient classification leaves three viable compact Lie",
        "types: `su(3)+su(2)+u(1)`, `su(2)+su(2)+u(1)^6`, and `u(1)^12`.",
        "The declared finite response construction realizes the first and its",
        "source-binding gate passes on the twelve-port carrier. Laboratory",
        "gauge-current identification remains open, and the bare source reduct",
        "admits both",
        f"`{'` and `'.join(no_go['inequivalent_current_completions'])}`.",
        "",
        "Matter completion is likewise non-unique: the same reduct admits",
        f"`{'` and `'.join(no_go['inequivalent_matter_completions'])}`. The",
        "rank-15 projector is verified, but the sterile-singlet completion is not",
        "excluded. Seven alternative deltahedra lack a compatible producer and are",
        "left unknown rather than treated as negative results.",
        "",
        "## W3c — RSCC zero-`w²` ablation gate",
        "",
        "| Model | Maximum residual | Raw residual sum |",
        "|---|---:|---:|",
        (
            f"| Full RSCC | "
            f"`{_fmt_decimal(rscc['full_rscc']['max_abs_relative_error_percent'], 10)}%` | "
            f"`{_fmt_decimal(rscc['full_rscc']['raw_diagonal_residual_sum'], 10)}` |"
        ),
        (
            f"| Zero-`w²`, zero-`δ_g` | "
            f"`{_fmt_decimal(rscc['zero_w2_zero_delta_g_ablation']['max_abs_relative_error_percent'], 10)}%` | "
            f"`{_fmt_decimal(rscc['zero_w2_zero_delta_g_ablation']['raw_diagonal_residual_sum'], 10)}` |"
        ),
        "",
        f"Negative control beats full model: **{str(rscc['ablation_beats_full_model']).lower()}**.",
        "The checker exits nonzero if this result is absent from the page header.",
        "",
        "## W3d — bounded quark Clebsch assignment audit",
        "",
        "The full six-permutation menu of `(1, 1/3, 3)` over",
        "`(b/τ, s/μ, d/e)` is recomputed through the current five-Yukawa plus",
        "common-light-down approximation. The discrepancy metric is explicitly",
        "retrospective and target-informed: minimize the worst absolute-mass",
        "relative error, then the light-ratio relative error against the FLAG",
        "`N_f=2+1+1` row. It was not preregistered and selects no physical order.",
        "",
        "| b/τ | s/μ | d/e | Worst mass error | FLAG light-ratio error | m_s/m_d | Result |",
        "|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in clebsch["rows"]:
        assignment = row["assignment"]
        metric = row["retrospective_metric"]
        result = (
            "**current assumed; uniquely least discrepant; FLAG rejected**"
            if row["retrospective_unique_least_discrepant"]
            else "FLAG rejected"
        )
        lines.append(
            f"| `{assignment['b_over_tau']}` | `{assignment['s_over_mu']}` | "
            f"`{assignment['d_over_e']}` | "
            f"`{_fmt_decimal(metric['worst_absolute_mass_relative_error'], 8)}` | "
            f"`{_fmt_decimal(metric['light_ratio_absolute_relative_error_vs_flag_nf_2+1+1'], 8)}` | "
            f"`{_fmt_decimal(row['conditional_outputs']['ms_over_md'], 10)}` | "
            f"{result} |"
        )
    identity = clebsch["scale_qualified_rg_identity"]
    flag = clebsch["flag_2024_compare_only"]
    pairing = clebsch["pairing_and_weight_set_provenance"]
    koide = clebsch["koide_and_sqrt_mass_invariant"]
    lines += [
        "",
        "All six assignments miss the conservative FLAG gate. The current",
        "assumed order is uniquely least discrepant only under the retrospective",
        "target-informed metric above. Physical generation order remains open",
        "under `GENERATION_REGISTER_ORDER`; relation existence, independent",
        "coefficient identification, and the third-generation factor are open too.",
        "",
        "The matter receipt separately proves one invariant line for each of",
        "`Q/Sbar/d_c` and `L/Sbar/e_c`, and zero for the forbidden",
        "`Q/S/d_c` control. Those exact conditional facts constrain candidate",
        "channel pairing; they **do not equate the two independent Yukawa",
        "coefficients**. Likewise, F1 measure balance and F2 register",
        "faithfulness are conditional declared algebra, not source-derived",
        "physical laws.",
        "",
        "### Scale-qualified RG identity",
        "",
        f"`{identity['exact_algebraic_statement']}`.",
        "The factor `1/9` is exactly `(1/3)/3`; the common down-lane Yukawa",
        "and same-scale QCD factors cancel. The implementation verifies the",
        "corresponding factor-ratio identity for all six permutations.",
        "",
        f"- Adopted assumed-order output: `{identity['adopted_output_ms_over_md']}`",
        f"- Register-scale right-hand side: `{identity['adopted_register_scale_rhs']}`",
        f"- Relative gap: `{identity['adopted_register_scale_relative_gap']}`",
        f"- Low-scale MCPR-display `(m_mu/m_e)/9`: `{identity['literal_low_scale_mcpr_input_rhs']}`",
        f"- Low-scale substitution relative gap: `{identity['literal_low_scale_relative_gap']}`",
        "",
        "The unqualified low-scale identity is **not exact** in the declared",
        "implementation: flavor-dependent charged-lepton self-running changes",
        "`y_mu/y_e` before `μ_U`. The RG-protected statement therefore requires",
        "the register-scale lepton ratio. Common down-lane and same-scale QCD",
        "transport cancel, but an arbitrary generation-dependent threshold would",
        "break the identity and define a modified lane.",
        "",
        "### FLAG 2024 compare-only rejection of the declared route",
        "",
        "The comparison consumes the two published dimensionless ratios",
        "`m_s/m_ud` and `m_u/m_d` and derives",
        "`m_s/m_d = (m_s/m_ud)(1+m_u/m_d)/2`. The input covariance is",
        "unavailable. The table therefore reports independent propagation and",
        "the maximally positively correlated (`ρ=+1`) case separately. The",
        "retrospective `5σ` rejection gate uses the larger, conservative `ρ=+1`",
        "uncertainty. The prediction preexisted this audit, but the significance",
        "rule was not preregistered.",
        "",
        "| FLAG row | Derived m_s/m_d | Independent σ (gap) | ρ=+1 σ (gap) | Gate |",
        "|---|---:|---:|---:|---|",
    ]
    for row in flag["rows"]:
        independent = row["independent_propagation"]
        positive = row[
            "maximally_positive_correlation_propagation"
        ]
        gate = row["conservative_rejection_gate"]
        lines.append(
            f"| `N_f={row['nf']}` | "
            f"`{_fmt_decimal(row['derived_ms_over_md'], 10)}` | "
            f"`{_fmt_decimal(independent['standard_uncertainty'], 9)}` "
            f"(`{_fmt_decimal(independent['gap_sigma'], 6)}σ`) | "
            f"`{_fmt_decimal(positive['standard_uncertainty'], 9)}` "
            f"(`{_fmt_decimal(positive['gap_sigma'], 6)}σ`) | "
            f"**{'rejected' if gate['rejected'] else 'not rejected'}** |"
        )
    lines += [
        "",
        "The FLAG values enter only the explicitly retrospective discrepancy",
        "metric and rejection gate. They do not enter the lane's arithmetic",
        "solve or select a physical Clebsch order. No OPH theory uncertainty is",
        "supplied or invented. The current conditional lane output is rejected",
        "for both declared `N_f` rows.",
        "",
        "### Koide / square-root-mass invariant provenance",
        "",
        f"- MCPR invariant: `{koide['mcpr_sqrt_mass_invariant']}`",
        f"- Formula: `{koide['mcpr_formula']}` with `chi={koide['mcpr_chi']}`",
        f"- Distance from `2/3`: `{koide['distance_from_two_thirds']}`",
        f"- Exact `2/3`: **{str(koide['mcpr_is_exact_two_thirds']).lower()}**",
        "",
        koide["classification"],
        "",
        "## W-F6 — bounded correction-expression grammar null",
        "",
        expression_null["source_grammar_reconstruction_status"],
        "",
        "The replacement grammar is exact and deliberately shallow. Its atoms are",
        "`P, phi, pi, e, sqrt(P), sqrt(pi), ln(2), P-phi, sqrt(5), 1, 2, 3, 4, 5`.",
        "It includes every atom, unordered atom pairs under `+` and `*`, ordered",
        "distinct atom pairs under `-`, and all ordered atom pairs under `/`.",
        "There is no recursive composition. This produces",
        f"**{expression_null['grammar']['syntactic_expression_count']} syntactic expressions**,",
        f"of which **{expression_null['grammar']['expressions_inside_random_target_interval']}**",
        "lie in the random-target interval. Syntactic aliases are retained; they",
        "cannot increase a yes/no any-hit fraction.",
        "",
        "The calibration uses 2,000 platform-stable SHA-256 counter-mode",
        "pseudorandom targets in the open interval `(0.6,0.8)`. Distance is",
        "relative: `|expression-target|/|target|`.",
        "",
        "| Relative tolerance | Random targets with at least one hit | Hit fraction |",
        "|---:|---:|---:|",
    ]
    for display, _ in EXPRESSION_NULL_TOLERANCES:
        row = expression_null["thresholds"][display]
        lines.append(
            f"| `{display}` | {row['random_target_hit_count']} / "
            f"{row['random_target_count']} | "
            f"`{row['random_target_hit_percent']}%` |"
        )
    lines += [
        "",
        "These are freshly computed rates for the declared 602-expression grammar,",
        "**not a reproduction** of the underspecified 280-expression rates in the",
        "audit memo. They grant no positive evidential weight.",
        "",
        "| Audit-supplied correction target | Nearest expression | Relative distance | 0.20% hit | 0.10% hit | 0.05% hit |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in expression_null["candidate_correction_diagnostics"]:
        hits = row["threshold_hits"]
        lines.append(
            f"| `{row['label']} = {row['audit_supplied_rounded_target']}` | "
            f"`{row['nearest_expression']} = {row['nearest_value']}` | "
            f"`{row['nearest_relative_distance_percent']}%` | "
            f"{str(hits['0.20%']).lower()} | "
            f"{str(hits['0.10%']).lower()} | "
            f"{str(hits['0.05%']).lower()} |"
        )
    lines += [
        "",
        expression_null["claim_entry_gate"]["verdict"],
        "Any future correction constant must ship an exact finite grammar, a",
        "relative tolerance, and the deterministic null-target hit rate before",
        "it can enter a claim. A short-expression near hit is not a derivation.",
        "",
        "## Reproducibility",
        "",
        "| Input | SHA-256 |",
        "|---|---|",
    ]
    for path, digest in sorted(report["input_sha256"].items()):
        lines.append(f"| `{path}` | `{digest}` |")
    lines.append("")
    text = "\n".join(lines)
    enforce_rscc_front_page_disclosure(report, text)
    return text


def enforce_rscc_front_page_disclosure(
    report: Mapping[str, Any],
    scorecard_text: str,
) -> None:
    """Fail if a winning RSCC ablation is hidden below the first W3 section."""

    if not report["rscc_ablation"]["ablation_beats_full_model"]:
        return
    front_page = scorecard_text.split("## W3a", 1)[0]
    require(
        RSCC_DISCLOSURE_MARKER in front_page,
        "RSCC ablation beats the full model but the scorecard header does not disclose it",
    )


def check_scorecard(path: Path, expected: str, report: Mapping[str, Any]) -> None:
    require(path.is_file(), f"null-model scorecard is missing: {path}")
    current = path.read_text(encoding="utf-8")
    enforce_rscc_front_page_disclosure(report, current)
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
                "selector probes, 6 Clebsch permutations all FLAG-rejected, "
                f"{report['expression_grammar_null_model']['random_target_sampler']['target_count']} "
                "W-F6 random targets, RSCC ablation disclosed"
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
