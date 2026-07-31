#!/usr/bin/env python3
"""Independent fail-closed checker for the scalar declared-chart pole receipt.

The checker validates the strict v2 schema and then replays the serialized
structural and interval-arithmetic evidence without importing either
certified producer.  It recomputes source/module digests, exact fixture and
sheet-vector facts, non-certifying diagnostic summaries, partition closure,
centered-form interval arithmetic, cone and endpoint-increment gates,
winding totals, strict interval-Newton inclusion margins, residue inversion,
and every precision-nesting comparison.  It deliberately does not
independently re-evaluate the loop functions; exact producer ``--verify``
remains a separate required trust-boundary check.

A pass certifies only a scalar simple pole on the recorded algebraic chart.
It does not certify the chart identity as a physical/unique second sheet,
a sign bridge to the separately written theorem convention, full matrix
rank-n-1 Laurent data, issue #593's combined row, BMHV restoration, a
physical-current pole, OPH-native provenance, units, or unitarity.
"""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from math import isqrt
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator
from mpmath import iv

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_PATH = ROOT / "outputs" / "certified_second_sheet_poles.json"
ZERO_EXCLUSION_PATH = ROOT / "outputs" / "certified_wz_contours.json"
VECTOR_PATH = ROOT / "outputs" / "fj_direct_vector_blocks.json"
SCHEMA_PATH = ROOT / "schemas" / "certified_second_sheet_poles_v2.schema.json"

PASS_STATUS = "SCALAR_DECLARED_CHART_POLE_CHECK_PASS"
SUCCESS_STATUS = "SCALAR_POLE_CERTIFIED_ON_DECLARED_ALGEBRAIC_CHART"

PRECISIONS = (128, 192, 256)
TREE_MASSES = {"W": Fraction(1, 9), "Z": Fraction(25, 144)}
WINDOWS = {
    "W": (Fraction(9, 100), Fraction(1, 9)),
    "Z": (Fraction(4, 25), Fraction(25, 144)),
}
POLE_BOXES = {
    "W": (
        (Fraction(1119, 10000), Fraction(14, 125)),
        (Fraction(-9, 10000), Fraction(-7, 10000)),
    ),
    "Z": (
        (Fraction(1742, 10000), Fraction(1746, 10000)),
        (Fraction(-15, 10000), Fraction(-11, 10000)),
    ),
}
CORRIDOR_BOXES = {
    "W": (
        (Fraction(9, 100), Fraction(113, 1000)),
        (Fraction(-1, 1000), Fraction(-1, 10**9)),
    ),
    "Z": (
        (Fraction(16, 100), Fraction(1746, 10000)),
        (Fraction(-16, 10000), Fraction(-1, 10**9)),
    ),
}
NEWTON_SEEDS = {
    "W": (
        (Fraction(111953, 10**6), Fraction(111955, 10**6)),
        (Fraction(-803, 10**6), Fraction(-802, 10**6)),
    ),
    "Z": (
        (Fraction(1744045, 10**7), Fraction(1744055, 10**7)),
        (Fraction(-12735, 10**7), Fraction(-12725, 10**7)),
    ),
}
EXPECTED_FIXTURE = {
    "g1": "1/4",
    "g2": "1/3",
    "lam": "1/8",
    "mfd1": "1/60",
    "mfd2": "1/25",
    "mfd3": "1/10",
    "mfe1": "1/80",
    "mfe2": "1/30",
    "mfe3": "1/15",
    "mfu1": "1/50",
    "mfu2": "1/20",
    "mfu3": "1/5",
    "mu2": "1/2",
    "mu_ren2": "1",
    "v": "2",
    "xi": "1",
}
EXPECTED_CKM = {
    f"V{i}{j}": "1" if i == j else "0"
    for i in (1, 2, 3)
    for j in (1, 2, 3)
}
EXPECTED_CORRECTIONS = {
    "W": "(312000*p2 + 90163)/1944000",
    "Z": "(292020000*p2 + 79268959)/1944000000",
}
EXPECTED_PAIRS = {
    "W": {
        (Fraction(0), Fraction(1, 9)),
        (Fraction(0), Fraction(1, 6400)),
        (Fraction(0), Fraction(1, 900)),
        (Fraction(0), Fraction(1, 225)),
        (Fraction(25, 144), Fraction(1, 9)),
        (Fraction(1, 9), Fraction(1)),
        (Fraction(1, 9), Fraction(25, 144)),
        (Fraction(1, 9), Fraction(0)),
        (Fraction(1, 2500), Fraction(1, 3600)),
        (Fraction(1, 400), Fraction(1, 625)),
        (Fraction(1, 25), Fraction(1, 100)),
    },
    "Z": {
        (Fraction(0), Fraction(0)),
        (Fraction(1, 6400), Fraction(1, 6400)),
        (Fraction(1, 3600), Fraction(1, 3600)),
        (Fraction(1, 2500), Fraction(1, 2500)),
        (Fraction(1, 900), Fraction(1, 900)),
        (Fraction(1, 625), Fraction(1, 625)),
        (Fraction(1, 400), Fraction(1, 400)),
        (Fraction(1, 225), Fraction(1, 225)),
        (Fraction(1, 100), Fraction(1, 100)),
        (Fraction(1, 25), Fraction(1, 25)),
        (Fraction(1, 9), Fraction(1, 9)),
        (Fraction(25, 144), Fraction(1)),
    },
}
EXPECTED_DENOMINATOR_COUNTS = {"W": 67, "Z": 40}
EXPECTED_DENOMINATOR_IDENTITY_DIGESTS = {
    "W": "sha256:fa9fb5908600498006f162972af47c161549968d7042360dce07972f90a7f2f9",
    "Z": "sha256:532c45f89ec092a287778796fdcd1b1d87e5b18a1f83149734f4c2f3589fd6cc",
}
EXPECTED_A0_MASSES = {
    Fraction(0),
    Fraction(1),
    Fraction(1, 100),
    Fraction(1, 225),
    Fraction(1, 25),
    Fraction(1, 2500),
    Fraction(1, 3600),
    Fraction(1, 400),
    Fraction(1, 625),
    Fraction(1, 6400),
    Fraction(1, 9),
    Fraction(1, 900),
    Fraction(25, 144),
}
EXPECTED_SHEET_STATEMENT = (
    "the pole boxes lie strictly in the open lower half plane on a "
    "declared algebraic chart with an explicit per-channel sheet "
    "vector; in the W chart W-gamma has the declared action principal, "
    "meaning no added chart-relative correction rather than a "
    "certified physical-sheet assignment; "
    "the one-mass base chart is mass-exchange symmetric and every "
    "correction is relative to its recorded base chart; "
    "finite-delta probes are non-certifying diagnostics, and no "
    "independent continuation identity or physical/unique resonance "
    "sheet is certified"
)
CERTIFIED_READING = (
    "the scalar transverse inverse propagator has exactly one simple zero "
    "in the declared lower-half pole box on the declared algebraic chart; "
    "this does not certify a unique or physical resonance sheet or full "
    "matrix Laurent data"
)

RI = tuple[Fraction, Fraction]
CI = tuple[RI, RI]


class ReceiptValidationError(ValueError):
    """Raised when a receipt fails any schema, binding, or replay gate."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ReceiptValidationError(message)


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"


def sha256_bytes(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def parse_fraction(value: Any, label: str) -> Fraction:
    require(isinstance(value, str) and value, f"{label} is not a rational string")
    try:
        parsed = Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ReceiptValidationError(f"{label} is not rational: {error}") from error
    require(str(parsed) == value, f"{label} is not canonically serialized")
    return parsed


def parse_diagnostic_number(value: Any, label: str) -> Fraction:
    """Parse a finite decimal/scientific diagnostic without promoting it."""

    require(isinstance(value, str) and value, f"{label} is not a numeric string")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ReceiptValidationError(
            f"{label} is not a finite decimal diagnostic: {error}"
        ) from error


def parse_ri(value: Any, label: str) -> RI:
    require(
        isinstance(value, dict) and set(value) == {"lo", "hi"},
        f"{label} is not a two-endpoint real interval",
    )
    result = (
        parse_fraction(value["lo"], f"{label}.lo"),
        parse_fraction(value["hi"], f"{label}.hi"),
    )
    require(result[0] <= result[1], f"{label} endpoints are reversed")
    return result


def parse_rational_pair(value: Any, label: str) -> RI:
    require(
        isinstance(value, list) and len(value) == 2,
        f"{label} is not a rational endpoint pair",
    )
    result = (
        parse_fraction(value[0], f"{label}[0]"),
        parse_fraction(value[1], f"{label}[1]"),
    )
    require(result[0] <= result[1], f"{label} endpoints are reversed")
    return result


def parse_rational_box(value: Any, label: str) -> CI:
    require(
        isinstance(value, dict) and set(value) == {"re", "im"},
        f"{label} is not a rational box",
    )
    return (
        parse_rational_pair(value["re"], f"{label}.re"),
        parse_rational_pair(value["im"], f"{label}.im"),
    )


def parse_ci(value: Any, label: str) -> CI:
    require(isinstance(value, dict), f"{label} is not a complex interval")
    require(set(value) == {"re", "im"}, f"{label} fields are not exact")
    return parse_ri(value["re"], f"{label}.re"), parse_ri(
        value["im"], f"{label}.im"
    )


def r_add(a: RI, b: RI) -> RI:
    return a[0] + b[0], a[1] + b[1]


def r_sub(a: RI, b: RI) -> RI:
    return a[0] - b[1], a[1] - b[0]


def r_mul(a: RI, b: RI) -> RI:
    products = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return min(products), max(products)


def r_square(a: RI) -> RI:
    high = max(abs(a[0]), abs(a[1])) ** 2
    low = Fraction(0) if a[0] <= 0 <= a[1] else min(abs(a[0]), abs(a[1])) ** 2
    return low, high


def r_div_positive(a: RI, denominator: RI) -> RI:
    require(denominator[0] > 0, "interval division denominator is not positive")
    return r_mul(a, (1 / denominator[1], 1 / denominator[0]))


def c_add(a: CI, b: CI) -> CI:
    return r_add(a[0], b[0]), r_add(a[1], b[1])


def c_sub(a: CI, b: CI) -> CI:
    return r_sub(a[0], b[0]), r_sub(a[1], b[1])


def c_mul(a: CI, b: CI) -> CI:
    return (
        r_sub(r_mul(a[0], b[0]), r_mul(a[1], b[1])),
        r_add(r_mul(a[0], b[1]), r_mul(a[1], b[0])),
    )


def c_pow(value: CI, exponent: int) -> CI:
    require(exponent >= 0, "negative polynomial power")
    result: CI = ((Fraction(1), Fraction(1)), (Fraction(0), Fraction(0)))
    base = value
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = c_mul(result, base)
        base = c_mul(base, base)
        remaining >>= 1
    return result


def c_conj(value: CI) -> CI:
    return value[0], (-value[1][1], -value[1][0])


def c_div(a: CI, b: CI) -> CI:
    denominator = r_add(r_square(b[0]), r_square(b[1]))
    numerator = c_mul(a, c_conj(b))
    return (
        r_div_positive(numerator[0], denominator),
        r_div_positive(numerator[1], denominator),
    )


def ri_contains(outer: RI, inner: RI) -> bool:
    return outer[0] <= inner[0] and inner[1] <= outer[1]


def ci_contains(outer: CI, inner: CI) -> bool:
    return ri_contains(outer[0], inner[0]) and ri_contains(outer[1], inner[1])


def ci_nested(outer: CI, inner: CI) -> bool:
    return ci_contains(outer, inner)


def ci_contains_zero(value: CI) -> bool:
    return value[0][0] <= 0 <= value[0][1] and value[1][0] <= 0 <= value[1][1]


def ci_zero_exclusion_abs2_lower(value: CI) -> Fraction:
    def distance(interval: RI) -> Fraction:
        if interval[0] <= 0 <= interval[1]:
            return Fraction(0)
        return min(abs(interval[0]), abs(interval[1]))

    return distance(value[0]) ** 2 + distance(value[1]) ** 2


def excludes_principal_cut(value: CI) -> bool:
    if ci_contains_zero(value):
        return False
    if value[0][0] > 0:
        return True
    return value[1][0] > 0 or value[1][1] < 0


def c_intersection(a: CI, b: CI) -> CI | None:
    result = (
        (max(a[0][0], b[0][0]), min(a[0][1], b[0][1])),
        (max(a[1][0], b[1][0]), min(a[1][1], b[1][1])),
    )
    if result[0][0] > result[0][1] or result[1][0] > result[1][1]:
        return None
    return result


def _iv_fraction(value: Fraction) -> Any:
    return iv.mpf(value.numerator) / iv.mpf(value.denominator)


def _mpf_tuple_fraction(raw: tuple[int, int, int, int]) -> Fraction:
    sign, mantissa, exponent, _bits = raw
    value = Fraction(mantissa)
    value = value * 2**exponent if exponent >= 0 else value / 2 ** (-exponent)
    return -value if sign else value


def _ri_to_iv(value: RI) -> Any:
    return iv.mpf(
        [_iv_fraction(value[0]).a, _iv_fraction(value[1]).b]
    )


def _iv_to_ri(value: Any) -> RI:
    return (
        _mpf_tuple_fraction(value._mpi_[0]),
        _mpf_tuple_fraction(value._mpi_[1]),
    )


def directed_arg(value: CI, precision: int) -> RI:
    require(excludes_principal_cut(value), "argument rectangle meets its cut")
    iv.prec = precision
    re = iv.mpf(
        [_iv_fraction(value[0][0]).a, _iv_fraction(value[0][1]).b]
    )
    im = iv.mpf(
        [_iv_fraction(value[1][0]).a, _iv_fraction(value[1][1]).b]
    )
    angle = iv.atan2(im, re)
    return (
        _mpf_tuple_fraction(angle._mpi_[0]),
        _mpf_tuple_fraction(angle._mpi_[1]),
    )


def c_real(value: Fraction) -> CI:
    return ((value, value), (Fraction(0), Fraction(0)))


def c_neg(value: CI) -> CI:
    return (
        (-value[0][1], -value[0][0]),
        (-value[1][1], -value[1][0]),
    )


def principal_sqrt(value: CI, precision: int) -> CI:
    """Independent interval replay of the declared principal square root."""

    angle = directed_arg(value, precision)
    iv.prec = precision
    radius = iv.sqrt(
        iv.sqrt(_ri_to_iv(r_add(r_square(value[0]), r_square(value[1]))))
    )
    half_angle = _ri_to_iv(angle) / iv.mpf(2)
    return (
        _iv_to_ri(radius * iv.cos(half_angle)),
        _iv_to_ri(radius * iv.sin(half_angle)),
    )


def expected_chart_correction(
    pair: tuple[Fraction, Fraction],
    segment_box: CI,
    precision: int,
) -> tuple[CI, CI]:
    """Replay the mass-symmetric chart correction and its derivative."""

    iv.prec = precision
    two_pi = _iv_to_ri(iv.mpf(2) * iv.pi)
    two_pi_i: CI = ((Fraction(0), Fraction(0)), two_pi)
    m1, m2 = pair
    zero: CI = c_real(Fraction(0))
    if m1 == 0 and m2 == 0:
        return zero, zero
    if m1 == 0 or m2 == 0:
        mass = m2 if m1 == 0 else m1
        correction = c_mul(
            two_pi_i,
            c_sub(c_real(Fraction(1)), c_div(c_real(mass), segment_box)),
        )
        slope = c_div(
            c_mul(two_pi_i, c_real(mass)),
            c_mul(segment_box, segment_box),
        )
        return correction, slope
    b = c_neg(c_add(segment_box, c_real(m1 - m2)))
    discriminant = c_sub(
        c_mul(b, b),
        c_mul(c_real(4 * m1), segment_box),
    )
    square_root = principal_sqrt(discriminant, precision)
    correction = c_div(c_mul(two_pi_i, square_root), segment_box)
    square_root_slope = c_div(
        c_sub(segment_box, c_real(m1 + m2)),
        square_root,
    )
    slope = c_div(
        c_mul(
            two_pi_i,
            c_sub(
                c_mul(square_root_slope, segment_box),
                square_root,
            ),
        ),
        c_mul(segment_box, segment_box),
    )
    return correction, slope


def rational_sqrt(value: Fraction) -> Fraction:
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    require(
        numerator * numerator == value.numerator
        and denominator * denominator == value.denominator,
        f"mass {value} has no exact rational square root",
    )
    return Fraction(numerator, denominator)


def expected_base_and_correction(m1: Fraction, m2: Fraction) -> tuple[str, str]:
    if m1 == 0 and m2 == 0:
        return "both_massless_upper_half_continued", "0"
    if m1 == 0 or m2 == 0:
        mass = m2 if m1 == 0 else m1
        return (
            "one_mass_symmetric_lower_principal",
            f"2*pi*i*(1-({mass})/s)",
        )
    return "two_mass_principal_root_chart", "2*pi*i*sqrt(lambda(s))/s"


def check_schema(receipt: dict[str, Any], schema: dict[str, Any]) -> None:
    errors = sorted(
        Draft202012Validator(schema).iter_errors(receipt),
        key=lambda error: error.json_path,
    )
    if errors:
        error = errors[0]
        raise ReceiptValidationError(
            f"schema {error.json_path}: {error.message}"
        )


def check_bindings(receipt: dict[str, Any]) -> None:
    recomputed = sha256_bytes(
        canonical_json(
            {key: value for key, value in receipt.items() if key != "receipt_sha256"}
        ).encode("utf-8")
    )
    require(receipt["receipt_sha256"] == recomputed, "receipt digest mismatch")

    vector_raw = VECTOR_PATH.read_bytes()
    vector = json.loads(vector_raw)
    source = receipt["source_subject"]
    require(
        source
        == {
            "relative_path": "outputs/fj_direct_vector_blocks.json",
            "artifact_role": "direct_FJ_vector_blocks_input",
            "schema": "fj_direct_vector_blocks.v1",
            "target": "FJ_DIRECT_1",
            "units": (
                "loop measure i/(16 pi^2) stripped; Delta is the single "
                "1/eps pole unit"
            ),
            "bytes": len(vector_raw),
            "sha256": sha256_bytes(vector_raw),
        },
        "source subject does not bind the exact vector artifact",
    )
    require(vector.get("schema") == source["schema"], "source schema drift")
    require(vector.get("target") == source["target"], "source target drift")
    require(vector.get("units") == source["units"], "source units drift")

    pins = receipt["pins"]
    expected_pins = {
        "vector_blocks_sha256": sha256_bytes(vector_raw),
        "diagnostic_producer_sha256": sha256_bytes(
            (ROOT / "producers" / "wz_pole_receipts.py").read_bytes()
        ),
        "interval_arithmetic_module_sha256": sha256_bytes(
            (ROOT / "producers" / "complex_interval.py").read_bytes()
        ),
        "principal_zero_exclusion_producer_sha256": sha256_bytes(
            (ROOT / "producers" / "certified_wz_contours.py").read_bytes()
        ),
        "producer_module_sha256": sha256_bytes(
            (ROOT / "producers" / "certified_second_sheet_poles.py").read_bytes()
        ),
        "checker_module_sha256": sha256_bytes(Path(__file__).read_bytes()),
        "principal_zero_exclusion_receipt_sha256": sha256_bytes(
            ZERO_EXCLUSION_PATH.read_bytes()
        ),
    }
    require(pins == expected_pins, "module or principal-receipt pin drift")
    require(receipt["fixture"] == EXPECTED_FIXTURE, "fixture drift")
    require(receipt["ckm_fixture"] == EXPECTED_CKM, "identity CKM fixture drift")
    require(
        receipt["dimensional_prefactor_finite_correction"]
        == EXPECTED_CORRECTIONS,
        "finite dimensional-prefactor correction drift",
    )


def check_scope_and_gates(receipt: dict[str, Any]) -> None:
    require(
        receipt["sheet_statement"] == EXPECTED_SHEET_STATEMENT,
        "sheet statement drift",
    )
    require(
        receipt["claim_scope"]
        == {
            "certified_object": (
                "scalar_transverse_inverse_propagator_on_declared_"
                "algebraic_chart"
            ),
            "declared_chart_only": True,
            "continuation_identity_independently_certified": False,
            "standard_second_sheet_identification_certified": False,
            "full_matrix_rank_n_minus_1_laurent_certified": False,
            "issue_593_pole_laurent_row_discharged": False,
            "finite_delta_diagnostics_are_noncertifying": True,
            "independent_numerical_replay_certified": False,
            "engine_inverse_propagator_convention": (
                "G(s)=s-m_tree^2-Pi_engine(s)"
            ),
            "theorem_self_energy_sign_bridge_certified": False,
            "issue_593_precision_ladder_row_discharged": False,
            "issue_593_independent_third_verifier_row_discharged": False,
            "issue_593_full_acceptance_satisfied": False,
        },
        "claim scope is not the exact fail-closed scalar statement",
    )
    require(
        receipt["serialized_gates"]
        == {
            "precisions_bits": list(PRECISIONS),
            "initial_segments_per_edge": 8,
            "max_subdivision_depth": 12,
            "holomorphy_method": "exact_rational_cut_exclusion_lower_half_plane",
            "boundary_method": "centered_form_segment_enclosure",
            "arg_width_gate": "157/100",
            "winding_tolerance": "157/100",
            "newton_steps": 4,
            "probe_residual_gate": "1/1000",
            "newton_strict_interior_required": True,
            "interval_encoding": "exact_dyadic_rational_endpoints",
        },
        "serialized gate contract drift",
    )


def check_continuation(receipt: dict[str, Any], name: str) -> set[tuple[Fraction, Fraction]]:
    continuation = receipt["declared_continuation"][name]
    window = tuple(
        parse_fraction(value, f"{name}.window") for value in continuation["window"]
    )
    require(window == WINDOWS[name], f"{name} declared window drift")
    require(
        parse_fraction(continuation["tree_mass_sq_context"], f"{name}.tree")
        == TREE_MASSES[name],
        f"{name} tree-mass context drift",
    )
    require(
        continuation["chart_definition"]
        == {
            "identity_status": "declared_not_independently_certified",
            "correction_reference": "relative_to_declared_base_chart",
            "one_mass_mass_exchange_symmetric": True,
            "window_is_declared_input": True,
        },
        f"{name} chart definition drift",
    )
    require(window[1] < POLE_BOXES[name][0][0], f"{name} window meets pole box")

    pairs: set[tuple[Fraction, Fraction]] = set()
    crossed: set[tuple[Fraction, Fraction]] = set()
    expected_vector = []
    previous_pair: tuple[Fraction, Fraction] | None = None
    for index, channel in enumerate(continuation["channels"]):
        m1 = parse_fraction(channel["m1"], f"{name}.channels[{index}].m1")
        m2 = parse_fraction(channel["m2"], f"{name}.channels[{index}].m2")
        pair = (m1, m2)
        require(pair not in pairs, f"{name} repeats channel {pair}")
        require(
            previous_pair is None or pair > previous_pair,
            f"{name} channels are not in canonical pair order",
        )
        previous_pair = pair
        pairs.add(pair)
        threshold = (rational_sqrt(m1) + rational_sqrt(m2)) ** 2
        require(
            parse_fraction(channel["threshold"], f"{name}.threshold")
            == threshold,
            f"{name} threshold does not replay for {pair}",
        )
        require(
            not (window[0] < threshold < window[1]),
            f"{name} threshold lies inside declared window",
        )
        action = "crossed" if threshold <= window[0] else "principal"
        require(channel["sheet_action"] == action, f"{name} sheet action drift")
        base, correction = expected_base_and_correction(m1, m2)
        require(channel["base_chart"] == base, f"{name} base chart drift")
        require(
            channel["crossing_correction"] == correction,
            f"{name} crossing correction drift for {pair}",
        )
        require(
            channel["applied_correction"]
            == (correction if action == "crossed" else "0"),
            f"{name} applied correction drift for {pair}",
        )
        require(
            channel["correction_reference"]
            == "relative_to_declared_base_chart",
            f"{name} correction reference drift",
        )
        if action == "crossed":
            crossed.add(pair)
        expected_vector.append(
            {"m1": str(m1), "m2": str(m2), "sheet_action": action}
        )
    require(pairs == EXPECTED_PAIRS[name], f"{name} channel census drift")
    require(
        max(
            (rational_sqrt(m1) + rational_sqrt(m2)) ** 2
            for m1, m2 in crossed
        )
        == window[0],
        f"{name} crossed frontier does not equal window lower endpoint",
    )
    require(
        continuation["sheet_vector"] == expected_vector,
        f"{name} exact sheet vector drift",
    )

    diagnostics = continuation["consistency_diagnostics"]
    require(
        diagnostics["role"] == "non_certifying_finite_delta_diagnostic"
        and diagnostics["gates_scalar_certificate"] is False,
        f"{name} finite-delta diagnostics were promoted",
    )
    gate = parse_fraction(diagnostics["residual_gate"], f"{name}.diagnostic_gate")
    require(gate == Fraction(1, 1000), f"{name} diagnostic gate drift")
    expected_points = {
        window[0] + (window[1] - window[0]) * Fraction(1, 3),
        window[0] + (window[1] - window[0]) * Fraction(2, 3),
    }
    seen_rows = set()
    aggregate = True
    for index, row in enumerate(diagnostics["rows"]):
        pair = (
            parse_fraction(row["m1"], f"{name}.diagnostics[{index}].m1"),
            parse_fraction(row["m2"], f"{name}.diagnostics[{index}].m2"),
        )
        point = parse_fraction(row["x"], f"{name}.diagnostics[{index}].x")
        require(pair in crossed and point in expected_points, f"{name} stray diagnostic")
        require((pair, point) not in seen_rows, f"{name} duplicate diagnostic")
        seen_rows.add((pair, point))
        require(
            [parse_fraction(value, f"{name}.delta") for value in row["deltas"]]
            == [Fraction(1, 10**6), Fraction(1, 10**8)],
            f"{name} diagnostic delta ladder drift",
        )
        residuals = [
            parse_diagnostic_number(value, f"{name}.diagnostic residual")
            for value in row["residuals_coarse_to_fine"]
        ]
        within = residuals[-1] < gate and residuals[-1] <= residuals[0]
        require(
            row["within_diagnostic_gate"] is within,
            f"{name} diagnostic aggregate flag does not replay",
        )
        expected_kind = "vanishing_jump" if pair == (0, 0) else "declared_addition_match"
        require(row["kind"] == expected_kind, f"{name} diagnostic kind drift")
        aggregate = aggregate and within
    require(
        seen_rows == {(pair, point) for pair in crossed for point in expected_points},
        f"{name} diagnostic census is incomplete",
    )
    require(
        diagnostics["all_within_diagnostic_gate"] is aggregate,
        f"{name} diagnostic summary does not replay",
    )
    return crossed


def check_holomorphy(
    certificate: dict[str, Any],
    box: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
    expected_pairs: set[tuple[Fraction, Fraction]],
    crossed: set[tuple[Fraction, Fraction]],
    expected_denominator_count: int,
    particle: str,
    label: str,
) -> tuple[bool, dict[str, CI]]:
    (re_lo, re_hi), (im_lo, im_hi) = box
    require(
        certificate["method"]
        == "exact_rational_cut_exclusion_lower_half_plane",
        f"{label} holomorphy method drift",
    )
    base = certificate["base_facts"]
    require(
        base["box_in_open_lower_half_plane"] is (im_hi < 0)
        and base["box_real_part_positive"] is (re_lo > 0),
        f"{label} base holomorphy facts drift",
    )
    box_interval: CI = ((re_lo, re_hi), (im_lo, im_hi))
    denominator_enclosures: dict[str, CI] = {}
    identity_payload: list[dict[str, Any]] = []
    previous_identifier = ""
    witnesses = certificate["coefficient_denominator_witnesses"]
    require(
        len(witnesses) == expected_denominator_count,
        f"{label} coefficient-denominator census drift",
    )
    for index, witness in enumerate(witnesses):
        coefficients = [
            (
                int(power),
                parse_fraction(value, f"{label}.denominator[{index}].coefficient"),
            )
            for power, value in witness["coefficients"]
        ]
        require(
            coefficients
            and all(power >= 0 for power, _value in coefficients)
            and [power for power, _value in coefficients]
            == sorted({power for power, _value in coefficients}),
            f"{label} denominator polynomial is not canonical",
        )
        require(
            witness["coefficients"]
            == [[power, str(value)] for power, value in coefficients],
            f"{label} denominator coefficients are not canonical",
        )
        identifier = sha256_bytes(
            canonical_json(
                [[power, str(value)] for power, value in coefficients]
            ).encode("utf-8")
        )
        require(
            witness["denominator_id"] == identifier
            and identifier not in denominator_enclosures,
            f"{label} coefficient-denominator id drift",
        )
        require(
            not previous_identifier or identifier > previous_identifier,
            f"{label} denominator witnesses are not strictly id-sorted",
        )
        previous_identifier = identifier
        identity_payload.append(
            {
                "denominator_id": identifier,
                "coefficients": [
                    [power, str(value)] for power, value in coefficients
                ],
            }
        )
        evaluated: CI = (
            (Fraction(0), Fraction(0)),
            (Fraction(0), Fraction(0)),
        )
        for power, coefficient in coefficients:
            factor: CI = (
                (coefficient, coefficient),
                (Fraction(0), Fraction(0)),
            )
            evaluated = c_add(evaluated, c_mul(factor, c_pow(box_interval, power)))
        enclosure = parse_ci(witness["enclosure"], f"{label}.denominator.enclosure")
        require(
            ci_contains(enclosure, evaluated),
            f"{label} coefficient-denominator enclosure does not replay",
        )
        zero_slack = ci_zero_exclusion_abs2_lower(enclosure)
        require(
            witness["excludes_zero"] is (zero_slack > 0)
            and parse_fraction(
                witness["zero_exclusion_abs2_lower"],
                f"{label}.denominator.zero_slack",
            )
            == zero_slack,
            f"{label} coefficient-denominator zero witness drift",
        )
        denominator_enclosures[identifier] = enclosure
    require(
        sha256_bytes(canonical_json(identity_payload).encode("utf-8"))
        == EXPECTED_DENOMINATOR_IDENTITY_DIGESTS[particle],
        f"{label} denominator identity set is not frozen-source bound",
    )
    denominators_ok = bool(
        denominator_enclosures
        and all(not ci_contains_zero(value) for value in denominator_enclosures.values())
    )
    require(
        certificate["coefficient_denominators_exclude_zero_on_box"]
        is denominators_ok,
        f"{label} coefficient-denominator aggregate drift",
    )

    seen = set()
    all_holomorphic = im_hi < 0 and re_lo > 0
    for row in certificate["loop_charts"]:
        require(row.get("kind") != "denominator", f"{label} denominator failure")
        pair = (
            parse_fraction(row["m1"], f"{label}.m1"),
            parse_fraction(row["m2"], f"{label}.m2"),
        )
        require(pair not in seen, f"{label} repeats loop chart")
        seen.add(pair)
        if pair == (0, 0):
            expected = im_hi < 0
            require(
                row["chart"] == "both_masses_zero"
                and row["certificate"]
                == "log(s/mu2) only; box in open lower half plane",
                f"{label} massless chart documentation drift",
            )
        elif 0 in pair:
            expected = im_hi < 0 and re_lo > 0
            base_chart, _formula = expected_base_and_correction(*pair)
            sign = "positive"
            require(
                row["chart"] == base_chart
                and row["root_chart_log_argument_imaginary_sign"] == sign
                and row["certificate"]
                == (
                    "declared mass-exchange-symmetric one-mass algebraic "
                    f"base chart; root-chart arguments have strictly {sign} "
                    "imaginary part and chart-relative rational corrections "
                    "are holomorphic because the box excludes s=0"
                ),
                f"{label} mass-symmetric one-mass chart drift",
            )
        else:
            symmetry_line = pair[0] + pair[1]
            expected = (
                (symmetry_line < re_lo or symmetry_line > re_hi)
                and im_hi < 0
                and re_lo > 0
            )
            require(
                parse_fraction(
                    row["discriminant_symmetry_line"],
                    f"{label}.symmetry_line",
                )
                == symmetry_line
                and row[
                    "discriminant_symmetry_line_outside_box_re_range"
                ]
                is (symmetry_line < re_lo or symmetry_line > re_hi),
                f"{label} discriminant-symmetry witness drift",
            )
            require(
                row["chart"] == "two_massive"
                and row["certificate"]
                == (
                    "discriminant imaginary part nonvanishing off the "
                    "discriminant-symmetry line Re(s)=m1+m2 (not the "
                    "physical threshold); Feynman roots never real"
                ),
                f"{label} two-mass chart documentation drift",
            )
        require(row["holomorphic"] is expected, f"{label} chart verdict drift")
        require(
            ("chart_correction_holomorphic" in row) is (pair in crossed),
            f"{label} correction-holomorphy field census drift",
        )
        if pair in crossed:
            require(
                row["chart_correction_holomorphic"] is expected,
                f"{label} correction holomorphy drift",
            )
        all_holomorphic = all_holomorphic and expected
    require(seen == expected_pairs, f"{label} loop-chart census drift")
    aggregate_holomorphic = bool(all_holomorphic and denominators_ok)
    require(
        certificate["holomorphic_on_box"] is aggregate_holomorphic,
        f"{label} holomorphy summary does not replay",
    )
    return aggregate_holomorphic, denominator_enclosures


def _is_power_of_two(value: int) -> bool:
    return value > 0 and value & (value - 1) == 0


def segment_edge_and_depth(
    start: tuple[Fraction, Fraction],
    end: tuple[Fraction, Fraction],
    pole_box: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
) -> tuple[int, int] | None:
    """Replay a CCW dyadic child of the eight-cell-per-edge boundary."""

    (re_lo, re_hi), (im_lo, im_hi) = pole_box
    if start[1] == end[1] == im_lo and re_lo <= start[0] < end[0] <= re_hi:
        edge = 0
        edge_length = re_hi - re_lo
        start_offset = start[0] - re_lo
        segment_length = end[0] - start[0]
    elif start[0] == end[0] == re_hi and im_lo <= start[1] < end[1] <= im_hi:
        edge = 1
        edge_length = im_hi - im_lo
        start_offset = start[1] - im_lo
        segment_length = end[1] - start[1]
    elif start[1] == end[1] == im_hi and re_lo <= end[0] < start[0] <= re_hi:
        edge = 2
        edge_length = re_hi - re_lo
        start_offset = re_hi - start[0]
        segment_length = start[0] - end[0]
    elif start[0] == end[0] == re_lo and im_lo <= end[1] < start[1] <= im_hi:
        edge = 3
        edge_length = im_hi - im_lo
        start_offset = im_hi - start[1]
        segment_length = start[1] - end[1]
    else:
        return None
    ratio = edge_length / (8 * segment_length)
    if ratio.denominator != 1 or not _is_power_of_two(ratio.numerator):
        return None
    if (start_offset / segment_length).denominator != 1:
        return None
    return edge, ratio.numerator.bit_length() - 1


def check_winding(
    winding: dict[str, Any],
    pole_box: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
    precision: int,
    label: str,
) -> tuple[bool, list[dict[str, Any]], RI, list[str]]:
    gate = Fraction(157, 100)
    tolerance = Fraction(157, 100)
    require(
        winding["method"] == "centered_form_segment_enclosure"
        and winding["partition"]
        == ("adaptive" if precision == 128 else "replayed_base_partition")
        and winding["reason"] is None,
        f"{label} winding method/partition metadata drift",
    )
    records = winding["segment_evidence"]
    require(records and len(records) == winding["segments"], f"{label} empty partition")

    coordinates = []
    parsed_records: list[dict[str, Any]] = []
    increment_sum: RI = (Fraction(0), Fraction(0))
    signed_area = Fraction(0)
    previous_end = None
    first_start = None
    previous_end_value: CI | None = None
    first_start_value: CI | None = None
    seen_segment_ids: set[str] = set()
    edge_groups: list[list[tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]]] = [
        [],
        [],
        [],
        [],
    ]
    previous_edge = -1
    actual_max_depth = 0
    (re_lo, re_hi), (im_lo, im_hi) = pole_box
    for expected_index, row in enumerate(records):
        require(
            row["segment_id"].startswith("sha256:"),
            f"{label} segment id is malformed",
        )
        require(
            row["segment_id"] not in seen_segment_ids,
            f"{label} repeats a segment id",
        )
        seen_segment_ids.add(row["segment_id"])
        start = (
            parse_fraction(row["start"][0], f"{label}.start.re"),
            parse_fraction(row["start"][1], f"{label}.start.im"),
        )
        end = (
            parse_fraction(row["end"][0], f"{label}.end.re"),
            parse_fraction(row["end"][1], f"{label}.end.im"),
        )
        require(start != end, f"{label} has zero-length segment")
        classification = segment_edge_and_depth(start, end, pole_box)
        require(
            classification is not None,
            f"{label} segment is not a CCW dyadic boundary child",
        )
        edge, depth = classification
        require(depth <= 12, f"{label} segment exceeds subdivision depth cap")
        require(
            edge >= previous_edge,
            f"{label} boundary repeats an edge after a later edge",
        )
        previous_edge = edge
        actual_max_depth = max(actual_max_depth, depth)
        edge_groups[edge].append((start, end))
        if previous_end is not None:
            require(start == previous_end, f"{label} partition is not chained")
        else:
            first_start = start
        previous_end = end
        signed_area += start[0] * end[1] - end[0] * start[1]
        center = (
            parse_fraction(row["midpoint"][0], f"{label}.center.re"),
            parse_fraction(row["midpoint"][1], f"{label}.center.im"),
        )
        require(
            center == ((start[0] + end[0]) / 2, (start[1] + end[1]) / 2),
            f"{label} midpoint is not exact",
        )
        expected_offset: CI = (
            (
                min(start[0], end[0]) - center[0],
                max(start[0], end[0]) - center[0],
            ),
            (
                min(start[1], end[1]) - center[1],
                max(start[1], end[1]) - center[1],
            ),
        )
        offset = parse_ci(row["offset"], f"{label}.offset")
        require(
            ci_contains(offset, expected_offset),
            f"{label} offset box misses exact geometry",
        )
        center_value = parse_ci(row["center_value"], f"{label}.center_value")
        derivative = parse_ci(
            row["derivative_hull"], f"{label}.derivative"
        )
        image = parse_ci(row["image"], f"{label}.image")
        require(
            ci_contains(image, c_add(center_value, c_mul(derivative, offset))),
            f"{label} centered-form image does not replay",
        )
        require(
            ci_contains(image, center_value),
            f"{label} midpoint value escapes its segment image",
        )
        require(not ci_contains_zero(image), f"{label} image contains zero")
        require(
            parse_fraction(
                row["image_zero_exclusion_abs2_lower"],
                f"{label}.image_zero_slack",
            )
            == ci_zero_exclusion_abs2_lower(image)
            > 0,
            f"{label} image zero-exclusion witness drift",
        )
        rotated = parse_ci(row["rotated_image"], f"{label}.rotated")
        require(
            ci_contains(rotated, c_mul(image, c_conj(center_value))),
            f"{label} rotated enclosure does not replay",
        )
        require(
            parse_fraction(
                row["rotated_zero_exclusion_abs2_lower"],
                f"{label}.rotated_zero_slack",
            )
            == ci_zero_exclusion_abs2_lower(rotated)
            > 0,
            f"{label} rotated zero-exclusion witness drift",
        )
        require(excludes_principal_cut(rotated), f"{label} rotated chart meets cut")
        rotated_arg = parse_ri(
            row["rotated_argument"], f"{label}.rotated_arg"
        )
        require(
            ri_contains(rotated_arg, directed_arg(rotated, precision)),
            f"{label} rotated argument is not independently enclosed",
        )
        argument_width = rotated_arg[1] - rotated_arg[0]
        require(argument_width <= gate, f"{label} cone gate fails")
        require(
            parse_fraction(row["rotated_argument_width"], f"{label}.arg_width")
            == argument_width
            and parse_fraction(
                row["argument_width_slack"], f"{label}.arg_slack"
            )
            == gate - argument_width,
            f"{label} argument-width evidence drift",
        )
        start_value = parse_ci(row["start_value"], f"{label}.start_value")
        end_value = parse_ci(row["end_value"], f"{label}.end_value")
        require(
            ci_contains(image, start_value) and ci_contains(image, end_value),
            f"{label} endpoint value escapes its segment image",
        )
        if previous_end_value is None:
            first_start_value = start_value
        else:
            require(
                start_value == previous_end_value,
                f"{label} endpoint values do not form one closed walk",
            )
        previous_end_value = end_value
        ratio = parse_ci(
            row["endpoint_ratio"], f"{label}.endpoint_ratio"
        )
        require(
            ci_contains(ratio, c_mul(end_value, c_conj(start_value))),
            f"{label} endpoint ratio does not replay",
        )
        require(excludes_principal_cut(ratio), f"{label} endpoint ratio meets cut")
        increment = parse_ri(
            row["endpoint_increment"], f"{label}.increment"
        )
        require(
            ri_contains(increment, directed_arg(ratio, precision)),
            f"{label} endpoint increment is not independently enclosed",
        )
        require(
            abs(increment[0]) <= gate and abs(increment[1]) <= gate,
            f"{label} endpoint increment gate fails",
        )
        require(
            parse_fraction(
                row["endpoint_increment_slack"],
                f"{label}.endpoint_increment_slack",
            )
            == gate - max(abs(increment[0]), abs(increment[1])),
            f"{label} endpoint increment slack drift",
        )
        increment_sum = r_add(increment_sum, increment)
        expected_segment_id = sha256_bytes(
            canonical_json(
                {"start": row["start"], "end": row["end"]}
            ).encode("utf-8")
        )
        require(
            row["segment_id"] == expected_segment_id,
            f"{label} segment id mismatch",
        )
        coordinates.append(
            {
                "segment_id": row["segment_id"],
                "start": row["start"],
                "end": row["end"],
            }
        )
        parsed_records.append(
            {
                "segment_id": row["segment_id"],
                "segment_box": (
                    (min(start[0], end[0]), max(start[0], end[0])),
                    (min(start[1], end[1]), max(start[1], end[1])),
                ),
                "center_value": center_value,
                "derivative_hull": derivative,
                "offset": offset,
                "image": image,
                "rotated_image": rotated,
                "rotated_argument": rotated_arg,
                "start_value": start_value,
                "end_value": end_value,
                "endpoint_ratio": ratio,
                "endpoint_increment": increment,
            }
        )
    require(previous_end == first_start, f"{label} partition is not closed")
    require(
        first_start == (re_lo, im_lo),
        f"{label} partition does not start at the lower-left corner",
    )
    edge_corners = (
        ((re_lo, im_lo), (re_hi, im_lo)),
        ((re_hi, im_lo), (re_hi, im_hi)),
        ((re_hi, im_hi), (re_lo, im_hi)),
        ((re_lo, im_hi), (re_lo, im_lo)),
    )
    require(
        all(
            group
            and group[0][0] == expected_start
            and group[-1][1] == expected_end
            for group, (expected_start, expected_end) in zip(
                edge_groups, edge_corners
            )
        ),
        f"{label} partition does not cover every edge exactly once",
    )
    require(
        winding["max_depth_used"]
        == (actual_max_depth if precision == 128 else 0),
        f"{label} max_depth_used does not replay from dyadic geometry",
    )
    require(
        previous_end_value == first_start_value,
        f"{label} endpoint-value walk does not close",
    )
    require(signed_area > 0, f"{label} boundary orientation is not positive")
    require(
        winding["partition_sha256"]
        == sha256_bytes(canonical_json(coordinates).encode("utf-8")),
        f"{label} partition digest mismatch",
    )

    total = parse_ri(winding["total_variation_interval"], f"{label}.total")
    require(ri_contains(total, increment_sum), f"{label} total does not contain sum")
    iv.prec = precision
    two_pi_raw = iv.pi * iv.mpf(2)
    independent_two_pi = (
        _mpf_tuple_fraction(two_pi_raw._mpi_[0]),
        _mpf_tuple_fraction(two_pi_raw._mpi_[1]),
    )
    candidate = winding["winding"]
    require(candidate == 1, f"{label} winding is not one")
    expected_residual = r_sub(
        total,
        (
            independent_two_pi[0] * candidate,
            independent_two_pi[1] * candidate,
        ),
    )
    residual = parse_ri(winding["winding_residual"], f"{label}.residual")
    require(
        ri_contains(residual, expected_residual),
        f"{label} winding residual does not replay",
    )
    within = residual[0] > -tolerance and residual[1] < tolerance
    require(
        parse_fraction(
            winding["winding_tolerance_slack"],
            f"{label}.winding_tolerance_slack",
        )
        == tolerance - max(abs(residual[0]), abs(residual[1])),
        f"{label} winding tolerance slack drift",
    )
    require(winding["certified"] is within, f"{label} winding flag does not replay")
    return within, parsed_records, total, [row["segment_id"] for row in records]


def strict_inclusion(source: CI, image: CI) -> tuple[dict[str, Fraction], Fraction]:
    margins = {
        "re_lower": image[0][0] - source[0][0],
        "re_upper": source[0][1] - image[0][1],
        "im_lower": image[1][0] - source[1][0],
        "im_upper": source[1][1] - image[1][1],
    }
    return margins, min(margins.values())


def check_declared_chart_quantities(
    rows: list[dict[str, Any]],
    segment_records: list[dict[str, Any]],
    crossed: set[tuple[Fraction, Fraction]],
    precision: int,
    label: str,
) -> dict[str, CI]:
    """Replay the explicit base/correction/combined chart decomposition."""

    segment_boxes = {
        record["segment_id"]: record["segment_box"]
        for record in segment_records
    }
    segment_ids = list(segment_boxes)
    expected_keys = {
        (segment_id, pair)
        for segment_id in segment_ids
        for pair in crossed
    }
    seen = set()
    flattened: dict[str, CI] = {}
    for index, row in enumerate(rows):
        segment_id = row["segment_id"]
        pair = (
            parse_fraction(row["m1"], f"{label}.chart[{index}].m1"),
            parse_fraction(row["m2"], f"{label}.chart[{index}].m2"),
        )
        require(
            (segment_id, pair) in expected_keys
            and (segment_id, pair) not in seen,
            f"{label} chart-quantity census drift",
        )
        seen.add((segment_id, pair))
        base_chart, formula = expected_base_and_correction(*pair)
        require(
            row["base_chart"] == base_chart
            and row["crossing_correction_formula"] == formula,
            f"{label} chart-quantity label drift",
        )
        base = parse_ci(row["base_B0"], f"{label}.base_B0")
        correction = parse_ci(
            row["chart_correction_B0"], f"{label}.correction_B0"
        )
        combined = parse_ci(row["combined_B0"], f"{label}.combined_B0")
        base_slope = parse_ci(row["base_B0p"], f"{label}.base_B0p")
        correction_slope = parse_ci(
            row["chart_correction_B0p"], f"{label}.correction_B0p"
        )
        combined_slope = parse_ci(
            row["combined_B0p"], f"{label}.combined_B0p"
        )
        require(
            ci_contains(combined, c_add(base, correction)),
            f"{label} combined B0 does not contain base plus correction",
        )
        require(
            ci_contains(combined_slope, c_add(base_slope, correction_slope)),
            f"{label} combined B0p does not contain derivative sum",
        )
        expected_correction, expected_correction_slope = (
            expected_chart_correction(
                pair,
                segment_boxes[segment_id],
                precision,
            )
        )
        require(
            ci_contains(correction, expected_correction),
            f"{label} chart correction does not replay from its formula",
        )
        require(
            ci_contains(correction_slope, expected_correction_slope),
            f"{label} chart-correction derivative does not replay",
        )
        rendered = f"{pair[0]},{pair[1]}"
        flattened.update(
            {
                f"segment:{segment_id}:base_chart:B0({rendered})": base,
                f"segment:{segment_id}:chart_correction:B0({rendered})": correction,
                f"segment:{segment_id}:combined:B0({rendered})": combined,
                f"segment:{segment_id}:base_chart:B0p({rendered})": base_slope,
                f"segment:{segment_id}:chart_correction:B0p({rendered})": correction_slope,
                f"segment:{segment_id}:combined:B0p({rendered})": combined_slope,
            }
        )
    require(seen == expected_keys and seen, f"{label} chart evidence is incomplete")
    return flattened


def expected_quantity_labels(
    name: str,
    segment_ids: list[str],
    crossed: set[tuple[Fraction, Fraction]],
) -> set[str]:
    """Exact semantic census for every serialized nesting quantity."""

    labels = {
        "probe:center:inverse_propagator",
        "probe:center:inverse_propagator_derivative",
    }
    labels.update(
        f"probe:center:integral:A0({mass})"
        for mass in EXPECTED_A0_MASSES
    )
    for pair in EXPECTED_PAIRS[name]:
        rendered = f"{pair[0]},{pair[1]}"
        labels.add(f"probe:center:base_chart:B0({rendered})")
        labels.add(f"probe:center:base_chart:B0p({rendered})")
        if pair in crossed:
            labels.update(
                {
                    f"probe:center:chart_correction:B0({rendered})",
                    f"probe:center:combined:B0({rendered})",
                    f"probe:center:chart_correction:B0p({rendered})",
                    f"probe:center:combined:B0p({rendered})",
                }
            )
        for segment_id in segment_ids:
            labels.add(
                f"segment:{segment_id}:base_chart:B0({rendered})"
            )
            labels.add(
                f"segment:{segment_id}:base_chart:B0p({rendered})"
            )
            if pair in crossed:
                labels.update(
                    {
                        f"segment:{segment_id}:chart_correction:B0({rendered})",
                        f"segment:{segment_id}:combined:B0({rendered})",
                        f"segment:{segment_id}:chart_correction:B0p({rendered})",
                        f"segment:{segment_id}:combined:B0p({rendered})",
                    }
                )
    return labels


def check_newton(
    newton: dict[str, Any], name: str, precision: int, label: str
) -> tuple[bool, dict[str, CI]]:
    require(precision in PRECISIONS, f"{label} Newton precision drift")
    seed = parse_rational_box(newton["seed_box"], f"{label}.seed")
    require(seed == NEWTON_SEEDS[name], f"{label} Newton seed drift")
    iterations = newton["iterations"]
    require(
        1 <= len(iterations) <= 4,
        f"{label} Newton iteration census drift",
    )
    source = seed
    first_strict = None
    first_evidence = None
    for expected_index, row in enumerate(iterations):
        require(row["step_index"] == expected_index, f"{label} Newton index drift")
        recorded_source = parse_ci(row["source_box"], f"{label}.source")
        if expected_index == 0:
            require(
                ci_contains(recorded_source, source),
                f"{label} initial Newton source misses exact seed",
            )
        else:
            require(
                recorded_source == source,
                f"{label} Newton source-chain drift",
            )
        source = recorded_source
        midpoint = (
            parse_fraction(row["midpoint"][0], f"{label}.midpoint.re"),
            parse_fraction(row["midpoint"][1], f"{label}.midpoint.im"),
        )
        require(
            midpoint
            == (
                (source[0][0] + source[0][1]) / 2,
                (source[1][0] + source[1][1]) / 2,
            ),
            f"{label} Newton midpoint is not exact",
        )
        midpoint_box: CI = ((midpoint[0], midpoint[0]), (midpoint[1], midpoint[1]))
        g_mid = parse_ci(row["midpoint_value"], f"{label}.g_mid")
        derivative = parse_ci(
            row["derivative_on_source"], f"{label}.derivative"
        )
        require(not ci_contains_zero(derivative), f"{label} Newton derivative contains zero")
        image = parse_ci(row["newton_image"], f"{label}.newton_image")
        require(
            ci_contains(image, c_sub(midpoint_box, c_div(g_mid, derivative))),
            f"{label} Newton image arithmetic does not replay",
        )
        margins, minimum = strict_inclusion(source, image)
        inclusion = row["strict_inclusion"]
        require(
            {
                key: parse_fraction(value, f"{label}.margin.{key}")
                for key, value in inclusion["margins"].items()
            }
            == margins,
            f"{label} strict-inclusion margins drift",
        )
        require(
            parse_fraction(inclusion["minimum_margin"], f"{label}.minimum_margin")
            == minimum,
            f"{label} minimum margin drift",
        )
        is_strict = minimum > 0
        require(inclusion["strict"] is is_strict, f"{label} strict flag drift")
        if is_strict and first_strict is None:
            first_strict = expected_index
            first_evidence = inclusion
        intersection = c_intersection(source, image)
        require(intersection is not None, f"{label} Newton intersection empty")
        recorded_intersection = parse_ci(
            row["intersection_box"], f"{label}.intersection"
        )
        require(
            recorded_intersection == intersection,
            f"{label} Newton intersection drift",
        )
        source = recorded_intersection
    require(first_strict is not None, f"{label} has no strict Newton inclusion")
    require(
        newton["proof_step_index"] == first_strict
        and newton["strict_inclusion_margin"] == first_evidence,
        f"{label} strict Newton proof witness drift",
    )
    require(
        newton["strict_interior_inclusion_certified"] is True
        and newton["contracted"] is True,
        f"{label} Newton certification does not require strict inclusion",
    )
    pole = parse_ci(newton["pole_ball"], f"{label}.pole_ball")
    require(pole == source, f"{label} final pole ball is not final intersection")
    derivative = parse_ci(newton["derivative_ball"], f"{label}.derivative_ball")
    residual = parse_ci(newton["residual_ball"], f"{label}.residual_ball")
    residue = parse_ci(newton["residue_ball"], f"{label}.residue_ball")
    denominator_excludes = not ci_contains_zero(derivative)
    residual_contains = ci_contains_zero(residual)
    require(
        newton["denominator_ball_excludes_zero"] is denominator_excludes,
        f"{label} denominator exclusion flag drift",
    )
    expected_residue = c_div(
        ((Fraction(1), Fraction(1)), (Fraction(0), Fraction(0))),
        derivative,
    )
    require(ci_contains(residue, expected_residue), f"{label} residue inversion drift")
    null = newton["null_vectors"]
    require(
        null["projection"] == "transverse_scalar_block_1x1"
        and null["full_matrix_rank_n_minus_1_certified"] is False
        and null["left_residual_contains_zero"] is residual_contains
        and null["right_residual_contains_zero"] is residual_contains
        and null["laurent_denominator_excludes_zero"] is denominator_excludes,
        f"{label} scalar Laurent evidence drift",
    )
    return (
        bool(denominator_excludes and residual_contains),
        {
            "pole_ball": pole,
            "derivative_ball": derivative,
            "residual_ball": residual,
            "residue_ball": residue,
        },
    )


def validate_receipt(
    receipt: dict[str, Any], schema: dict[str, Any] | None = None
) -> None:
    """Validate an in-memory receipt, raising ``ReceiptValidationError``."""

    if schema is None:
        schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    check_schema(receipt, schema)
    check_bindings(receipt)
    check_scope_and_gates(receipt)

    row_evidence: dict[str, dict[int, dict[str, Any]]] = {"W": {}, "Z": {}}
    all_certified = True
    for name in ("W", "Z"):
        crossed = check_continuation(receipt, name)
        rows = receipt["results"][name]
        for precision in PRECISIONS:
            row = rows[str(precision)]
            label = f"{name}@{precision}"
            require(row["particle"] == name, f"{label} particle drift")
            require(row["precision_bits"] == precision, f"{label} precision drift")
            require(
                row["pole_box"]["half_plane"] == "lower",
                f"{label} pole-box half-plane drift",
            )
            pole_box = parse_rational_box(
                {
                    "re": row["pole_box"]["re"],
                    "im": row["pole_box"]["im"],
                },
                f"{label}.pole_box",
            )
            expected_pole: CI = POLE_BOXES[name]
            require(pole_box == expected_pole, f"{label} pole box drift")
            require(
                parse_fraction(row["tree_mass_sq"], f"{label}.tree")
                == TREE_MASSES[name],
                f"{label} tree mass drift",
            )
            interior_ok, interior_denominators = check_holomorphy(
                row["interior_holomorphy"],
                POLE_BOXES[name],
                EXPECTED_PAIRS[name],
                crossed,
                EXPECTED_DENOMINATOR_COUNTS[name],
                name,
                f"{label}.interior",
            )
            corridor_ok, corridor_denominators = check_holomorphy(
                row["corridor_holomorphy"],
                CORRIDOR_BOXES[name],
                EXPECTED_PAIRS[name],
                crossed,
                EXPECTED_DENOMINATOR_COUNTS[name],
                name,
                f"{label}.corridor",
            )
            winding_ok, segment_records, total, segment_ids = check_winding(
                row["boundary_winding"], POLE_BOXES[name], precision, label
            )
            chart_quantities = check_declared_chart_quantities(
                row["declared_chart_segment_quantities"],
                segment_records,
                crossed,
                precision,
                label,
            )
            newton_ok, balls = check_newton(
                row["interval_newton"], name, precision, label
            )
            probes = {
                probe_name: parse_ci(value, f"{label}.probe.{probe_name}")
                for probe_name, value in row["quantity_enclosures"].items()
            }
            expected_labels = expected_quantity_labels(
                name, segment_ids, crossed
            )
            require(
                set(probes) == expected_labels,
                f"{label} quantity-enclosure semantic census drift",
            )
            require(
                all(
                    chart_label in probes
                    and probes[chart_label] == chart_interval
                    for chart_label, chart_interval in chart_quantities.items()
                ),
                f"{label} chart quantities are not bound into nesting evidence",
            )
            require(
                all(
                    (
                        "base_chart" in quantity_label
                        or "chart_correction" in quantity_label
                        or "combined" in quantity_label
                    )
                    for quantity_label in probes
                    if ":B0(" in quantity_label or ":B0p(" in quantity_label
                ),
                f"{label} ambiguous base B0/B0p quantity label",
            )
            simple = (
                interior_ok
                and winding_ok
                and row["boundary_winding"]["winding"] == 1
            )
            expected = interior_ok and corridor_ok and simple and newton_ok
            require(
                row["simple_scalar_root_certified"] is simple,
                f"{label} simple-scalar-root flag drift",
            )
            require(
                row["scalar_pole_certified"] is expected,
                f"{label} scalar-pole flag drift",
            )
            require(
                row["reading"]
                == (CERTIFIED_READING if expected else "not certified"),
                f"{label} claim-bearing reading drift",
            )
            all_certified = all_certified and expected
            row_evidence[name][precision] = {
                "segments": segment_records,
                "total": total,
                "quantities": probes,
                "balls": balls,
                "denominators": {
                    **{
                        f"pole_box:{identifier}": enclosure
                        for identifier, enclosure in interior_denominators.items()
                    },
                    **{
                        f"corridor_box:{identifier}": enclosure
                        for identifier, enclosure in corridor_denominators.items()
                    },
                },
                "partition_sha256": row["boundary_winding"][
                    "partition_sha256"
                ],
            }

        comparisons = []
        aggregate_quantity: dict[str, bool] = {}
        aggregate_balls = {
            label: True
            for label in ("pole_ball", "derivative_ball", "residual_ball", "residue_ball")
        }
        segment_field_map = {
            "center_values_nested": ("center_value",),
            "derivative_hulls_nested": ("derivative_hull",),
            "offsets_nested": ("offset",),
            "images_nested": ("image",),
            "rotated_images_nested": ("rotated_image",),
            "rotated_arguments_nested": ("rotated_argument",),
            "endpoint_values_nested": ("start_value", "end_value"),
            "endpoint_ratios_nested": ("endpoint_ratio",),
            "endpoint_increments_nested": ("endpoint_increment",),
        }
        real_segment_fields = {
            "rotated_argument",
            "endpoint_increment",
        }
        aggregate_segment_fields = {
            summary_label: True for summary_label in segment_field_map
        }
        aggregate_total = True
        aggregate_partitions = True
        aggregate_denominators = True
        quantity_key_sets_match = True
        denominator_key_sets_match = True
        for outer_precision, inner_precision in zip(PRECISIONS, PRECISIONS[1:]):
            outer = row_evidence[name][outer_precision]
            inner = row_evidence[name][inner_precision]
            quantity_pair_keys_match = (
                set(outer["quantities"]) == set(inner["quantities"])
            )
            denominator_pair_keys_match = (
                set(outer["denominators"]) == set(inner["denominators"])
            )
            total_nested = ri_contains(outer["total"], inner["total"])
            partition_identical = (
                outer["partition_sha256"] == inner["partition_sha256"]
                and [record["segment_id"] for record in outer["segments"]]
                == [record["segment_id"] for record in inner["segments"]]
            )
            segment_field_pair: dict[str, bool] = {}
            for summary_label, fields in segment_field_map.items():
                segment_field_pair[summary_label] = bool(
                    len(outer["segments"]) == len(inner["segments"])
                    and all(
                        (
                            ri_contains(
                                outer_record[field],
                                inner_record[field],
                            )
                            if field in real_segment_fields
                            else ci_nested(
                                outer_record[field],
                                inner_record[field],
                            )
                        )
                        for outer_record, inner_record in zip(
                            outer["segments"], inner["segments"]
                        )
                        for field in fields
                    )
                )
            segment_nested = bool(
                partition_identical and all(segment_field_pair.values())
            )
            quantity_nested = {
                label: ci_nested(
                    outer["quantities"][label],
                    inner["quantities"][label],
                )
                for label in sorted(outer["quantities"])
            } if quantity_pair_keys_match else {}
            balls_nested = {
                label: ci_nested(outer["balls"][label], inner["balls"][label])
                for label in aggregate_balls
            }
            denominators_nested = bool(
                denominator_pair_keys_match
                and all(
                    ci_nested(
                        outer["denominators"][identifier],
                        inner["denominators"][identifier],
                    )
                    for identifier in outer["denominators"]
                )
            )
            comparisons.append(
                {
                    "outer_precision": outer_precision,
                    "inner_precision": inner_precision,
                    "total_variation_nested": total_nested,
                    "partition_identical": partition_identical,
                    "segment_enclosures_nested": segment_nested,
                    "quantity_enclosure_nesting": quantity_nested,
                    "newton_ball_nesting": balls_nested,
                    "coefficient_denominators_nested": denominators_nested,
                }
            )
            aggregate_total = aggregate_total and total_nested
            aggregate_partitions = aggregate_partitions and partition_identical
            quantity_key_sets_match = (
                quantity_key_sets_match and quantity_pair_keys_match
            )
            denominator_key_sets_match = (
                denominator_key_sets_match and denominator_pair_keys_match
            )
            aggregate_denominators = (
                aggregate_denominators and denominators_nested
            )
            for summary_label, nested in segment_field_pair.items():
                aggregate_segment_fields[summary_label] = (
                    aggregate_segment_fields[summary_label] and nested
                )
            for label, nested in quantity_nested.items():
                aggregate_quantity[label] = aggregate_quantity.get(label, True) and nested
            for label, nested in balls_nested.items():
                aggregate_balls[label] = aggregate_balls[label] and nested

        nesting = receipt["precision_nesting"][name]
        require(nesting["comparison_pairs"] == comparisons, f"{name} comparison replay drift")
        require(
            nesting["enclosures_nested_with_precision"] is aggregate_total,
            f"{name} total nesting summary drift",
        )
        require(
            nesting["per_quantity_enclosure_nesting"] == aggregate_quantity,
            f"{name} quantity nesting summary drift",
        )
        require(
            nesting["quantity_key_sets_match"] is quantity_key_sets_match
            and nesting["quantity_enclosures_all_nested"]
            is (
                quantity_key_sets_match
                and bool(aggregate_quantity)
                and all(aggregate_quantity.values())
            ),
            f"{name} quantity key-set/aggregate nesting drift",
        )
        expected_segment_summary = {
            "segments": len(row_evidence[name][128]["segments"]),
            "ladders_compared": 2,
            "partition_ids_match": aggregate_partitions,
            **aggregate_segment_fields,
            "all_nested": bool(
                aggregate_partitions
                and all(aggregate_segment_fields.values())
            ),
        }
        require(
            nesting["per_segment_enclosure_nesting"]
            == expected_segment_summary,
            f"{name} segment nesting summary drift",
        )
        require(
            nesting["newton_ball_nesting"] == aggregate_balls,
            f"{name} Newton nesting summary drift",
        )
        require(
            nesting["coefficient_denominator_nesting"]
            == {
                "records": 2 * EXPECTED_DENOMINATOR_COUNTS[name],
                "key_sets_match": denominator_key_sets_match,
                "all_nested": aggregate_denominators,
            },
            f"{name} coefficient-denominator nesting summary drift",
        )
        nesting_ok = (
            aggregate_total
            and aggregate_partitions
            and all(aggregate_segment_fields.values())
            and quantity_key_sets_match
            and all(aggregate_quantity.values())
            and all(aggregate_balls.values())
            and denominator_key_sets_match
            and aggregate_denominators
        )
        all_certified = all_certified and nesting_ok

    expected_status = SUCCESS_STATUS if all_certified else "CERTIFICATION_INCOMPLETE"
    require(receipt["status"] == expected_status, "status does not follow from evidence")
    expected_promotion = {
        "scalar_pole_on_declared_algebraic_chart_certified": all_certified,
        "simple_scalar_root_certified": all_certified,
        "scalar_laurent_denominator_ball_certified": all_certified,
        "scalar_residue_ball_certified": all_certified,
        "continuation_identity_independently_certified": False,
        "standard_second_sheet_identification_certified": False,
        "theorem_self_energy_sign_bridge_certified": False,
        "full_matrix_rank_n_minus_1_laurent_certified": False,
        "issue_593_precision_ladder_row_discharged": False,
        "issue_593_independent_third_verifier_row_discharged": False,
        "issue_593_pole_laurent_row_discharged": False,
        "issue_593_full_acceptance_satisfied": False,
        "bmhv_restoration_certified": False,
        "physical_current_claim": False,
        "oph_native": False,
        "unit_claim": False,
        "unitarity_claim": False,
    }
    require(receipt["promotion"] == expected_promotion, "promotion scope drift")
    require(
        all_certified,
        "positive scalar declared-chart pole certificate is not established",
    )


def main() -> int:
    try:
        receipt = json.loads(ARTIFACT_PATH.read_text(encoding="utf-8"))
        validate_receipt(receipt)
    except (OSError, json.JSONDecodeError, ReceiptValidationError) as error:
        print(f"SCALAR_DECLARED_CHART_POLE_CHECK_FAIL: {error}", file=sys.stderr)
        return 1
    print(PASS_STATUS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
