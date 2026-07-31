#!/usr/bin/env python3
"""Fail-closed evidence checker for the principal-sheet W/Z receipt.

This checker is deliberately independent of the numerical producer: it does
not import ``certified_wz_contours``, ``complex_interval``, or
``wz_pole_receipts``.  It validates the strict v3 schema, the receipt and
source digests, the frozen external-SM fixture, the exact finite
dimensional-prefactor corrections, the complete chart/quantity/partition
keysets, and the arithmetic relations exposed by the serialized witnesses.
In particular, it recomputes coefficient-denominator exclusion, boundary
coverage, every recorded segment inequality, the winding-zero residual, and
128/192/256-bit enclosure nesting from exact rational endpoints.

It does *not* independently evaluate the loop functions or prove that the
recorded interval boxes enclose the physical expressions.  A pass is therefore
a structural/evidence validation of the committed auxiliary principal-sheet
zero-exclusion receipt, not the independent numerical third verifier required
by issue #593.  It is not a resonance-root, Laurent-denominator,
theorem-sign-bridge, physical-current, physical-unit, or OPH-native
certificate.
"""

from __future__ import annotations

import copy
import hashlib
import json
import re
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

from jsonschema import Draft202012Validator
from mpmath import iv
from mpmath.libmp import to_rational

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_PATH = ROOT / "outputs" / "certified_wz_contours.json"
VECTOR_PATH = ROOT / "outputs" / "fj_direct_vector_blocks.json"
SCHEMA_PATH = ROOT / "schemas" / "certified_wz_contours_v3.schema.json"

PASS_STATUS = (
    "CERTIFIED_CONTOUR_EVIDENCE_CHECK_PASS"
    "__NO_INDEPENDENT_PHYSICS_REEVALUATION"
)
FAIL_STATUS = "CERTIFIED_CONTOUR_EVIDENCE_CHECK_FAIL"

PRECISIONS = ("128", "192", "256")
INITIAL_SEGMENTS_PER_EDGE = 8
MAX_SUBDIVISION_DEPTH = 12
ARG_WIDTH_GATE = Fraction(157, 100)
WINDING_TOLERANCE = Fraction(157, 100)
EXPECTED_FIXTURE = {
    "g1": "1/4",
    "g2": "1/3",
    "v": "2",
    "lam": "1/8",
    "xi": "1",
    "mu2": "1/2",
    "mfu1": "1/50",
    "mfu2": "1/20",
    "mfu3": "1/5",
    "mfd1": "1/60",
    "mfd2": "1/25",
    "mfd3": "1/10",
    "mfe1": "1/80",
    "mfe2": "1/30",
    "mfe3": "1/15",
    "mu_ren2": "1",
}
EXPECTED_CKM = {
    f"V{i}{j}": str(int(i == j))
    for i in (1, 2, 3)
    for j in (1, 2, 3)
}
EXPECTED_VECTOR_UNITS = (
    "loop measure i/(16 pi^2) stripped; "
    "Delta is the single 1/eps pole unit"
)
EXPECTED_CORRECTION_TEXT = {
    "W": "(312000*p2 + 90163)/1944000",
    "Z": "(292020000*p2 + 79268959)/1944000000",
}
EXPECTED_CORRECTION_COEFFICIENTS = {
    "W": {1: Fraction(13, 81), 0: Fraction(90163, 1944000)},
    "Z": {
        1: Fraction(4867, 32400),
        0: Fraction(79268959, 1944000000),
    },
}
BOXES = {
    "W": {
        "re": (Fraction(109, 1000), Fraction(114, 1000)),
        "im": (Fraction(1, 5000), Fraction(1, 500)),
    },
    "Z": {
        "re": (Fraction(172, 1000), Fraction(176, 1000)),
        "im": (Fraction(1, 5000), Fraction(3, 2000)),
    },
}
TREE_MASSES = {"W": Fraction(1, 9), "Z": Fraction(25, 144)}
EXPECTED_DENOMINATOR_COUNTS = {"W": 67, "Z": 40}
EXPECTED_DENOMINATOR_IDENTITY_DIGESTS = {
    "W": "sha256:fa9fb5908600498006f162972af47c161549968d7042360dce07972f90a7f2f9",
    "Z": "sha256:532c45f89ec092a287778796fdcd1b1d87e5b18a1f83149734f4c2f3589fd6cc",
}
EXPECTED_SEGMENT_COUNTS = {"W": 49, "Z": 41}
EXPECTED_QUANTITY_COUNTS = {"W": 1115, "Z": 1023}
EXPECTED_A0_MASSES = {
    "0",
    "1",
    "1/100",
    "1/225",
    "1/25",
    "1/2500",
    "1/3600",
    "1/400",
    "1/625",
    "1/6400",
    "1/9",
    "1/900",
    "25/144",
}
EXPECTED_CHART_PAIRS = {
    "W": {
        ("0", "1/9"),
        ("0", "1/6400"),
        ("0", "1/900"),
        ("0", "1/225"),
        ("25/144", "1/9"),
        ("1/9", "1"),
        ("1/9", "25/144"),
        ("1/9", "0"),
        ("1/2500", "1/3600"),
        ("1/400", "1/625"),
        ("1/25", "1/100"),
    },
    "Z": {
        ("0", "0"),
        ("25/144", "1"),
        ("1/9", "1/9"),
        ("1/3600", "1/3600"),
        ("1/625", "1/625"),
        ("1/100", "1/100"),
        ("1/6400", "1/6400"),
        ("1/900", "1/900"),
        ("1/225", "1/225"),
        ("1/2500", "1/2500"),
        ("1/400", "1/400"),
        ("1/25", "1/25"),
    },
}
CORRECTION_RE = re.compile(
    r"^\((-?[0-9]+)\*p2 ([+-]) ([0-9]+)\)/([1-9][0-9]*)$"
)


RealInterval = tuple[Fraction, Fraction]
ComplexInterval = tuple[RealInterval, RealInterval]
Point = tuple[Fraction, Fraction]


class ReceiptValidationError(ValueError):
    """Raised by :func:`validate_receipt` with all fail-closed problems."""

    def __init__(self, problems: Sequence[str]) -> None:
        self.problems = list(problems)
        super().__init__("; ".join(self.problems))


def sha256_bytes(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"


def canonical_digest(value: Any) -> str:
    return sha256_bytes(canonical_json(value).encode("utf-8"))


def _fraction(value: Any, context: str, problems: list[str]) -> Fraction:
    try:
        return Fraction(str(value))
    except (ValueError, ZeroDivisionError):
        problems.append(f"{context}: invalid exact rational {value!r}")
        return Fraction(0)


def _real_interval(
    value: Mapping[str, Any], context: str, problems: list[str]
) -> RealInterval:
    lo = _fraction(value.get("lo"), f"{context}.lo", problems)
    hi = _fraction(value.get("hi"), f"{context}.hi", problems)
    if lo > hi:
        problems.append(f"{context}: lower endpoint exceeds upper endpoint")
    return lo, hi


def _complex_interval(
    value: Mapping[str, Any], context: str, problems: list[str]
) -> ComplexInterval:
    return (
        _real_interval(value.get("re", {}), f"{context}.re", problems),
        _real_interval(value.get("im", {}), f"{context}.im", problems),
    )


def _point(value: Sequence[Any], context: str, problems: list[str]) -> Point:
    if len(value) != 2:
        problems.append(f"{context}: point does not have two coordinates")
        return Fraction(0), Fraction(0)
    return (
        _fraction(value[0], f"{context}[0]", problems),
        _fraction(value[1], f"{context}[1]", problems),
    )


def _ri_nested(outer: RealInterval, inner: RealInterval) -> bool:
    return outer[0] <= inner[0] and inner[1] <= outer[1]


def _ci_nested(outer: ComplexInterval, inner: ComplexInterval) -> bool:
    return _ri_nested(outer[0], inner[0]) and _ri_nested(outer[1], inner[1])


def _ci_contains_zero(value: ComplexInterval) -> bool:
    return (
        value[0][0] <= 0 <= value[0][1]
        and value[1][0] <= 0 <= value[1][1]
    )


def _square_lower(interval: RealInterval) -> Fraction:
    lo, hi = interval
    if lo <= 0 <= hi:
        return Fraction(0)
    return min(lo * lo, hi * hi)


def _abs2_lower(value: ComplexInterval) -> Fraction:
    return _square_lower(value[0]) + _square_lower(value[1])


def _atan2_replay(
    value: ComplexInterval,
    precision_bits: int,
    context: str,
    problems: list[str],
) -> RealInterval:
    """Independently replay ``arg`` at the row's declared ball precision."""

    old_precision = iv.prec
    try:
        iv.prec = precision_bits
        real = iv.mpf([str(value[0][0]), str(value[0][1])])
        imag = iv.mpf([str(value[1][0]), str(value[1][1])])
        angle = iv.atan2(imag, real)
        lower_raw, upper_raw = angle._mpi_
        lower_num, lower_den = to_rational(lower_raw)
        upper_num, upper_den = to_rational(upper_raw)
        return Fraction(lower_num, lower_den), Fraction(upper_num, upper_den)
    except (ValueError, ZeroDivisionError) as exc:
        problems.append(f"{context}: independent atan2 replay failed: {exc}")
        return Fraction(0), Fraction(0)
    finally:
        iv.prec = old_precision


def _ri_add(left: RealInterval, right: RealInterval) -> RealInterval:
    return left[0] + right[0], left[1] + right[1]


def _ri_mul(left: RealInterval, right: RealInterval) -> RealInterval:
    products = (
        left[0] * right[0],
        left[0] * right[1],
        left[1] * right[0],
        left[1] * right[1],
    )
    return min(products), max(products)


def _ci_add(left: ComplexInterval, right: ComplexInterval) -> ComplexInterval:
    return _ri_add(left[0], right[0]), _ri_add(left[1], right[1])


def _ci_mul(left: ComplexInterval, right: ComplexInterval) -> ComplexInterval:
    ac = _ri_mul(left[0], right[0])
    bd = _ri_mul(left[1], right[1])
    ad = _ri_mul(left[0], right[1])
    bc = _ri_mul(left[1], right[0])
    return (
        (ac[0] - bd[1], ac[1] - bd[0]),
        (ad[0] + bc[0], ad[1] + bc[1]),
    )


def _ci_conjugate(value: ComplexInterval) -> ComplexInterval:
    return value[0], (-value[1][1], -value[1][0])


def _ci_pow(value: ComplexInterval, power: int) -> ComplexInterval:
    result: ComplexInterval = ((Fraction(1), Fraction(1)), (Fraction(0), Fraction(0)))
    for _ in range(power):
        result = _ci_mul(result, value)
    return result


def _ci_scale(value: ComplexInterval, coefficient: Fraction) -> ComplexInterval:
    scale = (coefficient, coefficient)
    return _ri_mul(value[0], scale), _ri_mul(value[1], scale)


def _evaluate_polynomial(
    coefficients: Sequence[Sequence[Any]],
    box: Mapping[str, tuple[Fraction, Fraction]],
    context: str,
    problems: list[str],
) -> ComplexInterval:
    s_box: ComplexInterval = (box["re"], box["im"])
    total: ComplexInterval = (
        (Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(0)),
    )
    powers: set[int] = set()
    previous_power = -1
    for index, term in enumerate(coefficients):
        if len(term) != 2:
            problems.append(f"{context}[{index}]: malformed polynomial term")
            continue
        try:
            power = int(term[0])
        except (TypeError, ValueError):
            problems.append(f"{context}[{index}]: invalid power")
            continue
        coefficient = _fraction(
            term[1], f"{context}[{index}].coefficient", problems
        )
        if power < 0 or power in powers or power <= previous_power:
            problems.append(
                f"{context}: powers must be unique and strictly increasing"
            )
        powers.add(power)
        previous_power = power
        total = _ci_add(total, _ci_scale(_ci_pow(s_box, power), coefficient))
    return total


def _parse_correction(
    text: str, name: str, problems: list[str]
) -> dict[int, Fraction]:
    match = CORRECTION_RE.fullmatch(text)
    if match is None:
        problems.append(f"{name}: finite correction is not canonical linear form")
        return {}
    slope = Fraction(int(match.group(1)), int(match.group(4)))
    constant_num = int(match.group(3))
    if match.group(2) == "-":
        constant_num = -constant_num
    constant = Fraction(constant_num, int(match.group(4)))
    return {1: slope, 0: constant}


def _expected_quantity_keys(name: str, segment_ids: Sequence[str]) -> set[str]:
    pairs = EXPECTED_CHART_PAIRS[name]
    keys = {
        f"probe:center:integral:A0({mass})" for mass in EXPECTED_A0_MASSES
    }
    keys.update(
        f"probe:center:integral:B0({m1},{m2})" for m1, m2 in pairs
    )
    keys.update(
        f"probe:center:derivative:B0p({m1},{m2})" for m1, m2 in pairs
    )
    keys.update(
        {
            "probe:center:inverse_propagator",
            "probe:center:inverse_propagator_derivative",
        }
    )
    for segment_id in segment_ids:
        keys.update(
            f"segment:{segment_id}:integral:B0({m1},{m2})"
            for m1, m2 in pairs
        )
        keys.update(
            f"segment:{segment_id}:derivative:B0p({m1},{m2})"
            for m1, m2 in pairs
        )
    return keys


def _check_source_digests(
    receipt: Mapping[str, Any],
    vector_path: Path,
    source_paths: Mapping[str, Path],
    problems: list[str],
) -> None:
    pins = receipt["pins"]
    all_paths = {"vector_blocks_sha256": vector_path, **source_paths}
    for key, path in all_paths.items():
        try:
            actual = sha256_bytes(path.read_bytes())
        except OSError as exc:
            problems.append(f"pin {key}: cannot read {path}: {exc}")
            continue
        if pins[key] != actual:
            problems.append(f"pin {key} does not match {path.name}")
    try:
        vector_raw = vector_path.read_bytes()
        vector = json.loads(vector_raw.decode("utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        problems.append(f"vector source cannot be decoded: {exc}")
        return
    if vector.get("schema") != "fj_direct_vector_blocks.v1":
        problems.append("vector source schema is not fj_direct_vector_blocks.v1")
    if vector.get("target") != "FJ_DIRECT_1":
        problems.append("vector source target is not FJ_DIRECT_1")
    if vector.get("units") != EXPECTED_VECTOR_UNITS:
        problems.append("vector source units drifted")
    source_subject = receipt["source_subject"]
    expected_subject = {
        "artifact_role": "direct_FJ_vector_blocks_input",
        "relative_path": "outputs/fj_direct_vector_blocks.json",
        "schema": vector.get("schema"),
        "target": vector.get("target"),
        "units": vector.get("units"),
        "bytes": len(vector_raw),
        "sha256": sha256_bytes(vector_raw),
    }
    if source_subject != expected_subject:
        problems.append("source_subject does not match the resolved vector artifact")
    if receipt["pins"]["vector_blocks_sha256"] != source_subject["sha256"]:
        problems.append("source_subject hash and vector pin disagree")


def _check_corrections(receipt: Mapping[str, Any], problems: list[str]) -> None:
    corrections = receipt["dimensional_prefactor_finite_correction"]
    for name in ("W", "Z"):
        text = corrections[name]
        if text != EXPECTED_CORRECTION_TEXT[name]:
            problems.append(f"{name}: finite correction text drifted")
        normalized = _parse_correction(text, name, problems)
        if normalized != EXPECTED_CORRECTION_COEFFICIENTS[name]:
            problems.append(f"{name}: normalized finite correction drifted")


def _check_loop_charts(
    name: str,
    interior: Mapping[str, Any],
    box: Mapping[str, tuple[Fraction, Fraction]],
    context: str,
    problems: list[str],
) -> bool:
    charts = interior["loop_charts"]
    actual_pairs = [(row["m1"], row["m2"]) for row in charts]
    if len(actual_pairs) != len(set(actual_pairs)):
        problems.append(f"{context}: duplicate loop chart")
    if set(actual_pairs) != EXPECTED_CHART_PAIRS[name]:
        missing = EXPECTED_CHART_PAIRS[name] - set(actual_pairs)
        extra = set(actual_pairs) - EXPECTED_CHART_PAIRS[name]
        problems.append(
            f"{context}: loop chart keyset mismatch; "
            f"missing={sorted(missing)}, extra={sorted(extra)}"
        )
    all_holomorphic = True
    re_lo, re_hi = box["re"]
    im_lo, _ = box["im"]
    for index, row in enumerate(charts):
        row_context = f"{context}.loop_charts[{index}]"
        m1 = _fraction(row["m1"], f"{row_context}.m1", problems)
        m2 = _fraction(row["m2"], f"{row_context}.m2", problems)
        if m1 == 0 and m2 == 0:
            expected_chart = "both_masses_zero"
            expected_certificate = (
                "log(s/mu2) only; box in open upper half plane"
            )
            expected = im_lo > 0
        elif m1 == 0 or m2 == 0:
            expected_chart = "one_mass_zero"
            expected_sign = "negative"
            expected_certificate = (
                "explicit root chart; both root-chart logarithm "
                f"arguments have strictly {expected_sign} imaginary part "
                "on the box"
            )
            if row.get("root_chart_log_argument_imaginary_sign") != expected_sign:
                problems.append(
                    f"{row_context}: root-chart imaginary sign does not replay"
                )
            expected = im_lo > 0 and re_lo > 0
        else:
            expected_chart = "two_massive"
            expected_certificate = (
                "discriminant imaginary part nonvanishing off the "
                "discriminant-symmetry line Re(s)=m1+m2 (not the "
                "physical threshold); Feynman roots never real because "
                "the quadratic equals m1 at x=0 and m2 at x=1"
            )
            line = m1 + m2
            outside = line < re_lo or line > re_hi
            if _fraction(
                row.get("discriminant_symmetry_line"),
                f"{row_context}.discriminant_symmetry_line",
                problems,
            ) != line:
                problems.append(f"{row_context}: symmetry line does not replay")
            if (
                bool(
                    row.get(
                        "discriminant_symmetry_line_outside_box_re_range"
                    )
                )
                != outside
            ):
                problems.append(
                    f"{row_context}: symmetry-line exclusion does not replay"
                )
            expected = outside and im_lo > 0 and re_lo > 0
        if row["chart"] != expected_chart:
            problems.append(f"{row_context}: chart classification drifted")
        if row["certificate"] != expected_certificate:
            problems.append(f"{row_context}: chart certificate text drifted")
        if bool(row["holomorphic"]) != expected:
            problems.append(f"{row_context}: holomorphy verdict does not replay")
        all_holomorphic = all_holomorphic and expected
    return all_holomorphic and set(actual_pairs) == EXPECTED_CHART_PAIRS[name]


def _check_denominators(
    name: str,
    interior: Mapping[str, Any],
    box: Mapping[str, tuple[Fraction, Fraction]],
    context: str,
    problems: list[str],
) -> tuple[
    bool,
    tuple[tuple[str, tuple[tuple[int, str], ...]], ...],
    dict[str, ComplexInterval],
]:
    witnesses = interior["coefficient_denominator_witnesses"]
    if len(witnesses) != EXPECTED_DENOMINATOR_COUNTS[name]:
        problems.append(
            f"{context}: expected {EXPECTED_DENOMINATOR_COUNTS[name]} "
            f"coefficient denominators, got {len(witnesses)}"
        )
    identities: list[tuple[str, tuple[tuple[int, str], ...]]] = []
    enclosures: dict[str, ComplexInterval] = {}
    all_excluded = True
    previous_id = ""
    for index, witness in enumerate(witnesses):
        witness_context = f"{context}.coefficient_denominator_witnesses[{index}]"
        coefficients = witness["coefficients"]
        expected_id = canonical_digest(coefficients)
        denominator_id = witness["denominator_id"]
        if denominator_id != expected_id:
            problems.append(f"{witness_context}: denominator id does not replay")
        if previous_id and denominator_id <= previous_id:
            problems.append(
                f"{context}: denominator witnesses are not strictly id-sorted"
            )
        previous_id = denominator_id
        coefficient_identity = tuple(
            (int(term[0]), str(term[1])) for term in coefficients
        )
        identities.append((denominator_id, coefficient_identity))
        exact_enclosure = _evaluate_polynomial(
            coefficients, box, f"{witness_context}.coefficients", problems
        )
        recorded = _complex_interval(
            witness["enclosure"], f"{witness_context}.enclosure", problems
        )
        enclosures[denominator_id] = recorded
        if not _ci_nested(recorded, exact_enclosure):
            problems.append(
                f"{witness_context}: recorded enclosure misses exact "
                "rational interval evaluation"
            )
        computed_lower = _abs2_lower(recorded)
        declared_lower = _fraction(
            witness["zero_exclusion_abs2_lower"],
            f"{witness_context}.zero_exclusion_abs2_lower",
            problems,
        )
        exclusion = (
            not _ci_contains_zero(recorded)
            and declared_lower > 0
            and declared_lower <= computed_lower
        )
        if bool(witness["excludes_zero"]) != exclusion:
            problems.append(
                f"{witness_context}: coefficient-denominator exclusion "
                "does not follow from its witness"
            )
        all_excluded = all_excluded and exclusion
    if len({identity[0] for identity in identities}) != len(identities):
        problems.append(f"{context}: duplicate denominator id")
        all_excluded = False
    identity_payload = [
        {
            "denominator_id": denominator_id,
            "coefficients": [
                [power, coefficient]
                for power, coefficient in coefficient_identity
            ],
        }
        for denominator_id, coefficient_identity in identities
    ]
    if canonical_digest(identity_payload) != EXPECTED_DENOMINATOR_IDENTITY_DIGESTS[
        name
    ]:
        problems.append(
            f"{context}: denominator identity set is not bound to the "
            "frozen vector subject"
        )
        all_excluded = False
    if (
        bool(interior["coefficient_denominators_exclude_zero_on_box"])
        != all_excluded
    ):
        problems.append(
            f"{context}: coefficient-denominator aggregate does not replay"
        )
    return all_excluded, tuple(identities), enclosures


def _partition_payload(
    evidence: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    return [
        {
            "segment_id": segment["segment_id"],
            "start": segment["start"],
            "end": segment["end"],
        }
        for segment in evidence
    ]


def _is_power_of_two(value: int) -> bool:
    return value > 0 and value & (value - 1) == 0


def _segment_edge_and_depth(
    start: Point,
    end: Point,
    box: Mapping[str, tuple[Fraction, Fraction]],
) -> tuple[int, int] | None:
    """Classify one CCW box-edge segment and replay its dyadic depth.

    The adaptive producer starts from eight equal cells on each edge and
    only bisects cells.  In the edge's positive traversal coordinate, a
    depth-``d`` leaf therefore has length ``edge_length / (8 * 2**d)``
    and begins at an integer multiple of that length.
    """

    re_lo, re_hi = box["re"]
    im_lo, im_hi = box["im"]
    if (
        start[1] == end[1] == im_lo
        and re_lo <= start[0] < end[0] <= re_hi
    ):
        edge = 0
        edge_length = re_hi - re_lo
        start_offset = start[0] - re_lo
        segment_length = end[0] - start[0]
    elif (
        start[0] == end[0] == re_hi
        and im_lo <= start[1] < end[1] <= im_hi
    ):
        edge = 1
        edge_length = im_hi - im_lo
        start_offset = start[1] - im_lo
        segment_length = end[1] - start[1]
    elif (
        start[1] == end[1] == im_hi
        and re_lo <= end[0] < start[0] <= re_hi
    ):
        edge = 2
        edge_length = re_hi - re_lo
        start_offset = re_hi - start[0]
        segment_length = start[0] - end[0]
    elif (
        start[0] == end[0] == re_lo
        and im_lo <= end[1] < start[1] <= im_hi
    ):
        edge = 3
        edge_length = im_hi - im_lo
        start_offset = im_hi - start[1]
        segment_length = start[1] - end[1]
    else:
        return None

    depth_ratio = edge_length / (
        INITIAL_SEGMENTS_PER_EDGE * segment_length
    )
    if depth_ratio.denominator != 1:
        return None
    subdivision_factor = depth_ratio.numerator
    if not _is_power_of_two(subdivision_factor):
        return None
    depth = subdivision_factor.bit_length() - 1
    if (start_offset / segment_length).denominator != 1:
        return None
    return edge, depth


def _check_partition_coverage(
    name: str,
    evidence: Sequence[Mapping[str, Any]],
    context: str,
    problems: list[str],
) -> tuple[tuple[tuple[str, Point, Point], ...], int]:
    box = BOXES[name]
    re_lo, re_hi = box["re"]
    im_lo, im_hi = box["im"]
    parsed: list[tuple[str, Point, Point]] = []
    edge_depths: list[tuple[int, int] | None] = []
    ids: set[str] = set()
    for index, segment in enumerate(evidence):
        segment_context = f"{context}.segment_evidence[{index}]"
        segment_id = segment["segment_id"]
        start = _point(segment["start"], f"{segment_context}.start", problems)
        end = _point(segment["end"], f"{segment_context}.end", problems)
        midpoint = _point(
            segment["midpoint"], f"{segment_context}.midpoint", problems
        )
        expected_midpoint = (
            (start[0] + end[0]) / 2,
            (start[1] + end[1]) / 2,
        )
        if midpoint != expected_midpoint:
            problems.append(f"{segment_context}: midpoint does not replay")
        expected_segment_id = canonical_digest(
            {"start": segment["start"], "end": segment["end"]}
        )
        if segment_id != expected_segment_id:
            problems.append(f"{segment_context}: segment id does not replay")
        if segment_id in ids:
            problems.append(f"{context}: duplicate segment id {segment_id}")
        ids.add(segment_id)
        if start == end:
            problems.append(f"{segment_context}: zero-length segment")
        parsed.append((segment_id, start, end))
        edge_depths.append(_segment_edge_and_depth(start, end, box))
    if not parsed:
        problems.append(f"{context}: empty boundary partition")
        return tuple(), 0
    if parsed[0][1] != (re_lo, im_lo):
        problems.append(f"{context}: partition does not start at lower-left")
    edge_groups: list[list[tuple[Point, Point]]] = [[], [], [], []]
    previous_edge = -1
    for index, (_, start, end) in enumerate(parsed):
        if index and parsed[index - 1][2] != start:
            problems.append(f"{context}: partition is not contiguous at {index}")
        classification = edge_depths[index]
        if classification is None:
            problems.append(
                f"{context}: segment {index} is not a counter-clockwise "
                "dyadic child of the declared 8-per-edge box boundary"
            )
            continue
        edge, depth = classification
        if depth > MAX_SUBDIVISION_DEPTH:
            problems.append(
                f"{context}: segment {index} reaches subdivision depth "
                f"{depth}, above the cap {MAX_SUBDIVISION_DEPTH}"
            )
        if edge < previous_edge:
            problems.append(
                f"{context}: partition is not exactly one CCW traversal "
                "(an edge repeats after a later edge)"
            )
        previous_edge = max(previous_edge, edge)
        edge_groups[edge].append((start, end))
    if parsed[-1][2] != parsed[0][1]:
        problems.append(f"{context}: partition does not close")
    edge_corners = (
        ((re_lo, im_lo), (re_hi, im_lo)),
        ((re_hi, im_lo), (re_hi, im_hi)),
        ((re_hi, im_hi), (re_lo, im_hi)),
        ((re_lo, im_hi), (re_lo, im_lo)),
    )
    for edge, (expected_start, expected_end) in enumerate(edge_corners):
        group = edge_groups[edge]
        if not group:
            problems.append(
                f"{context}: partition omits edge {edge} of the box"
            )
            continue
        if group[0][0] != expected_start or group[-1][1] != expected_end:
            problems.append(
                f"{context}: edge {edge} is not covered exactly once "
                "from corner to corner"
            )
    actual_max_depth = max(
        (
            classification[1]
            for classification in edge_depths
            if classification is not None
        ),
        default=0,
    )
    return tuple(parsed), actual_max_depth


def _check_boundary(
    name: str,
    row: Mapping[str, Any],
    context: str,
    problems: list[str],
) -> tuple[
    bool,
    tuple[tuple[str, Point, Point], ...],
    dict[str, dict[str, ComplexInterval | RealInterval]],
    RealInterval,
]:
    boundary = row["boundary_winding"]
    precision_bits = int(row["precision_bits"])
    evidence = boundary["segment_evidence"]
    if boundary["segments"] != len(evidence):
        problems.append(f"{context}: segment count does not match evidence")
    if len(evidence) != EXPECTED_SEGMENT_COUNTS[name]:
        problems.append(
            f"{context}: expected {EXPECTED_SEGMENT_COUNTS[name]} segments, "
            f"got {len(evidence)}"
        )
    partition, actual_max_depth = _check_partition_coverage(
        name, evidence, context, problems
    )
    recorded_max_depth = boundary["max_depth_used"]
    expected_max_depth = actual_max_depth if precision_bits == 128 else 0
    if recorded_max_depth != expected_max_depth:
        problems.append(
            f"{context}: max_depth_used={recorded_max_depth} does not replay "
            f"from the partition (expected {expected_max_depth})"
        )
    expected_partition_digest = canonical_digest(_partition_payload(evidence))
    if boundary["partition_sha256"] != expected_partition_digest:
        problems.append(f"{context}: partition digest does not replay")
    all_segments = True
    segment_values: dict[str, dict[str, ComplexInterval | RealInterval]] = {}
    summed_increment: RealInterval = (Fraction(0), Fraction(0))
    for index, segment in enumerate(evidence):
        segment_context = f"{context}.segment_evidence[{index}]"
        segment_id = segment["segment_id"]
        image = _complex_interval(
            segment["image"], f"{segment_context}.image", problems
        )
        center = _complex_interval(
            segment["center_value"],
            f"{segment_context}.center_value",
            problems,
        )
        derivative = _complex_interval(
            segment["derivative_hull"],
            f"{segment_context}.derivative_hull",
            problems,
        )
        offset = _complex_interval(
            segment["offset"], f"{segment_context}.offset", problems
        )
        rotated = _complex_interval(
            segment["rotated_image"],
            f"{segment_context}.rotated_image",
            problems,
        )
        ratio = _complex_interval(
            segment["endpoint_ratio"],
            f"{segment_context}.endpoint_ratio",
            problems,
        )
        start_value = _complex_interval(
            segment["start_value"],
            f"{segment_context}.start_value",
            problems,
        )
        end_value = _complex_interval(
            segment["end_value"],
            f"{segment_context}.end_value",
            problems,
        )
        rotated_arg = _real_interval(
            segment["rotated_argument"],
            f"{segment_context}.rotated_argument",
            problems,
        )
        increment = _real_interval(
            segment["endpoint_increment"],
            f"{segment_context}.endpoint_increment",
            problems,
        )
        image_lower = _fraction(
            segment["image_zero_exclusion_abs2_lower"],
            f"{segment_context}.image_zero_exclusion_abs2_lower",
            problems,
        )
        rotated_lower = _fraction(
            segment["rotated_zero_exclusion_abs2_lower"],
            f"{segment_context}.rotated_zero_exclusion_abs2_lower",
            problems,
        )
        image_excluded = (
            image_lower > 0
            and image_lower <= _abs2_lower(image)
            and not _ci_contains_zero(image)
        )
        rotated_excluded = (
            rotated_lower > 0
            and rotated_lower <= _abs2_lower(rotated)
            and not _ci_contains_zero(rotated)
        )
        start_point = _point(
            segment["start"], f"{segment_context}.start", problems
        )
        end_point = _point(
            segment["end"], f"{segment_context}.end", problems
        )
        midpoint = (
            (start_point[0] + end_point[0]) / 2,
            (start_point[1] + end_point[1]) / 2,
        )
        exact_offset: ComplexInterval = (
            (
                min(start_point[0], end_point[0]) - midpoint[0],
                max(start_point[0], end_point[0]) - midpoint[0],
            ),
            (
                min(start_point[1], end_point[1]) - midpoint[1],
                max(start_point[1], end_point[1]) - midpoint[1],
            ),
        )
        if not _ci_nested(offset, exact_offset):
            problems.append(
                f"{segment_context}: offset misses exact segment geometry"
            )
        expected_image = _ci_add(center, _ci_mul(derivative, offset))
        if not _ci_nested(image, expected_image):
            problems.append(
                f"{segment_context}: centered image containment does not replay"
            )
        if not _ci_nested(image, center):
            problems.append(
                f"{segment_context}: midpoint value escapes the segment image"
            )
        if not _ci_nested(image, start_value):
            problems.append(
                f"{segment_context}: start value escapes the segment image"
            )
        if not _ci_nested(image, end_value):
            problems.append(
                f"{segment_context}: end value escapes the segment image"
            )
        expected_rotated = _ci_mul(image, _ci_conjugate(center))
        if not _ci_nested(rotated, expected_rotated):
            problems.append(
                f"{segment_context}: rotated-image containment does not replay"
            )
        expected_ratio = _ci_mul(end_value, _ci_conjugate(start_value))
        if not _ci_nested(ratio, expected_ratio):
            problems.append(
                f"{segment_context}: endpoint-ratio containment does not replay"
            )
        replayed_rotated_arg = _atan2_replay(
            rotated,
            precision_bits,
            f"{segment_context}.rotated_argument",
            problems,
        )
        rotated_arg_replayed = _ri_nested(
            rotated_arg, replayed_rotated_arg
        )
        if not rotated_arg_replayed:
            problems.append(
                f"{segment_context}: rotated argument misses independent "
                "atan2 replay"
            )
        replayed_increment = _atan2_replay(
            ratio,
            precision_bits,
            f"{segment_context}.endpoint_increment",
            problems,
        )
        increment_replayed = _ri_nested(increment, replayed_increment)
        if not increment_replayed:
            problems.append(
                f"{segment_context}: endpoint increment misses independent "
                "atan2 replay"
            )
        width = rotated_arg[1] - rotated_arg[0]
        recorded_width = _fraction(
            segment["rotated_argument_width"],
            f"{segment_context}.rotated_argument_width",
            problems,
        )
        width_slack = _fraction(
            segment["argument_width_slack"],
            f"{segment_context}.argument_width_slack",
            problems,
        )
        expected_width_slack = ARG_WIDTH_GATE - width
        if recorded_width != width:
            problems.append(f"{segment_context}: rotated argument width drifted")
        if width_slack != expected_width_slack:
            problems.append(f"{segment_context}: argument-width slack drifted")
        increment_extent = max(abs(increment[0]), abs(increment[1]))
        increment_slack = _fraction(
            segment["endpoint_increment_slack"],
            f"{segment_context}.endpoint_increment_slack",
            problems,
        )
        expected_increment_slack = ARG_WIDTH_GATE - increment_extent
        if increment_slack != expected_increment_slack:
            problems.append(f"{segment_context}: endpoint-increment slack drifted")
        segment_ok = (
            image_excluded
            and rotated_excluded
            and not _ci_contains_zero(ratio)
            and rotated_arg_replayed
            and increment_replayed
            and width >= 0
            and width_slack > 0
            and increment_slack > 0
        )
        all_segments = all_segments and segment_ok
        summed_increment = _ri_add(summed_increment, increment)
        segment_values[segment_id] = {
            "center_value": center,
            "derivative_hull": derivative,
            "offset": offset,
            "image": image,
            "rotated_image": rotated,
            "start_value": start_value,
            "end_value": end_value,
            "endpoint_ratio": ratio,
            "rotated_argument": rotated_arg,
            "endpoint_increment": increment,
        }
    ordered_ids = [segment["segment_id"] for segment in evidence]
    chain_continuous = bool(ordered_ids)
    for index, segment_id in enumerate(ordered_ids):
        next_id = ordered_ids[(index + 1) % len(ordered_ids)]
        if (
            segment_values[segment_id]["end_value"]
            != segment_values[next_id]["start_value"]
        ):
            problems.append(
                f"{context}: endpoint-value chain is discontinuous after "
                f"segment {segment_id}"
            )
            chain_continuous = False
    all_segments = all_segments and chain_continuous
    total = _real_interval(
        boundary["total_variation_interval"],
        f"{context}.total_variation_interval",
        problems,
    )
    if not _ri_nested(total, summed_increment):
        problems.append(
            f"{context}: total variation does not enclose the sum of "
            "segment increments"
        )
    residual = _real_interval(
        boundary["winding_residual"],
        f"{context}.winding_residual",
        problems,
    )
    winding = boundary["winding"]
    if winding != 0:
        problems.append(f"{context}: principal-sheet receipt must record winding 0")
    if winding == 0 and not _ri_nested(residual, total):
        problems.append(
            f"{context}: winding-zero residual does not enclose total"
        )
    tolerance_slack = _fraction(
        boundary["winding_tolerance_slack"],
        f"{context}.winding_tolerance_slack",
        problems,
    )
    expected_tolerance_slack = WINDING_TOLERANCE - max(
        abs(residual[0]), abs(residual[1])
    )
    if tolerance_slack != expected_tolerance_slack:
        problems.append(f"{context}: winding-tolerance slack drifted")
    expected_certified = bool(
        all_segments
        and winding == 0
        and tolerance_slack > 0
        and len(evidence) > 0
    )
    if bool(boundary["certified"]) != expected_certified:
        problems.append(f"{context}: boundary verdict does not replay")
    if expected_certified and boundary["reason"] is not None:
        problems.append(f"{context}: certified boundary has a failure reason")
    return expected_certified, partition, segment_values, total


def _check_quantity_keys(
    name: str,
    partition: Sequence[tuple[str, Point, Point]],
    quantities: Mapping[str, Any],
    context: str,
    problems: list[str],
) -> dict[str, ComplexInterval]:
    expected = _expected_quantity_keys(name, [item[0] for item in partition])
    actual = set(quantities)
    if len(actual) != EXPECTED_QUANTITY_COUNTS[name]:
        problems.append(
            f"{context}: expected {EXPECTED_QUANTITY_COUNTS[name]} quantities, "
            f"got {len(actual)}"
        )
    if actual != expected:
        missing = expected - actual
        extra = actual - expected
        problems.append(
            f"{context}: quantity keyset mismatch; "
            f"missing={len(missing)}, extra={len(extra)}"
        )
    return {
        key: _complex_interval(value, f"{context}.{key}", problems)
        for key, value in quantities.items()
    }


def _check_particle_rows(
    receipt: Mapping[str, Any], name: str, problems: list[str]
) -> tuple[bool, dict[str, Any]]:
    rows = receipt["results"][name]
    box = BOXES[name]
    row_evidence: dict[str, Any] = {}
    all_rows = True
    reference_denominators = None
    reference_charts = None
    for precision_text in PRECISIONS:
        row = rows[precision_text]
        context = f"results.{name}.{precision_text}"
        if row["particle"] != name:
            problems.append(f"{context}: particle label does not match parent")
        if row["precision_bits"] != int(precision_text):
            problems.append(f"{context}: precision label does not match key")
        recorded_box = row["box"]
        parsed_box = {
            "re": tuple(
                _fraction(v, f"{context}.box.re", problems)
                for v in recorded_box["re"]
            ),
            "im": tuple(
                _fraction(v, f"{context}.box.im", problems)
                for v in recorded_box["im"]
            ),
        }
        if parsed_box != box or recorded_box["half_plane"] != "upper":
            problems.append(f"{context}: declared box drifted")
        if (
            _fraction(row["tree_mass_sq"], f"{context}.tree_mass_sq", problems)
            != TREE_MASSES[name]
        ):
            problems.append(f"{context}: tree mass drifted")
        interior = row["interior_holomorphy"]
        expected_base = {
            "box_in_open_upper_half_plane": box["im"][0] > 0,
            "box_real_part_positive": box["re"][0] > 0,
        }
        if interior["base_facts"] != expected_base:
            problems.append(f"{context}: interior base facts do not replay")
        (
            denominators_ok,
            denominator_identity,
            denominator_enclosures,
        ) = _check_denominators(
            name, interior, box, f"{context}.interior_holomorphy", problems
        )
        charts_ok = _check_loop_charts(
            name, interior, box, f"{context}.interior_holomorphy", problems
        )
        chart_identity = tuple(
            (chart["m1"], chart["m2"], chart["chart"])
            for chart in interior["loop_charts"]
        )
        if reference_denominators is None:
            reference_denominators = denominator_identity
            reference_charts = chart_identity
        else:
            if denominator_identity != reference_denominators:
                problems.append(
                    f"{context}: denominator witness keyset/order differs "
                    "across precision"
                )
            if chart_identity != reference_charts:
                problems.append(
                    f"{context}: loop chart keyset/order differs across precision"
                )
        expected_holomorphic = bool(
            all(expected_base.values()) and denominators_ok and charts_ok
        )
        if bool(interior["holomorphic_on_box"]) != expected_holomorphic:
            problems.append(f"{context}: interior aggregate does not replay")
        boundary_ok, partition, segment_values, total = _check_boundary(
            name, row, f"{context}.boundary_winding", problems
        )
        expected_partition_kind = (
            "adaptive" if precision_text == "128" else "replayed_base_partition"
        )
        if row["boundary_winding"]["partition"] != expected_partition_kind:
            problems.append(
                f"{context}: partition kind is not base-adaptive then replayed"
            )
        quantities = _check_quantity_keys(
            name,
            partition,
            row["quantity_enclosures"],
            f"{context}.quantity_enclosures",
            problems,
        )
        expected_row = expected_holomorphic and boundary_ok
        if bool(row["zero_exclusion_certified"]) != expected_row:
            problems.append(f"{context}: row zero-exclusion verdict does not replay")
        all_rows = all_rows and expected_row
        row_evidence[precision_text] = {
            "partition": partition,
            "segments": segment_values,
            "total": total,
            "quantities": quantities,
            "denominators": denominator_enclosures,
        }
    return all_rows, row_evidence


def _check_nesting(
    receipt: Mapping[str, Any],
    name: str,
    evidence: Mapping[str, Any],
    problems: list[str],
) -> bool:
    context = f"precision_nesting.{name}"
    summary = receipt["precision_nesting"][name]
    rows = [evidence[precision] for precision in PRECISIONS]
    partitions_match = all(
        rows[index]["partition"] == rows[0]["partition"]
        for index in range(1, len(rows))
    )
    totals_nested = all(
        _ri_nested(rows[index]["total"], rows[index + 1]["total"])
        for index in range(len(rows) - 1)
    )
    quantity_key_sets_match = all(
        set(row["quantities"]) == set(rows[0]["quantities"]) for row in rows[1:]
    )
    quantity_keys = set(rows[0]["quantities"])
    per_quantity = {
        key: all(
            key in rows[index + 1]["quantities"]
            and _ci_nested(
                rows[index]["quantities"][key],
                rows[index + 1]["quantities"][key],
            )
            for index in range(len(rows) - 1)
        )
        for key in quantity_keys
    }
    quantities_nested = bool(
        quantity_key_sets_match and per_quantity and all(per_quantity.values())
    )
    segment_ids = [item[0] for item in rows[0]["partition"]]

    def segment_ci_nested(field: str) -> bool:
        return bool(
            partitions_match
            and all(
                _ci_nested(
                    rows[index]["segments"][segment_id][field],
                    rows[index + 1]["segments"][segment_id][field],
                )
                for index in range(len(rows) - 1)
                for segment_id in segment_ids
            )
        )

    def segment_ri_nested(field: str) -> bool:
        return bool(
            partitions_match
            and all(
                _ri_nested(
                    rows[index]["segments"][segment_id][field],
                    rows[index + 1]["segments"][segment_id][field],
                )
                for index in range(len(rows) - 1)
                for segment_id in segment_ids
            )
        )

    center_values_nested = segment_ci_nested("center_value")
    derivative_hulls_nested = segment_ci_nested("derivative_hull")
    offsets_nested = segment_ci_nested("offset")
    images_nested = segment_ci_nested("image")
    rotated_images_nested = segment_ci_nested("rotated_image")
    rotated_arguments_nested = segment_ri_nested("rotated_argument")
    start_values_nested = segment_ci_nested("start_value")
    end_values_nested = segment_ci_nested("end_value")
    endpoint_values_nested = start_values_nested and end_values_nested
    endpoint_ratios_nested = segment_ci_nested("endpoint_ratio")
    endpoint_increments_nested = segment_ri_nested("endpoint_increment")
    segment_all = bool(
        partitions_match
        and center_values_nested
        and derivative_hulls_nested
        and offsets_nested
        and images_nested
        and rotated_images_nested
        and rotated_arguments_nested
        and endpoint_values_nested
        and endpoint_ratios_nested
        and endpoint_increments_nested
    )
    denominator_key_sets_match = all(
        set(row["denominators"]) == set(rows[0]["denominators"])
        for row in rows[1:]
    ) and bool(rows[0]["denominators"])
    denominators_nested = bool(
        denominator_key_sets_match
        and all(
            _ci_nested(
                rows[index]["denominators"][denominator_id],
                rows[index + 1]["denominators"][denominator_id],
            )
            for index in range(len(rows) - 1)
            for denominator_id in rows[0]["denominators"]
        )
    )
    expected_per_quantity_keys = set(rows[0]["quantities"])
    recorded_per_quantity = summary["per_quantity_probe_nesting"]
    if set(recorded_per_quantity) != expected_per_quantity_keys:
        problems.append(f"{context}: per-quantity summary keyset mismatch")
    for key in expected_per_quantity_keys & set(recorded_per_quantity):
        if bool(recorded_per_quantity[key]) != per_quantity[key]:
            problems.append(f"{context}: nesting summary drifted for {key}")
    expected_summary = {
        "enclosures_nested_with_precision": totals_nested,
        "comparison": (
            "exact rational comparison of serialized binary endpoints"
        ),
        "partition_ids_match": partitions_match,
        "quantity_key_sets_match": quantity_key_sets_match,
        "quantity_enclosures_all_nested": quantities_nested,
    }
    for key, expected in expected_summary.items():
        if summary[key] != expected:
            problems.append(f"{context}: {key} does not replay")
    segment_summary = summary["per_segment_enclosure_nesting"]
    expected_segment_summary = {
        "segments": len(segment_ids),
        "ladders_compared": 2,
        "partition_ids_match": partitions_match,
        "center_values_nested": center_values_nested,
        "derivative_hulls_nested": derivative_hulls_nested,
        "offsets_nested": offsets_nested,
        "images_nested": images_nested,
        "rotated_images_nested": rotated_images_nested,
        "rotated_arguments_nested": rotated_arguments_nested,
        "endpoint_values_nested": endpoint_values_nested,
        "endpoint_ratios_nested": endpoint_ratios_nested,
        "endpoint_increments_nested": endpoint_increments_nested,
        "all_nested": segment_all,
    }
    for key, expected in expected_segment_summary.items():
        if segment_summary[key] != expected:
            problems.append(f"{context}: segment summary {key} does not replay")
    expected_probe = [
        str((BOXES[name]["re"][0] + BOXES[name]["re"][1]) / 2),
        str((BOXES[name]["im"][0] + BOXES[name]["im"][1]) / 2),
    ]
    if summary["probe_point"] != expected_probe:
        problems.append(f"{context}: probe point drifted")
    denominator_summary = summary["coefficient_denominator_nesting"]
    expected_denominator_summary = {
        "records": len(rows[0]["denominators"]),
        "key_sets_match": denominator_key_sets_match,
        "all_nested": denominators_nested,
    }
    if denominator_summary != expected_denominator_summary:
        problems.append(
            f"{context}: coefficient-denominator nesting summary does not replay"
        )
    return bool(
        totals_nested
        and partitions_match
        and quantities_nested
        and denominators_nested
        and segment_all
    )


def _check_claims(
    receipt: Mapping[str, Any], all_certified: bool, problems: list[str]
) -> None:
    scope = receipt["acceptance_scope"]
    expected_scope = {
        "auxiliary_principal_sheet_zero_exclusion_only": True,
        "independent_numerical_replay_certified": False,
        "engine_inverse_propagator_convention": (
            "G(s)=s-m_tree^2-Pi_engine(s)"
        ),
        "theorem_self_energy_sign_bridge_certified": False,
        "coefficient_denominator_witness_is_laurent_denominator": False,
        "root_enclosure_certified": False,
        "laurent_denominator_certified": False,
        "issue_593_precision_ladder_row_satisfied": False,
        "issue_593_independent_third_verifier_row_satisfied": False,
        "issue_593_root_laurent_row_satisfied": False,
        "issue_593_full_acceptance_satisfied": False,
    }
    if scope != expected_scope:
        problems.append("acceptance scope overstates or changes the receipt")
    promotion = receipt["promotion"]
    for key in (
        "complex_ball_certified",
        "sheet_certified_on_declared_boxes",
        "principal_sheet_zero_exclusion_certified",
        "root_count_certified_on_declared_boxes",
    ):
        if bool(promotion[key]) != all_certified:
            problems.append(f"promotion.{key} does not follow from evidence")
    for key in (
        "pole_enclosure_certified",
        "second_sheet_certified",
        "laurent_residue_certified",
        "bmhv_restoration_certified",
        "physical_current_claim",
        "oph_native",
        "unit_claim",
    ):
        if promotion[key] is not False:
            problems.append(f"promotion.{key} must remain false")
    expected_status = (
        "PRINCIPAL_SHEET_ZERO_EXCLUSION_CERTIFIED"
        if all_certified
        else "CERTIFICATION_INCOMPLETE"
    )
    if receipt["status"] != expected_status:
        problems.append("top-level status does not follow from evidence")


def check(
    payload: Mapping[str, Any] | None = None,
    *,
    artifact_path: Path = ARTIFACT_PATH,
    vector_path: Path = VECTOR_PATH,
    schema_path: Path = SCHEMA_PATH,
    source_paths: Mapping[str, Path] | None = None,
) -> dict[str, Any]:
    """Validate one v3 receipt and return a mutation-test-friendly verdict."""

    if payload is None:
        try:
            payload = json.loads(artifact_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            return {
                "status": FAIL_STATUS,
                "evidence_contract_valid": False,
                "independent_physics_reevaluation_performed": False,
                "problems": [f"artifact cannot be decoded: {exc}"],
            }
    receipt = copy.deepcopy(dict(payload))
    try:
        schema = json.loads(schema_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return {
            "status": FAIL_STATUS,
            "evidence_contract_valid": False,
            "independent_physics_reevaluation_performed": False,
            "problems": [f"schema cannot be decoded: {exc}"],
        }
    schema_errors = sorted(
        Draft202012Validator(schema).iter_errors(receipt),
        key=lambda error: list(error.absolute_path),
    )
    problems = [
        "schema at "
        + (".".join(str(part) for part in error.absolute_path) or "<root>")
        + f": {error.message}"
        for error in schema_errors
    ]
    if problems:
        return {
            "status": FAIL_STATUS,
            "evidence_contract_valid": False,
            "independent_physics_reevaluation_performed": False,
            "problems": problems,
        }
    declared_digest = receipt["receipt_sha256"]
    recomputed_digest = canonical_digest(
        {key: value for key, value in receipt.items() if key != "receipt_sha256"}
    )
    if declared_digest != recomputed_digest:
        problems.append("receipt digest mismatch")
    if receipt["fixture"] != EXPECTED_FIXTURE:
        problems.append("frozen rational fixture drifted")
    if receipt["ckm_fixture"] != EXPECTED_CKM:
        problems.append("frozen CKM fixture is not the exact identity matrix")
    if source_paths is None:
        source_paths = {
            "diagnostic_module_sha256": ROOT
            / "producers"
            / "wz_pole_receipts.py",
            "interval_module_sha256": ROOT
            / "producers"
            / "complex_interval.py",
            "producer_module_sha256": ROOT
            / "producers"
            / "certified_wz_contours.py",
        }
    _check_source_digests(
        receipt, vector_path, source_paths, problems
    )
    _check_corrections(receipt, problems)
    all_certified = True
    evidence_by_particle = {}
    for name in ("W", "Z"):
        rows_ok, evidence = _check_particle_rows(receipt, name, problems)
        evidence_by_particle[name] = evidence
        nesting_ok = _check_nesting(receipt, name, evidence, problems)
        all_certified = all_certified and rows_ok and nesting_ok
    _check_claims(receipt, all_certified, problems)
    if not all_certified:
        problems.append(
            "committed principal-sheet zero-exclusion certificate is incomplete"
        )
    return {
        "schema": "oph.certified_wz_contours.evidence_check.v1",
        "status": PASS_STATUS if not problems else FAIL_STATUS,
        "evidence_contract_valid": not problems,
        "principal_sheet_zero_exclusion_evidence_valid": bool(
            not problems and all_certified
        ),
        "independent_physics_reevaluation_performed": False,
        "issue_593_precision_ladder_row_satisfied": False,
        "issue_593_independent_third_verifier_row_satisfied": False,
        "issue_593_root_laurent_row_satisfied": False,
        "issue_593_full_acceptance_satisfied": False,
        "problems": problems,
    }


def validate_receipt(
    receipt: Mapping[str, Any],
    **kwargs: Any,
) -> dict[str, Any]:
    """Validate an in-memory receipt or raise ``ReceiptValidationError``.

    This API is intended for adversarial tests: callers may mutate a deep
    copy, recompute its self-digest, and establish that a semantic check—not
    merely the digest gate—rejects the forged evidence.
    """

    verdict = check(receipt, **kwargs)
    if not verdict["evidence_contract_valid"]:
        raise ReceiptValidationError(verdict["problems"])
    return verdict


def main() -> int:
    try:
        receipt = json.loads(ARTIFACT_PATH.read_text(encoding="utf-8"))
        validate_receipt(receipt)
    except (OSError, json.JSONDecodeError) as exc:
        print(f"{FAIL_STATUS}: artifact cannot be decoded: {exc}", file=sys.stderr)
        return 1
    except ReceiptValidationError as exc:
        for problem in exc.problems:
            print(f"{FAIL_STATUS}: {problem}", file=sys.stderr)
        return 1
    print(PASS_STATUS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
