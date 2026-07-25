#!/usr/bin/env python3
"""Build the exact icosahedral residual-axis angle-spectrum receipt.

The construction starts from the standard real icosahedron with vertices

    (0, +/-1, +/-phi), (+/-1, +/-phi, 0), (+/-phi, 0, +/-1).

Its vertex, face-centre, and edge-midpoint lines are respectively the
five-fold, three-fold, and two-fold unoriented residual axes.  All geometry is
constructed in Q(sqrt(5)); the receipt records an exact value of cos(angle)^2
for every distinct acute angle and a rounded numerical angle for readability.

The no-go conclusion is intentionally narrow.  It excludes only direct
identification of the Cabibbo angle with an angle between these real
three-dimensional residual axes.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "icosahedral_axis_angle_spectrum.json"
VUS_FIXTURE = ROOT / "particles" / "data" / "pdg_2024_vus_kmu2_fixture.json"


@dataclass(frozen=True)
class Q5:
    """An exact element ``a + b*sqrt(5)`` of Q(sqrt(5))."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    @staticmethod
    def coerce(value: Q5 | Fraction | int) -> Q5:
        if isinstance(value, Q5):
            return value
        return Q5(Fraction(value), Fraction(0))

    def __add__(self, other: Q5 | Fraction | int) -> Q5:
        rhs = self.coerce(other)
        return Q5(self.a + rhs.a, self.b + rhs.b)

    def __radd__(self, other: Q5 | Fraction | int) -> Q5:
        return self + other

    def __sub__(self, other: Q5 | Fraction | int) -> Q5:
        rhs = self.coerce(other)
        return Q5(self.a - rhs.a, self.b - rhs.b)

    def __rsub__(self, other: Q5 | Fraction | int) -> Q5:
        return self.coerce(other) - self

    def __neg__(self) -> Q5:
        return Q5(-self.a, -self.b)

    def __mul__(self, other: Q5 | Fraction | int) -> Q5:
        rhs = self.coerce(other)
        return Q5(
            self.a * rhs.a + 5 * self.b * rhs.b,
            self.a * rhs.b + self.b * rhs.a,
        )

    def __rmul__(self, other: Q5 | Fraction | int) -> Q5:
        return self * other

    def __truediv__(self, other: Q5 | Fraction | int) -> Q5:
        rhs = self.coerce(other)
        denominator = rhs.a * rhs.a - 5 * rhs.b * rhs.b
        if denominator == 0:
            raise ZeroDivisionError("division by zero in Q(sqrt(5))")
        numerator = self * Q5(rhs.a, -rhs.b)
        return Q5(numerator.a / denominator, numerator.b / denominator)

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0

    def sign(self) -> int:
        """Return the exact sign without a floating-point comparison."""

        if self.is_zero():
            return 0
        if self.a == 0:
            return 1 if self.b > 0 else -1
        if self.b == 0:
            return 1 if self.a > 0 else -1
        if self.a > 0 and self.b > 0:
            return 1
        if self.a < 0 and self.b < 0:
            return -1
        comparison = self.a * self.a - 5 * self.b * self.b
        if comparison == 0:
            raise AssertionError("sqrt(5) would be rational")
        if self.a > 0:
            return 1 if comparison > 0 else -1
        return -1 if comparison > 0 else 1

    def to_float(self) -> float:
        return float(self.a) + float(self.b) * math.sqrt(5.0)


ZERO = Q5()
ONE = Q5(Fraction(1))
PHI = Q5(Fraction(1, 2), Fraction(1, 2))
Vector = tuple[Q5, Q5, Q5]


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _q5_record(value: Q5) -> dict[str, str]:
    a_text = _fraction_text(value.a)
    b_text = _fraction_text(value.b)
    return {
        "rational_coefficient": a_text,
        "sqrt5_coefficient": b_text,
        "expression": f"({a_text}) + ({b_text})*sqrt(5)",
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_compare_only_vus(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    coordinate = payload.get("coordinate")
    boundary = payload.get("claim_boundary")
    source = payload.get("source")
    if (
        payload.get("artifact") != "oph_pdg_2024_vus_kmu2_fixture"
        or payload.get("status") != "COMPARE_ONLY_HAND_TRANSCRIBED_REFERENCE"
        or not isinstance(coordinate, dict)
        or not isinstance(boundary, dict)
        or not isinstance(source, dict)
        or coordinate.get("name") != "abs_Vus"
        or coordinate.get("determination") != "Kmu2_decay_constant_ratio"
        or boundary.get("comparison_only") is not True
        or boundary.get("used_to_construct_or_select_axes") is not False
        or boundary.get("oph_fit_or_selection_input") is not False
    ):
        raise ValueError("fail closed: invalid compare-only |V_us| fixture boundary")
    value = float(str(coordinate.get("value")))
    uncertainty = float(str(coordinate.get("standard_uncertainty")))
    if not (0.0 < value < 1.0 and 0.0 < uncertainty < value):
        raise ValueError("fail closed: invalid |V_us| value or uncertainty")
    return {
        "value": value,
        "standard_uncertainty": uncertainty,
        "published_notation": str(coordinate.get("published_notation")),
        "determination": str(coordinate.get("determination")),
        "source": source,
        "fixture_path": "code/particles/data/pdg_2024_vus_kmu2_fixture.json",
        "fixture_sha256": _sha256(path),
    }


def _add(left: Vector, right: Vector) -> Vector:
    return tuple(x + y for x, y in zip(left, right, strict=True))  # type: ignore[return-value]


def _sub(left: Vector, right: Vector) -> Vector:
    return tuple(x - y for x, y in zip(left, right, strict=True))  # type: ignore[return-value]


def _neg(vector: Vector) -> Vector:
    return tuple(-entry for entry in vector)  # type: ignore[return-value]


def _dot(left: Vector, right: Vector) -> Q5:
    return sum((x * y for x, y in zip(left, right, strict=True)), ZERO)


def _norm_squared(vector: Vector) -> Q5:
    return _dot(vector, vector)


def _axis_key(vector: Vector) -> tuple[Fraction, ...]:
    return tuple(component for entry in vector for component in (entry.a, entry.b))


def _canonical_axis(vector: Vector) -> Vector:
    if all(entry.is_zero() for entry in vector):
        raise ValueError("the zero vector does not define an axis")
    for entry in vector:
        sign = entry.sign()
        if sign > 0:
            return vector
        if sign < 0:
            return _neg(vector)
    raise AssertionError("unreachable")


def _collapse_antipodes(vectors: Iterable[Vector]) -> list[Vector]:
    axes = {_canonical_axis(vector) for vector in vectors}
    return sorted(axes, key=_axis_key)


def build_icosahedral_axes() -> dict[str, list[Vector]]:
    """Construct the 6, 10, and 15 unoriented symmetry-axis families."""

    vertices: list[Vector] = []
    for first_sign in (-1, 1):
        for second_sign in (-1, 1):
            first = Q5(Fraction(first_sign))
            second_phi = Q5(Fraction(second_sign)) * PHI
            vertices.extend(
                [
                    (ZERO, first, second_phi),
                    (first, second_phi, ZERO),
                    (second_phi, ZERO, first),
                ]
            )
    if len(set(vertices)) != 12:
        raise AssertionError("the standard coordinate construction must emit 12 vertices")

    edge_pairs: list[tuple[int, int]] = []
    for left, right in itertools.combinations(range(len(vertices)), 2):
        if _norm_squared(_sub(vertices[left], vertices[right])) == Q5(Fraction(4)):
            edge_pairs.append((left, right))
    if len(edge_pairs) != 30:
        raise AssertionError(f"expected 30 icosahedron edges, found {len(edge_pairs)}")

    edge_set = {tuple(sorted(pair)) for pair in edge_pairs}
    faces: list[tuple[int, int, int]] = []
    for triple in itertools.combinations(range(len(vertices)), 3):
        if all(tuple(sorted(pair)) in edge_set for pair in itertools.combinations(triple, 2)):
            faces.append(triple)
    if len(faces) != 20:
        raise AssertionError(f"expected 20 triangular faces, found {len(faces)}")

    face_centres = [
        _add(_add(vertices[first], vertices[second]), vertices[third])
        for first, second, third in faces
    ]
    edge_midpoints = [_add(vertices[first], vertices[second]) for first, second in edge_pairs]

    axes = {
        "fivefold": _collapse_antipodes(vertices),
        "threefold": _collapse_antipodes(face_centres),
        "twofold": _collapse_antipodes(edge_midpoints),
    }
    expected = {"fivefold": 6, "threefold": 10, "twofold": 15}
    actual = {family: len(vectors) for family, vectors in axes.items()}
    if actual != expected:
        raise AssertionError(f"unexpected unoriented-axis counts: {actual}")
    return axes


def _cosine_squared(left: Vector, right: Vector) -> Q5:
    value = (_dot(left, right) * _dot(left, right)) / (
        _norm_squared(left) * _norm_squared(right)
    )
    numeric = value.to_float()
    if numeric < -1e-12 or numeric > 1.0 + 1e-12:
        raise AssertionError(f"cosine squared outside [0,1]: {numeric}")
    return value


def _angle_degrees(cosine_squared: Q5) -> float:
    numeric = min(1.0, max(0.0, cosine_squared.to_float()))
    return math.degrees(math.acos(math.sqrt(numeric)))


def _family_pairs(
    left: Sequence[Vector],
    right: Sequence[Vector],
    same_family: bool,
) -> Iterable[tuple[Vector, Vector]]:
    if same_family:
        return itertools.combinations(left, 2)
    return itertools.product(left, right)


def _spectrum(
    left: Sequence[Vector],
    right: Sequence[Vector],
    same_family: bool,
) -> list[dict[str, object]]:
    multiplicities: dict[Q5, int] = {}
    for left_axis, right_axis in _family_pairs(left, right, same_family):
        cosine_squared = _cosine_squared(left_axis, right_axis)
        multiplicities[cosine_squared] = multiplicities.get(cosine_squared, 0) + 1

    entries = []
    for cosine_squared, multiplicity in multiplicities.items():
        exact_record = _q5_record(cosine_squared)
        entries.append(
            {
                "acute_angle_degrees": round(_angle_degrees(cosine_squared), 12),
                "acute_angle_exact": (
                    f"acos(sqrt({exact_record['expression']}))"
                ),
                "cosine_squared_exact": exact_record,
                "multiplicity": multiplicity,
            }
        )
    return sorted(entries, key=lambda entry: float(entry["acute_angle_degrees"]))


def build_artifact(
    vus_fixture_path: Path = VUS_FIXTURE,
) -> dict[str, object]:
    axes = build_icosahedral_axes()
    compare_only_vus = _load_compare_only_vus(vus_fixture_path)
    family_specs = [
        ("fivefold_x_fivefold", "fivefold", "fivefold", True),
        ("fivefold_x_threefold", "fivefold", "threefold", False),
        ("fivefold_x_twofold", "fivefold", "twofold", False),
        ("threefold_x_threefold", "threefold", "threefold", True),
        ("threefold_x_twofold", "threefold", "twofold", False),
        ("twofold_x_twofold", "twofold", "twofold", True),
    ]
    spectra: dict[str, dict[str, object]] = {}
    all_angles: list[float] = []
    for row_name, left_name, right_name, same_family in family_specs:
        entries = _spectrum(axes[left_name], axes[right_name], same_family)
        pair_count = (
            len(axes[left_name]) * (len(axes[left_name]) - 1) // 2
            if same_family
            else len(axes[left_name]) * len(axes[right_name])
        )
        if sum(int(entry["multiplicity"]) for entry in entries) != pair_count:
            raise AssertionError(f"multiplicity mismatch in {row_name}")
        spectra[row_name] = {
            "left_family": left_name,
            "right_family": right_name,
            "unoriented_pair_count": pair_count,
            "distinct_acute_angle_count": len(entries),
            "angles": entries,
        }
        all_angles.extend(float(entry["acute_angle_degrees"]) for entry in entries)

    golden_cosine_squared = (PHI * PHI) / (PHI * PHI + ONE)
    five_by_two_exact = {
        (
            str(entry["cosine_squared_exact"]["rational_coefficient"]),
            str(entry["cosine_squared_exact"]["sqrt5_coefficient"]),
        )
        for entry in spectra["fivefold_x_twofold"]["angles"]  # type: ignore[index]
    }
    golden_key = (
        _fraction_text(golden_cosine_squared.a),
        _fraction_text(golden_cosine_squared.b),
    )
    if golden_key not in five_by_two_exact:
        raise AssertionError("arctan(1/phi) is absent from the five-fold x two-fold spectrum")

    golden_angle = math.degrees(math.atan(1.0 / PHI.to_float()))
    emitted_five_by_two_angles = [
        float(entry["acute_angle_degrees"])
        for entry in spectra["fivefold_x_twofold"]["angles"]  # type: ignore[index]
    ]
    if min(abs(angle - golden_angle) for angle in emitted_five_by_two_angles) > 1e-10:
        raise AssertionError("golden-angle numerical self-test failed")

    nonzero_angles = [angle for angle in all_angles if angle > 1e-12]
    minimum_nonzero = min(nonzero_angles)
    expected_minimum = 20.905157447889
    if abs(minimum_nonzero - expected_minimum) > 1e-10:
        raise AssertionError(
            f"unexpected minimum nonzero axis angle: {minimum_nonzero}"
        )

    compare_only_abs_vus = float(compare_only_vus["value"])
    compare_only_abs_vus_uncertainty = float(
        compare_only_vus["standard_uncertainty"]
    )
    compare_only_angle = math.degrees(math.asin(compare_only_abs_vus))
    one_sigma_lower = compare_only_abs_vus - compare_only_abs_vus_uncertainty
    one_sigma_upper = compare_only_abs_vus + compare_only_abs_vus_uncertainty
    plus_five_sigma = (
        compare_only_abs_vus + 5.0 * compare_only_abs_vus_uncertainty
    )
    minimum_sine = math.sin(math.radians(minimum_nonzero))
    direct_match_available = compare_only_angle >= minimum_nonzero - 1e-12
    if direct_match_available:
        raise AssertionError("the direct residual-axis comparison unexpectedly became available")

    return {
        "artifact": "oph_icosahedral_axis_angle_spectrum",
        "schema_version": 1,
        "determinism": {
            "timestamp_omitted": True,
            "construction_field": "Q(sqrt(5))",
            "exact_spectrum_coordinate": "cosine_squared",
        },
        "construction": {
            "oriented_object_counts": {
                "vertices": 12,
                "triangular_faces": 20,
                "edges": 30,
            },
            "unoriented_axis_counts": {
                "fivefold_vertex_axes": len(axes["fivefold"]),
                "threefold_face_axes": len(axes["threefold"]),
                "twofold_edge_axes": len(axes["twofold"]),
                "total": sum(len(family) for family in axes.values()),
            },
            "antipodes_collapsed": True,
        },
        "pairwise_acute_angle_spectra": spectra,
        "method_self_test": {
            "identity": "angle = arctan(1/phi)",
            "expected_family_pair": "fivefold_x_twofold",
            "cosine_squared_exact": _q5_record(golden_cosine_squared),
            "expected_angle_degrees": round(golden_angle, 12),
            "present": True,
            "hard_assertion_passed": True,
        },
        "minimum_nonzero_axis_angle": {
            "angle_degrees": round(minimum_nonzero, 12),
            "family_pair": "threefold_x_twofold",
            "sine": round(minimum_sine, 12),
        },
        "compare_only_cabibbo_readback": {
            "input": {
                "name": "abs_Vus",
                "value": compare_only_abs_vus,
                "display_value": "0.2250",
                "standard_uncertainty": compare_only_abs_vus_uncertainty,
                "published_notation": compare_only_vus[
                    "published_notation"
                ],
                "determination": compare_only_vus["determination"],
                "role": "compare_only_not_used_to_construct_or_select_axes",
                "source": compare_only_vus["source"],
                "fixture_path": compare_only_vus["fixture_path"],
                "fixture_sha256": compare_only_vus["fixture_sha256"],
            },
            "angle_degrees_from_arcsin": round(compare_only_angle, 12),
            "one_standard_uncertainty_interval_abs_Vus": [
                round(one_sigma_lower, 12),
                round(one_sigma_upper, 12),
            ],
            "angle_degrees_at_one_standard_uncertainty_bounds": [
                round(math.degrees(math.asin(one_sigma_lower)), 12),
                round(math.degrees(math.asin(one_sigma_upper)), 12),
            ],
            "plus_five_standard_uncertainties_abs_Vus": round(
                plus_five_sigma, 12
            ),
            "minimum_axis_sine_minus_plus_five_standard_uncertainties": round(
                minimum_sine - plus_five_sigma, 12
            ),
            "direct_match_available_at_plus_five_standard_uncertainties": (
                plus_five_sigma >= minimum_sine - 1e-12
            ),
            "minimum_axis_angle_margin_degrees": round(
                minimum_nonzero - compare_only_angle, 12
            ),
            "minimum_axis_sine_minus_abs_Vus": round(
                minimum_sine - compare_only_abs_vus, 12
            ),
            "minimum_axis_sine_relative_excess": round(
                minimum_sine / compare_only_abs_vus - 1.0, 12
            ),
            "residual_mismatch_cabibbo_available": False,
        },
        "conclusion": {
            "status": "direct_real_3d_residual_axis_identification_excluded",
            "excluded_claim": (
                "The Cabibbo angle is directly equal to an acute angle between "
                "two real three-dimensional icosahedral residual axes."
            ),
            "scope_is_narrow": True,
            "not_excluded": [
                "all_A5_flavor_models",
                "spinorial_or_other_representations",
                "higher_order_symmetry_breaking",
                "arbitrary_overlap_geometry",
            ],
            "literature_claims_used_to_construct_axis_spectrum": False,
            "external_comparison_coordinate_used": True,
        },
    }


def _write_json_lf(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build the exact icosahedral residual-axis angle spectrum."
    )
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--vus-fixture", default=str(VUS_FIXTURE))
    args = parser.parse_args()
    output = Path(args.output)
    _write_json_lf(output, build_artifact(Path(args.vus_fixture)))
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
