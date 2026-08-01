#!/usr/bin/env python3
"""Independent standard-library verifier for the issue #641 receipt.

This module deliberately imports neither the producer nor repository receipt
helpers.  It reconstructs the finite field census and the full epsilon-pairing
spaces from the stored grammar, checks their exact Grassmann ranks, validates
all raw source and parent pins, and refuses any physical promotion.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RECEIPT = ROOT / "code/a5_closure/receipts/baryon_dimension_six_census.receipt.json"
EXPECTED_SCHEMA = "oph.baryon_dimension_six_census.v1"
EXPECTED_VERDICT = "EXHAUSTIVE_ADMITTED_BASIS__GAUGE_ALONE_DOES_NOT_PROTECT_PROTON"

ORDER = (
    "Q",
    "L",
    "u_R",
    "d_R",
    "e_R",
    "Q_bar",
    "L_bar",
    "u_R_bar",
    "d_R_bar",
    "e_R_bar",
)

EXPECTED_ROWS = {
    ("Q", "Q", "Q", "L"): "QQQL",
    ("Q", "Q", "u_R", "e_R"): "QQUE",
    ("Q", "L", "u_R", "d_R"): "DUQL",
    ("u_R", "u_R", "d_R", "e_R"): "DUUE",
}

CONJUGATE = {
    "Q": "Q_bar",
    "L": "L_bar",
    "u_R": "u_R_bar",
    "d_R": "d_R_bar",
    "e_R": "e_R_bar",
    "Q_bar": "Q",
    "L_bar": "L",
    "u_R_bar": "u_R",
    "d_R_bar": "d_R",
    "e_R_bar": "e_R",
}


class IndependentVerificationError(RuntimeError):
    pass


def ensure(condition: bool, message: str) -> None:
    if not condition:
        raise IndependentVerificationError(message)


def canonical_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    ).encode("utf-8")


def digest_bytes(value: bytes) -> str:
    return "sha256:" + hashlib.sha256(value).hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def normalized(row: Iterable[str]) -> tuple[str, ...]:
    positions = {name: index for index, name in enumerate(ORDER)}
    return tuple(sorted(row, key=positions.__getitem__))


def totals(row: Sequence[str], table: Mapping[str, Mapping[str, Any]]) -> tuple[int, ...]:
    return (
        sum(int(table[name]["q6"]) for name in row),
        sum(int(table[name]["b3"]) for name in row),
        sum(int(table[name]["lepton_number"]) for name in row),
        sum(int(table[name]["color_triality"]) for name in row),
        sum(int(table[name]["weak_doublet"]) for name in row),
        sum(table[name]["chirality"] == "left" for name in row),
    )


def admissible(row: Sequence[str], table: Mapping[str, Mapping[str, Any]]) -> bool:
    q6, b3, _ell, triality, weak, left = totals(row, table)
    return b3 != 0 and q6 == 0 and triality % 3 == 0 and weak % 2 == 0 and left % 2 == 0


def wedge_sign(generators: Sequence[tuple[Any, ...]]) -> tuple[tuple[tuple[Any, ...], ...] | None, int]:
    if len(set(generators)) != len(generators):
        return None, 0
    crossings = sum(
        generators[i] > generators[j]
        for i in range(len(generators))
        for j in range(i + 1, len(generators))
    )
    return tuple(sorted(generators)), (-1 if crossings % 2 else 1)


def e2(i: int, j: int) -> int:
    return {(0, 1): 1, (1, 0): -1}.get((i, j), 0)


def e3(i: int, j: int, k: int) -> int:
    if len({i, j, k}) != 3:
        return 0
    values = (i, j, k)
    crossings = sum(
        values[a] > values[b] for a in range(3) for b in range(a + 1, 3)
    )
    return -1 if crossings % 2 else 1


def matchings(slots: tuple[int, ...]) -> list[tuple[tuple[int, int], ...]]:
    if not slots:
        return [()]
    ensure(len(slots) % 2 == 0, "odd epsilon slot count")
    head = slots[0]
    answer = []
    for offset in range(1, len(slots)):
        partner = slots[offset]
        rest = slots[1:offset] + slots[offset + 1 :]
        answer.extend(((head, partner),) + tail for tail in matchings(rest))
    return answer


def components(field: str) -> Iterable[tuple[int, int, int]]:
    colors = range(3) if field in {"Q", "u_R", "d_R"} else (0,)
    weak = range(2) if field in {"Q", "L"} else (0,)
    return itertools.product(colors, weak, range(2))


def variable(field: str, color: int, weak: int, spin: int) -> tuple[Any, ...]:
    label: list[Any] = [field]
    if field in {"Q", "u_R", "d_R"}:
        label.append(color)
    if field in {"Q", "L"}:
        label.append(weak)
    label.append(spin)
    return tuple(label)


def expand_pattern(
    row: Sequence[str],
    weak_pairs: Sequence[tuple[int, int]],
    left_pairs: Sequence[tuple[int, int]],
    right_pairs: Sequence[tuple[int, int]],
) -> dict[tuple[tuple[Any, ...], ...], int]:
    polynomial: dict[tuple[tuple[Any, ...], ...], int] = {}
    colored = tuple(index for index, field in enumerate(row) if field in {"Q", "u_R", "d_R"})
    ensure(len(colored) == 3, "a BNV row must have three colored slots")
    for assignment in itertools.product(*(components(field) for field in row)):
        coefficient = e3(*(assignment[index][0] for index in colored))
        for a, b in weak_pairs:
            coefficient *= e2(assignment[a][1], assignment[b][1])
        for a, b in tuple(left_pairs) + tuple(right_pairs):
            coefficient *= e2(assignment[a][2], assignment[b][2])
        generators = [
            variable(field, *assignment[index]) for index, field in enumerate(row)
        ]
        monomial, sign = wedge_sign(generators)
        if coefficient and sign:
            assert monomial is not None
            polynomial[monomial] = polynomial.get(monomial, 0) + sign * coefficient
    return {key: value for key, value in polynomial.items() if value}


def exact_rank(polynomials: Sequence[Mapping[Any, int]]) -> int:
    columns = sorted(set().union(*(polynomial for polynomial in polynomials)))
    rows = [
        [Fraction(polynomial.get(column, 0)) for column in columns]
        for polynomial in polynomials
    ]
    rank = 0
    for column in range(len(columns)):
        pivot = next((row for row in range(rank, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][column]
        rows[rank] = [entry / scale for entry in rows[rank]]
        for row in range(len(rows)):
            if row == rank or rows[row][column] == 0:
                continue
            scale = rows[row][column]
            rows[row] = [entry - scale * basis for entry, basis in zip(rows[row], rows[rank])]
        rank += 1
    return rank


def reconstruct_contraction(row: Sequence[str]) -> tuple[list[dict[Any, int]], list[dict[str, Any]]]:
    weak = tuple(index for index, field in enumerate(row) if field in {"Q", "L"})
    left = weak
    right = tuple(index for index, field in enumerate(row) if field in {"u_R", "d_R", "e_R"})
    polynomials = []
    patterns = []
    for weak_pairs in matchings(weak):
        for left_pairs in matchings(left):
            for right_pairs in matchings(right):
                polynomial = expand_pattern(row, weak_pairs, left_pairs, right_pairs)
                polynomials.append(polynomial)
                patterns.append(
                    {
                        "weak_pairs": [list(pair) for pair in weak_pairs],
                        "left_spinor_pairs": [list(pair) for pair in left_pairs],
                        "right_spinor_pairs": [list(pair) for pair in right_pairs],
                        "nonzero": bool(polynomial),
                    }
                )
    return polynomials, patterns


def polynomial_summary(polynomial: Mapping[tuple[tuple[Any, ...], ...], int]) -> tuple[int, int, int, str, dict[str, Any]]:
    serialized = [
        {
            "generators": [list(generator) for generator in monomial],
            "coefficient": coefficient,
        }
        for monomial, coefficient in sorted(polynomial.items())
    ]
    gcd = 0
    for row in serialized:
        value = abs(int(row["coefficient"]))
        while value:
            gcd, value = value, gcd % value
    maximum = max(abs(int(row["coefficient"])) for row in serialized)
    return len(serialized), gcd, maximum, digest_bytes(canonical_bytes(serialized)), serialized[0]


def validate_pins(document: Mapping[str, Any]) -> None:
    source_paths = {
        "code/a5_closure/baryon_dimension_six_census.py",
        "code/a5_closure/verify_baryon_dimension_six_census_independent.py",
        "code/a5_closure/tests/test_baryon_dimension_six_census.py",
        "Lean/Screen/BaryonDimensionSix.lean",
    }
    pins = document["source_pins"]
    ensure({pin["path"] for pin in pins} == source_paths, "source pin inventory mismatch")
    for pin in pins:
        path = ROOT / pin["path"]
        ensure(path.is_file(), f"missing pinned source {pin['path']}")
        ensure(digest_bytes(path.read_bytes()) == pin["sha256"], f"source digest mismatch: {pin['path']}")

    parents = document["upstream_pins"]
    expected_parent_paths = {
        "matter_manifest": "code/a5_closure/manifests/super_tannakian_matter_reference.json",
        "matter_receipt": "code/a5_closure/receipts/super_tannakian_matter_reference.receipt.json",
        "axis_center_manifest": "code/a5_closure/manifests/axis_center_descent_reference.json",
        "axis_center_receipt": "code/a5_closure/receipts/axis_center_descent_reference.receipt.json",
    }
    for label, expected_path in expected_parent_paths.items():
        pin = parents[label]
        ensure(pin["path"] == expected_path, f"wrong {label} path")
        ensure(digest_bytes((ROOT / expected_path).read_bytes()) == pin["sha256"], f"wrong {label} digest")


def validate_upstream_semantics(document: Mapping[str, Any]) -> None:
    matter_path = ROOT / document["upstream_pins"]["matter_receipt"]["path"]
    axis_path = ROOT / document["upstream_pins"]["axis_center_receipt"]["path"]
    matter = read_json(matter_path)
    axis = read_json(axis_path)
    expected_source = {
        "Q": (6, "1/6"),
        "u_c": (3, "-2/3"),
        "d_c": (3, "1/3"),
        "L": (2, "-1/2"),
        "e_c": (1, "1"),
    }
    actual = {
        name: (int(row["dimension"]), str(row["charge"]))
        for name, row in matter["realized_package"]["fields"].items()
    }
    ensure(actual == expected_source, "upstream rank-fifteen table mismatch")
    gate = matter["physical_source_gate"]
    ensure(gate["matter_lift_source_bound"] is False, "matter source gate promoted")
    ensure(gate["declared_scalar_content_source_bound"] is False, "scalar source gate promoted")
    ensure(axis["kernel_on_realized_tensors"]["cyclic_generator"] == [1, 1, 1], "Z6 generator mismatch")


def verify_document(document: Mapping[str, Any], *, verify_files: bool = True) -> dict[str, Any]:
    ensure(document.get("schema") == EXPECTED_SCHEMA, "schema mismatch")
    ensure(document.get("issue") == 641, "issue mismatch")
    ensure(document.get("verdict") == EXPECTED_VERDICT, "verdict mismatch")
    claimed_digest = document.get("receipt_sha256")
    body = {key: value for key, value in document.items() if key != "receipt_sha256"}
    ensure(claimed_digest == digest_bytes(canonical_bytes(body)), "receipt digest mismatch")
    if verify_files:
        validate_pins(document)
        validate_upstream_semantics(document)

    grammar = document["grammar"]
    ensure(grammar["family_scope"] == "one_generation", "family scope widened")
    ensure(grammar["included_fields"] == list(ORDER), "field grammar mismatch")
    scalar = grammar["declared_scalar"]
    ensure(scalar["b3"] == 0 and scalar["lepton_number"] == 0, "scalar carries B or L")
    ensure("right-handed neutrino" in grammar["excluded_extensions"], "right-handed neutrino boundary absent")

    table = document["fields"]
    ensure(set(table) == set(ORDER) and len(table) == len(ORDER), "field table inventory mismatch")
    expected_numbers = {
        "Q": (1, 1, 0, 1, 1, "left"),
        "L": (-3, 0, 1, 0, 1, "left"),
        "u_R": (4, 1, 0, 1, 0, "right"),
        "d_R": (-2, 1, 0, 1, 0, "right"),
        "e_R": (-6, 0, 1, 0, 0, "right"),
        "Q_bar": (-1, -1, 0, -1, 1, "right"),
        "L_bar": (3, 0, -1, 0, 1, "right"),
        "u_R_bar": (-4, -1, 0, -1, 0, "left"),
        "d_R_bar": (2, -1, 0, -1, 0, "left"),
        "e_R_bar": (6, 0, -1, 0, 0, "left"),
    }
    for name, expected in expected_numbers.items():
        row = table[name]
        actual = (
            int(row["q6"]),
            int(row["b3"]),
            int(row["lepton_number"]),
            int(row["color_triality"]),
            int(row["weak_doublet"]),
            row["chirality"],
        )
        ensure(actual == expected, f"field numbers changed for {name}")
        phase = (2 * actual[3] + 3 * actual[4] + actual[0]) % 6
        ensure(phase == 0, f"field does not descend through Z6: {name}")
        ensure(document["diagonal_z6"]["field_phases_sixths"][name] == 0, f"stored Z6 phase changed for {name}")

    pairs = list(itertools.combinations_with_replacement(ORDER, 2))
    ensure(not [row for row in pairs if admissible(row, table)], "two-fermion BNV survivor found")
    four = list(itertools.combinations_with_replacement(ORDER, 4))
    admitted = [row for row in four if admissible(row, table)]
    expected_positive = set(EXPECTED_ROWS)
    expected_all = expected_positive | {
        normalized(CONJUGATE[field] for field in row) for row in expected_positive
    }
    ensure(set(admitted) == expected_all and len(admitted) == 8, "four-fermion census mismatch")

    stored_rows = document["census"]["oriented_monomials"]
    ensure({tuple(row["fields"]) for row in stored_rows} == expected_all, "stored operator rows mismatch")
    ensure(document["census"]["four_fermion_multisets_checked"] == 715, "enumeration domain mismatch")
    ensure(document["census"]["hermitian_conjugacy_class_count"] == 4, "class count mismatch")
    for row in stored_rows:
        q6, b3, ell, _triality, _weak, _left = totals(row["fields"], table)
        ensure(q6 == 0 and abs(b3) == 3 and abs(ell) == 1 and b3 - 3 * ell == 0, "operator charge mismatch")

    stored_contractions = {row["class"]: row for row in document["explicit_nonzero_contractions"]}
    ensure(set(stored_contractions) == set(EXPECTED_ROWS.values()), "contraction class set mismatch")
    representatives = []
    for row, name in EXPECTED_ROWS.items():
        polynomials, patterns = reconstruct_contraction(row)
        nonzero = [polynomial for polynomial in polynomials if polynomial]
        representatives.append(nonzero[0])
        ensure(nonzero, f"{name} vanishes after Pauli antisymmetry")
        ensure(exact_rank(polynomials) == 1, f"{name} invariant contraction rank is not one")
        count, gcd, maximum, digest, witness = polynomial_summary(nonzero[0])
        stored = stored_contractions[name]
        ensure(stored["ordered_field_slots"] == list(row), f"{name} slot order mismatch")
        ensure(stored["epsilon_pairing_patterns_enumerated"] == len(polynomials), f"{name} pattern count mismatch")
        ensure(stored["nonzero_pairing_patterns"] == len(nonzero), f"{name} nonzero pattern count mismatch")
        ensure(stored["grassmann_span_rank_after_schouten_fierz_and_pauli_relations"] == 1, f"{name} rank receipt mismatch")
        ensure(stored["pairing_patterns"] == patterns, f"{name} pairing inventory mismatch")
        ensure(stored["nonzero_grassmann_monomials"] == count, f"{name} monomial count mismatch")
        ensure(stored["coefficient_gcd"] == gcd, f"{name} coefficient gcd mismatch")
        ensure(stored["maximum_absolute_coefficient"] == maximum, f"{name} coefficient maximum mismatch")
        ensure(stored["polynomial_sha256"] == digest, f"{name} polynomial digest mismatch")
        ensure(stored["first_nonzero_witness"] == witness, f"{name} witness mismatch")
        ensure(stored["one_generation_representative_survives"] is True, f"{name} survival flag missing")

    global_basis = document["global_contraction_basis"]
    ensure(exact_rank(representatives) == 4, "four representatives are not independent")
    ensure(global_basis["exact_grassmann_span_rank"] == 4, "stored global contraction rank mismatch")
    ensure(
        set(global_basis["representative_classes"]) == set(EXPECTED_ROWS.values()),
        "stored global class inventory mismatch",
    )

    boundary = document["physical_boundary"]
    ensure(boundary and all(value is False for value in boundary.values()), "physical boundary promoted")
    ensure(all(document["controls"].values()), "stored adversarial control failed")
    return {
        "status": "PASS",
        "verdict": EXPECTED_VERDICT,
        "oriented_monomials": 8,
        "independent_classes": 4,
        "global_contraction_rank": 4,
        "qqql_survives": True,
        "duue_survives": True,
    }


def audit(path: Path = DEFAULT_RECEIPT) -> dict[str, Any]:
    return verify_document(read_json(path), verify_files=True)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path, nargs="?", default=DEFAULT_RECEIPT)
    args = parser.parse_args(argv)
    print(json.dumps(audit(args.receipt), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
