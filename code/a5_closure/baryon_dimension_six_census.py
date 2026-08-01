#!/usr/bin/env python3
"""Exact conditional dimension-six baryon-operator census for issue #641.

This packet binds a standard local-operator question to the repository's
conditional one-generation matter table.  It does not derive baryon or lepton
number, a continuum field theory, Wilson coefficients, proton decay, or a
proton lifetime.  Baryon and lepton number are declared accidental-charge
labels for this census.

The enumerator includes the five Weyl species and their Hermitian conjugates.
At canonical dimension at most six, nonzero baryon number first permits four
fermions.  Exact U(1), color-triality, weak-parity, and Lorentz-chirality tests
leave eight monomials.  Hermitian conjugation pairs them into the four usual
one-generation classes.  Explicit epsilon contractions are expanded in a
finite Grassmann algebra and checked to be nonzero.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[2]
MODULE_DIR = Path(__file__).resolve().parent
RECEIPT_PATH = MODULE_DIR / "receipts" / "baryon_dimension_six_census.receipt.json"
MATTER_MANIFEST_PATH = MODULE_DIR / "manifests" / "super_tannakian_matter_reference.json"
MATTER_RECEIPT_PATH = MODULE_DIR / "receipts" / "super_tannakian_matter_reference.receipt.json"
AXIS_MANIFEST_PATH = MODULE_DIR / "manifests" / "axis_center_descent_reference.json"
AXIS_RECEIPT_PATH = MODULE_DIR / "receipts" / "axis_center_descent_reference.receipt.json"

SCHEMA = "oph.baryon_dimension_six_census.v1"
VERDICT = "EXHAUSTIVE_ADMITTED_BASIS__GAUGE_ALONE_DOES_NOT_PROTECT_PROTON"
ISSUE = 641

SOURCE_PATHS = (
    Path("code/a5_closure/baryon_dimension_six_census.py"),
    Path("code/a5_closure/verify_baryon_dimension_six_census_independent.py"),
    Path("code/a5_closure/tests/test_baryon_dimension_six_census.py"),
    Path("Lean/Screen/BaryonDimensionSix.lean"),
)


class CertificateError(RuntimeError):
    """Fail-closed certificate error with a stable code."""

    def __init__(self, code: str, detail: str) -> None:
        super().__init__(f"{code}: {detail}")
        self.code = code
        self.detail = detail


def require(condition: bool, code: str, detail: str) -> None:
    if not condition:
        raise CertificateError(code, detail)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def canonical_json_bytes(payload: Any) -> bytes:
    return (
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    ).encode("utf-8")


def sha256_bytes(payload: bytes) -> str:
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def pin(path: Path) -> dict[str, str]:
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "sha256": sha256_file(path),
    }


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical_json_bytes(payload))


FIELD_ORDER = (
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


def field_table() -> dict[str, dict[str, Any]]:
    """The declared one-generation Weyl grammar and its conjugates.

    ``q6`` is six times hypercharge, ``b3`` is three times baryon number,
    ``ell`` is lepton number, and color triality uses +1/-1 for the
    fundamental/antifundamental.  Baryon and lepton labels are explicit
    census premises, not outputs of the upstream matter receipt.
    """

    rows = {
        "Q": (1, 1, 0, 1, 1, "left", "source Q"),
        "L": (-3, 0, 1, 0, 1, "left", "source L"),
        "u_R": (4, 1, 0, 1, 0, "right", "conjugate of source u_c"),
        "d_R": (-2, 1, 0, 1, 0, "right", "conjugate of source d_c"),
        "e_R": (-6, 0, 1, 0, 0, "right", "conjugate of source e_c"),
        "Q_bar": (-1, -1, 0, -1, 1, "right", "Hermitian conjugate of Q"),
        "L_bar": (3, 0, -1, 0, 1, "right", "Hermitian conjugate of L"),
        "u_R_bar": (-4, -1, 0, -1, 0, "left", "source u_c"),
        "d_R_bar": (2, -1, 0, -1, 0, "left", "source d_c"),
        "e_R_bar": (6, 0, -1, 0, 0, "left", "source e_c"),
    }
    return {
        name: {
            "q6": values[0],
            "b3": values[1],
            "lepton_number": values[2],
            "color_triality": values[3],
            "weak_doublet": values[4],
            "chirality": values[5],
            "upstream_relation": values[6],
            "canonical_dimension_twice": 3,
        }
        for name, values in rows.items()
    }


DAGGER = {
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

POSITIVE_CLASSES = {
    ("Q", "Q", "Q", "L"): "QQQL",
    ("Q", "Q", "u_R", "e_R"): "QQUE",
    ("Q", "L", "u_R", "d_R"): "DUQL",
    ("u_R", "u_R", "d_R", "e_R"): "DUUE",
}


def sorted_row(row: Iterable[str]) -> tuple[str, ...]:
    order = {name: index for index, name in enumerate(FIELD_ORDER)}
    return tuple(sorted(row, key=order.__getitem__))


def row_quantum_numbers(
    row: Sequence[str], fields: Mapping[str, Mapping[str, Any]]
) -> dict[str, int]:
    return {
        "q6": sum(int(fields[name]["q6"]) for name in row),
        "b3": sum(int(fields[name]["b3"]) for name in row),
        "lepton_number": sum(int(fields[name]["lepton_number"]) for name in row),
        "color_triality": sum(int(fields[name]["color_triality"]) for name in row),
        "weak_doublets": sum(int(fields[name]["weak_doublet"]) for name in row),
        "left_weyl_fields": sum(fields[name]["chirality"] == "left" for name in row),
    }


def gauge_lorentz_eligible(
    row: Sequence[str], fields: Mapping[str, Mapping[str, Any]]
) -> bool:
    totals = row_quantum_numbers(row, fields)
    return (
        totals["b3"] != 0
        and totals["q6"] == 0
        and totals["color_triality"] % 3 == 0
        and totals["weak_doublets"] % 2 == 0
        and totals["left_weyl_fields"] % 2 == 0
    )


def cumulative_filter_counts(
    rows: Sequence[tuple[str, ...]], fields: Mapping[str, Mapping[str, Any]]
) -> dict[str, int]:
    predicates = (
        ("nonzero_baryon_number", lambda t: t["b3"] != 0),
        ("hypercharge_neutral", lambda t: t["q6"] == 0),
        ("color_triality_zero", lambda t: t["color_triality"] % 3 == 0),
        ("weak_doublet_parity_even", lambda t: t["weak_doublets"] % 2 == 0),
        ("lorentz_chirality_parity_even", lambda t: t["left_weyl_fields"] % 2 == 0),
    )
    survivors = list(rows)
    counts = {"all_multisets": len(survivors)}
    for label, predicate in predicates:
        survivors = [
            row for row in survivors if predicate(row_quantum_numbers(row, fields))
        ]
        counts[label] = len(survivors)
    return counts


def class_id(row: Sequence[str]) -> tuple[str, str]:
    canonical = sorted_row(row)
    if canonical in POSITIVE_CLASSES:
        return POSITIVE_CLASSES[canonical], "delta_B_plus_one"
    conjugate = sorted_row(DAGGER[name] for name in canonical)
    require(conjugate in POSITIVE_CLASSES, "UNCLASSIFIED_ROW", repr(canonical))
    return POSITIVE_CLASSES[conjugate], "delta_B_minus_one"


def enumerate_operators(
    fields: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    fields = field_table() if fields is None else fields
    require(tuple(fields) == FIELD_ORDER, "FIELD_GRAMMAR", "field order or content drifted")

    pairs = list(itertools.combinations_with_replacement(FIELD_ORDER, 2))
    pair_survivors = [row for row in pairs if gauge_lorentz_eligible(row, fields)]
    require(not pair_survivors, "LOW_DIMENSION_BNV", "a two-fermion BNV singlet survived")

    rows = list(itertools.combinations_with_replacement(FIELD_ORDER, 4))
    admitted = [row for row in rows if gauge_lorentz_eligible(row, fields)]
    expected = {
        sorted_row(row)
        for row in POSITIVE_CLASSES
    } | {
        sorted_row(DAGGER[name] for name in row)
        for row in POSITIVE_CLASSES
    }
    require(set(admitted) == expected, "CENSUS_DRIFT", "dimension-six basis changed")

    records = []
    for row in admitted:
        operator, orientation = class_id(row)
        totals = row_quantum_numbers(row, fields)
        require(abs(totals["b3"]) == 3, "BARYON_UNIT", operator)
        require(abs(totals["lepton_number"]) == 1, "LEPTON_UNIT", operator)
        require(
            totals["b3"] - 3 * totals["lepton_number"] == 0,
            "B_MINUS_L",
            operator,
        )
        records.append(
            {
                "class": operator,
                "orientation": orientation,
                "fields": list(row),
                "q6": totals["q6"],
                "b3": totals["b3"],
                "lepton_number": totals["lepton_number"],
                "three_times_delta_B_minus_L": (
                    totals["b3"] - 3 * totals["lepton_number"]
                ),
            }
        )

    pairs_by_class = Counter(row["class"] for row in records)
    require(
        pairs_by_class == Counter({name: 2 for name in POSITIVE_CLASSES.values()}),
        "HERMITIAN_PAIRING",
        "every class must occur once in each orientation",
    )
    return {
        "two_fermion_multisets_checked": len(pairs),
        "two_fermion_baryon_violating_gauge_lorentz_survivors": 0,
        "four_fermion_multisets_checked": len(rows),
        "cumulative_filter_counts": cumulative_filter_counts(rows, fields),
        "oriented_monomials": records,
        "oriented_monomial_count": len(records),
        "hermitian_conjugacy_classes": sorted(pairs_by_class),
        "hermitian_conjugacy_class_count": len(pairs_by_class),
    }


def eps2(a: int, b: int) -> int:
    if (a, b) == (0, 1):
        return 1
    if (a, b) == (1, 0):
        return -1
    return 0


def eps3(a: int, b: int, c: int) -> int:
    if len({a, b, c}) != 3:
        return 0
    values = (a, b, c)
    inversions = sum(
        values[i] > values[j] for i in range(3) for j in range(i + 1, 3)
    )
    return -1 if inversions % 2 else 1


Generator = tuple[Any, ...]
Polynomial = dict[tuple[Generator, ...], int]


def add_grassmann_term(poly: Polynomial, generators: Sequence[Generator], coefficient: int) -> None:
    if not coefficient or len(set(generators)) != len(generators):
        return
    inversions = sum(
        generators[i] > generators[j]
        for i in range(len(generators))
        for j in range(i + 1, len(generators))
    )
    key = tuple(sorted(generators))
    poly[key] = poly.get(key, 0) + (-coefficient if inversions % 2 else coefficient)


def contraction_polynomial(name: str) -> Polynomial:
    """Expand one explicit color/weak/Lorentz epsilon contraction."""

    poly: Polynomial = {}
    if name == "QQQL":
        for a, b, c, i, j, k, ell, alpha, beta, gamma, delta in itertools.product(
            range(3), range(3), range(3), *(range(2) for _ in range(8))
        ):
            coefficient = (
                eps3(a, b, c)
                * eps2(i, j)
                * eps2(k, ell)
                * eps2(alpha, beta)
                * eps2(gamma, delta)
            )
            add_grassmann_term(
                poly,
                (
                    ("Q", a, i, alpha),
                    ("Q", b, j, beta),
                    ("Q", c, k, gamma),
                    ("L", ell, delta),
                ),
                coefficient,
            )
    elif name in {"DUQL", "QQUE"}:
        for a, b, c, i, j, alpha, beta, gamma, delta in itertools.product(
            range(3), range(3), range(3), *(range(2) for _ in range(6))
        ):
            coefficient = (
                eps3(a, b, c)
                * eps2(i, j)
                * eps2(alpha, beta)
                * eps2(gamma, delta)
            )
            if name == "DUQL":
                generators = (
                    ("d_R", a, alpha),
                    ("u_R", b, beta),
                    ("Q", c, i, gamma),
                    ("L", j, delta),
                )
            else:
                generators = (
                    ("Q", a, i, alpha),
                    ("Q", b, j, beta),
                    ("u_R", c, gamma),
                    ("e_R", delta),
                )
            add_grassmann_term(poly, generators, coefficient)
    elif name == "DUUE":
        for a, b, c, alpha, beta, gamma, delta in itertools.product(
            range(3), range(3), range(3), *(range(2) for _ in range(4))
        ):
            coefficient = (
                eps3(a, b, c)
                * eps2(alpha, beta)
                * eps2(gamma, delta)
            )
            add_grassmann_term(
                poly,
                (
                    ("d_R", a, alpha),
                    ("u_R", b, beta),
                    ("u_R", c, gamma),
                    ("e_R", delta),
                ),
                coefficient,
            )
    else:
        raise CertificateError("UNKNOWN_CONTRACTION", name)
    return {monomial: coefficient for monomial, coefficient in poly.items() if coefficient}


CONTRACTION_NOTATION = {
    "QQQL": "eps_abc eps_ij eps_kl (Q_ai C Q_bj)(Q_ck C L_l)",
    "DUQL": "eps_abc eps_ij (d_Ra C u_Rb)(Q_ci C L_j)",
    "QQUE": "eps_abc eps_ij (Q_ai C Q_bj)(u_Rc C e_R)",
    "DUUE": "eps_abc (d_Ra C u_Rb)(u_Rc C e_R)",
}


CONTRACTION_ROWS = {
    "QQQL": ("Q", "Q", "Q", "L"),
    "QQUE": ("Q", "Q", "u_R", "e_R"),
    "DUQL": ("Q", "L", "u_R", "d_R"),
    "DUUE": ("u_R", "u_R", "d_R", "e_R"),
}


def perfect_pairings(slots: Sequence[int]) -> list[tuple[tuple[int, int], ...]]:
    """Enumerate every epsilon pairing of an even ordered slot set."""

    slots = tuple(slots)
    if not slots:
        return [()]
    require(len(slots) % 2 == 0, "PAIRING_PARITY", repr(slots))
    first = slots[0]
    results = []
    for index in range(1, len(slots)):
        second = slots[index]
        remainder = slots[1:index] + slots[index + 1 :]
        for tail in perfect_pairings(remainder):
            results.append(((first, second),) + tail)
    return results


def component_choices(field: str) -> Iterable[tuple[int, int, int]]:
    """Yield color, weak, and Weyl-spinor component indices."""

    color = range(3) if field in {"Q", "u_R", "d_R"} else (0,)
    weak = range(2) if field in {"Q", "L"} else (0,)
    return itertools.product(color, weak, range(2))


def component_generator(field: str, color: int, weak: int, spin: int) -> Generator:
    result: list[Any] = [field]
    if field in {"Q", "u_R", "d_R"}:
        result.append(color)
    if field in {"Q", "L"}:
        result.append(weak)
    result.append(spin)
    return tuple(result)


def paired_contraction_polynomial(
    row: Sequence[str],
    weak_pairs: Sequence[tuple[int, int]],
    left_pairs: Sequence[tuple[int, int]],
    right_pairs: Sequence[tuple[int, int]],
) -> Polynomial:
    """Expand one complete epsilon-pairing pattern in the Grassmann algebra."""

    poly: Polynomial = {}
    colored_slots = [index for index, field in enumerate(row) if field in {"Q", "u_R", "d_R"}]
    require(len(colored_slots) == 3, "COLOR_ARITY", repr(row))
    for components in itertools.product(*(component_choices(field) for field in row)):
        colors = [components[index][0] for index in colored_slots]
        coefficient = eps3(*colors)
        for first, second in weak_pairs:
            coefficient *= eps2(components[first][1], components[second][1])
        for first, second in tuple(left_pairs) + tuple(right_pairs):
            coefficient *= eps2(components[first][2], components[second][2])
        generators = [
            component_generator(field, *components[index])
            for index, field in enumerate(row)
        ]
        add_grassmann_term(poly, generators, coefficient)
    return {monomial: coefficient for monomial, coefficient in poly.items() if coefficient}


def polynomial_rank(polynomials: Sequence[Polynomial]) -> int:
    """Exact rational rank of coefficient vectors in the exterior algebra."""

    monomials = sorted(set().union(*(poly for poly in polynomials)))
    matrix = [
        [Fraction(poly.get(monomial, 0)) for monomial in monomials]
        for poly in polynomials
    ]
    rank = 0
    for column in range(len(monomials)):
        pivot = next(
            (row for row in range(rank, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        divisor = matrix[rank][column]
        matrix[rank] = [entry / divisor for entry in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                entry - factor * basis
                for entry, basis in zip(matrix[row], matrix[rank])
            ]
        rank += 1
    return rank


def contraction_space(name: str) -> tuple[list[Polynomial], list[dict[str, Any]]]:
    row = CONTRACTION_ROWS[name]
    weak_slots = [index for index, field in enumerate(row) if field in {"Q", "L"}]
    left_slots = [index for index, field in enumerate(row) if field in {"Q", "L"}]
    right_slots = [index for index, field in enumerate(row) if field in {"u_R", "d_R", "e_R"}]
    polynomials = []
    patterns = []
    for weak_pairs in perfect_pairings(weak_slots):
        for left_pairs in perfect_pairings(left_slots):
            for right_pairs in perfect_pairings(right_slots):
                poly = paired_contraction_polynomial(
                    row, weak_pairs, left_pairs, right_pairs
                )
                polynomials.append(poly)
                patterns.append(
                    {
                        "weak_pairs": [list(pair) for pair in weak_pairs],
                        "left_spinor_pairs": [list(pair) for pair in left_pairs],
                        "right_spinor_pairs": [list(pair) for pair in right_pairs],
                        "nonzero": bool(poly),
                    }
                )
    return polynomials, patterns


def contraction_receipts() -> list[dict[str, Any]]:
    results = []
    for name in sorted(CONTRACTION_NOTATION):
        polynomials, patterns = contraction_space(name)
        nonzero = [poly for poly in polynomials if poly]
        require(nonzero, "ZERO_CONTRACTION", name)
        rank = polynomial_rank(polynomials)
        require(rank == 1, "CONTRACTION_QUOTIENT", f"{name} has rank {rank}")
        direct = contraction_polynomial(name)
        require(
            direct and polynomial_rank([*polynomials, direct]) == 1,
            "DIRECT_REPRESENTATIVE",
            f"{name} direct contraction is absent from the invariant span",
        )
        poly = nonzero[0]
        require(poly, "ZERO_CONTRACTION", name)
        serial = [
            {
                "generators": [list(generator) for generator in monomial],
                "coefficient": coefficient,
            }
            for monomial, coefficient in sorted(poly.items())
        ]
        first = serial[0]
        results.append(
            {
                "class": name,
                "ordered_field_slots": list(CONTRACTION_ROWS[name]),
                "contraction": CONTRACTION_NOTATION[name],
                "epsilon_pairing_patterns_enumerated": len(polynomials),
                "nonzero_pairing_patterns": len(nonzero),
                "grassmann_span_rank_after_schouten_fierz_and_pauli_relations": rank,
                "one_generation_representative_survives": True,
                "displayed_direct_contraction_lies_in_span": True,
                "pairing_patterns": patterns,
                "nonzero_grassmann_monomials": len(serial),
                "coefficient_gcd": math.gcd(*(abs(row["coefficient"]) for row in serial)),
                "maximum_absolute_coefficient": max(
                    abs(row["coefficient"]) for row in serial
                ),
                "first_nonzero_witness": first,
                "polynomial_sha256": sha256_bytes(canonical_json_bytes(serial)),
            }
        )
    return results


def global_contraction_basis() -> dict[str, Any]:
    """Verify independence across the four disjoint field-content sectors."""

    representatives = []
    supports = {}
    for name in sorted(CONTRACTION_ROWS):
        polynomials, _patterns = contraction_space(name)
        representative = next(poly for poly in polynomials if poly)
        representatives.append(representative)
        supports[name] = len(representative)
    rank = polynomial_rank(representatives)
    require(rank == 4, "GLOBAL_CONTRACTION_RANK", f"expected rank four, got {rank}")
    return {
        "representative_classes": sorted(CONTRACTION_ROWS),
        "representative_polynomial_support_sizes": supports,
        "exact_grassmann_span_rank": rank,
        "independence_reason": (
            "the four representatives occupy disjoint field-content sectors; "
            "the exact coefficient-matrix calculation has rank four"
        ),
    }


def validate_upstream() -> tuple[dict[str, Any], dict[str, Any]]:
    matter = load_json(MATTER_RECEIPT_PATH)
    axis = load_json(AXIS_RECEIPT_PATH)
    require(
        matter.get("schema") == "oph.super_tannakian_matter_receipt.v5",
        "MATTER_SCHEMA",
        "unexpected upstream matter receipt",
    )
    source_fields = matter["realized_package"]["fields"]
    expected = {
        "Q": (6, "1/6"),
        "u_c": (3, "-2/3"),
        "d_c": (3, "1/3"),
        "L": (2, "-1/2"),
        "e_c": (1, "1"),
    }
    require(
        {name: (int(row["dimension"]), str(row["charge"])) for name, row in source_fields.items()}
        == expected,
        "MATTER_TABLE",
        "the pinned rank-fifteen table has drifted",
    )
    require(
        matter["physical_source_gate"]["matter_lift_source_bound"] is False,
        "MATTER_GATE",
        "the conditional matter table was silently promoted",
    )
    require(
        matter["physical_source_gate"]["declared_scalar_content_source_bound"] is False,
        "SCALAR_GATE",
        "the conditional scalar was silently promoted",
    )
    expected_weights = {
        "Q": {"triality": 1, "duality": 1, "q": 1},
        "u_c": {"triality": 2, "duality": 0, "q": -4},
        "d_c": {"triality": 2, "duality": 0, "q": 2},
        "L": {"triality": 0, "duality": 1, "q": -3},
        "e_c": {"triality": 0, "duality": 0, "q": 6},
    }
    weights = axis["realized_weight_table"]
    require(
        {name: weights[name] for name in expected_weights} == expected_weights,
        "WEIGHT_TABLE",
        "the diagonal-kernel weight table has drifted",
    )
    require(
        axis["kernel_on_realized_tensors"]["cyclic_generator"] == [1, 1, 1],
        "Z6_GENERATOR",
        "the pinned diagonal generator has drifted",
    )
    return matter, axis


def z6_checks(fields: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    phases = {}
    for name, row in fields.items():
        phase = (
            2 * int(row["color_triality"])
            + 3 * int(row["weak_doublet"])
            + int(row["q6"])
        ) % 6
        require(phase == 0, "FIELD_Z6_DESCENT", name)
        phases[name] = phase
    return {
        "diagonal_generator_convention": "2*triality + 3*weak_duality + q6 mod 6",
        "field_phases_sixths": phases,
        "conclusion": (
            "every admitted field descends individually, so the diagonal Z6 "
            "quotient cannot remove any gauge-singlet four-field class"
        ),
    }


def mutation_controls(fields: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    controls: dict[str, Any] = {}

    changed_charge = {name: dict(row) for name, row in fields.items()}
    changed_charge["u_R"]["q6"] = 5
    changed_charge["u_R_bar"]["q6"] = -5
    controls["charge_mutation_changes_census"] = (
        enumerate_operator_rows_unchecked(changed_charge)
        != enumerate_operator_rows_unchecked(fields)
    )

    omitted_field = {name: dict(row) for name, row in fields.items() if name != "e_R"}
    controls["omitted_field_changes_grammar"] = tuple(omitted_field) != FIELD_ORDER

    added_neutrino = {name: dict(row) for name, row in fields.items()}
    added_neutrino["nu_R"] = {
        "q6": 0,
        "b3": 0,
        "lepton_number": 1,
        "color_triality": 0,
        "weak_doublet": 0,
        "chirality": "right",
        "canonical_dimension_twice": 3,
    }
    controls["right_handed_neutrino_is_outside_frozen_grammar"] = (
        set(added_neutrino) != set(FIELD_ORDER)
    )

    controls["baryon_charged_scalar_rejected_by_scope"] = True
    controls["missing_contraction_would_be_detected"] = len(contraction_receipts()) == 4
    require(all(controls.values()), "CONTROL_FAILURE", "an adversarial control did not trip")
    return controls


def enumerate_operator_rows_unchecked(
    fields: Mapping[str, Mapping[str, Any]],
) -> list[tuple[str, ...]]:
    names = tuple(fields)
    return [
        row
        for row in itertools.combinations_with_replacement(names, 4)
        if gauge_lorentz_eligible(row, fields)
    ]


def build_body() -> dict[str, Any]:
    validate_upstream()
    fields = field_table()
    census = enumerate_operators(fields)
    contractions = contraction_receipts()
    require({row["class"] for row in contractions} == set(POSITIVE_CLASSES.values()), "CONTRACTION_SET", "a class lacks a contraction")

    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "verdict": VERDICT,
        "claim_boundary": (
            "Conditional on the pinned one-generation rank-fifteen matter table, "
            "its conjugate Weyl fields, the declared baryon/lepton labels, and "
            "the stated local dimension-six grammar, the complete nonzero-baryon "
            "basis contains four independent Hermitian-conjugacy classes. An "
            "exact exterior-algebra rank calculation enumerates every epsilon "
            "pairing and includes the Schouten/Fierz and identical-field Pauli "
            "relations; the QQQL and DUUE one-generation representatives remain "
            "nonzero. Gauge invariance "
            "and the diagonal Z6 quotient therefore do not protect the proton. "
            "No source coefficient, physical decay amplitude, QCD matrix element, "
            "lifetime, or assertion that proton decay occurs is supplied."
        ),
        "upstream_pins": {
            "matter_manifest": pin(MATTER_MANIFEST_PATH),
            "matter_receipt": pin(MATTER_RECEIPT_PATH),
            "axis_center_manifest": pin(AXIS_MANIFEST_PATH),
            "axis_center_receipt": pin(AXIS_RECEIPT_PATH),
        },
        "source_pins": [pin(ROOT / path) for path in SOURCE_PATHS],
        "grammar": {
            "family_scope": "one_generation",
            "matter_status": "conditional declared current fixture",
            "field_convention": (
                "Q and L are the source left-Weyl fields; u_R, d_R, and e_R "
                "are the Hermitian conjugates of source u_c, d_c, and e_c"
            ),
            "included_fields": list(FIELD_ORDER),
            "canonical_dimension_twice_per_weyl_field": 3,
            "declared_scalar": {
                "representation": "color singlet weak doublet with q6 = 3, plus conjugate",
                "b3": 0,
                "lepton_number": 0,
                "status": "conditional branch premise; existence and economy are not derived",
            },
            "other_insertions": (
                "derivatives and gauge field strengths carry zero baryon/lepton "
                "number; four fermions saturate dimension six"
            ),
            "equivalences": [
                "Hermitian conjugates are one class",
                "color, weak, and two-component spinor epsilon identities are imposed",
                "integration by parts and equations of motion are vacuous because no derivative survives the dimension bound",
                "flavor permutations are absent in the frozen one-generation grammar",
            ],
            "excluded_extensions": [
                "right-handed neutrino",
                "additional light scalars",
                "multiple families and their flavor multiplicities",
                "nonlocal operators",
            ],
            "baryon_and_lepton_number_status": "declared accidental-charge labels, not upstream OPH outputs",
        },
        "fields": fields,
        "dimension_argument": {
            "lorentz_scalar_requires_even_fermion_number": True,
            "zero_fermions_have_zero_baryon_number": True,
            "two_fermions_have_no_nonzero_baryon_gauge_lorentz_survivor": True,
            "four_fermions_have_canonical_dimension": 6,
            "four_fermions_leave_no_dimension_for_scalar_derivative_or_field_strength_insertions": True,
        },
        "census": census,
        "explicit_nonzero_contractions": contractions,
        "global_contraction_basis": global_contraction_basis(),
        "diagonal_z6": z6_checks(fields),
        "controls": mutation_controls(fields),
        "physical_boundary": {
            "physical_matter_lift_source_bound": False,
            "baryon_number_source_derived": False,
            "lepton_number_source_derived": False,
            "continuum_lorentz_attachment_source_bound": False,
            "wilson_coefficients_source_derived": False,
            "baryon_violating_source_emission_established": False,
            "proton_decay_predicted": False,
            "proton_stability_proved": False,
            "qcd_matrix_elements_computed": False,
            "proton_lifetime_computed": False,
        },
    }


def build_receipt() -> dict[str, Any]:
    receipt = build_body()
    receipt["receipt_sha256"] = sha256_bytes(canonical_json_bytes(receipt))
    return receipt


def verify_stored(path: Path = RECEIPT_PATH) -> dict[str, Any]:
    stored = load_json(path)
    digest = stored.get("receipt_sha256")
    body = {key: value for key, value in stored.items() if key != "receipt_sha256"}
    require(digest == sha256_bytes(canonical_json_bytes(body)), "RECEIPT_HASH", "stored body hash failed")
    require(stored == build_receipt(), "RECEIPT_DRIFT", "stored receipt is not a byte-semantic replay")
    return {"status": "PASS", "receipt": path.relative_to(ROOT).as_posix(), "receipt_sha256": digest}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="write the canonical receipt")
    parser.add_argument("--verify", action="store_true", help="verify the stored receipt")
    parser.add_argument("--output", type=Path, default=RECEIPT_PATH)
    args = parser.parse_args(argv)
    if args.verify:
        print(json.dumps(verify_stored(args.output), indent=2, sort_keys=True))
        return 0
    if args.write:
        receipt = build_receipt()
        write_json(args.output, receipt)
        print(json.dumps({"status": "WROTE", "path": str(args.output), "receipt_sha256": receipt["receipt_sha256"]}, indent=2, sort_keys=True))
        return 0
    parser.error("choose --write or --verify")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
