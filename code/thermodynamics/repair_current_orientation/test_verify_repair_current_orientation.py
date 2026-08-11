from __future__ import annotations

import copy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

import pytest


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from verify_repair_current_orientation import (  # noqa: E402
    DEFAULT_PAYLOAD,
    EXPECTED_COUNTS,
    EXPECTED_NORMALIZED_CYCLE,
    EXPECTED_NORMALIZED_PAIR,
    EXPECTED_PAYLOAD_SHA256,
    VerificationError,
    select_cycle,
    select_pair,
    symmetric_control,
    transpose,
    verify_file,
    verify_lean_binding,
    verify_normalized_designations,
    verify_payload,
    verify_symmetric_control,
    verify_transposition_semantics,
)


@pytest.fixture
def payload() -> dict:
    return json.loads(DEFAULT_PAYLOAD.read_text(encoding="utf-8"))


def test_vendored_schema_v3_payload_verifies_exactly() -> None:
    report = verify_file()
    assert report["status"] == "PASS"
    assert report["payload_sha256"] == EXPECTED_PAYLOAD_SHA256
    assert report["analysis_status"] == "post_hoc_diagnostic"
    assert report["eligible_as_validation"] is False
    assert report["normalized_pair_maximizer"]["buckets"] == [2, 3]
    assert report["normalized_cycle_maximizer"]["buckets"] == [3, 4, 5]


def test_vendored_raw_bytes_have_the_pinned_sha() -> None:
    assert hashlib.sha256(DEFAULT_PAYLOAD.read_bytes()).hexdigest() == (
        EXPECTED_PAYLOAD_SHA256
    )


def test_lean_count_literal_is_bound_to_the_vendored_payload() -> None:
    verify_lean_binding()


def test_lean_count_literal_mutation_is_rejected(tmp_path: Path) -> None:
    source = (HERE.parents[2] / "Lean" / "Thermodynamics" / "RepairCurrentOrientation.lean").read_text(
        encoding="utf-8"
    )
    path = tmp_path / "RepairCurrentOrientation.lean"
    path.write_text(source.replace("1343, 0, 46", "1342, 0, 46", 1), encoding="utf-8")
    with pytest.raises(VerificationError, match="differs from the vendored payload"):
        verify_lean_binding(path)


def test_unresigned_raw_byte_mutation_is_rejected(tmp_path: Path) -> None:
    mutant = DEFAULT_PAYLOAD.read_bytes().replace(
        b'"post_hoc_diagnostic"', b'"post_hoc_diagnostic "'
    )
    path = tmp_path / "mutant.json"
    path.write_bytes(mutant)
    with pytest.raises(VerificationError, match="SHA-256"):
        verify_file(path)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    (
        ("analysis_status", "preregistered_validation", "analysis status"),
        ("statistic_preregistered", True, "statistic falsely"),
        ("designation_rule_preregistered", True, "designation rule falsely"),
        ("eligible_as_validation", True, "promoted to validation"),
    ),
)
def test_epistemic_promotion_mutations_fail_closed(
    payload: dict, field: str, value: object, message: str
) -> None:
    mutant = copy.deepcopy(payload)
    mutant[field] = value
    with pytest.raises(VerificationError, match=message):
        verify_payload(mutant)


def test_pinned_raw_input_sha_mutation_is_rejected(payload: dict) -> None:
    mutant = copy.deepcopy(payload)
    mutant["provenance"]["pinned_input_sha256"]["observer_views.jsonl"] = "0" * 64
    with pytest.raises(VerificationError, match="provenance metadata"):
        verify_payload(mutant)


def test_extra_top_level_field_is_rejected(payload: dict) -> None:
    mutant = copy.deepcopy(payload)
    mutant["physical_orientation"] = True
    with pytest.raises(VerificationError, match="top-level payload keys"):
        verify_payload(mutant)


def test_count_table_mutation_is_rejected(payload: dict) -> None:
    mutant = copy.deepcopy(payload)
    mutant["ordered_counts"][1][2] += 1
    with pytest.raises(VerificationError, match="pinned ordered count table"):
        verify_payload(mutant)


def test_transition_total_mutation_is_rejected(payload: dict) -> None:
    mutant = copy.deepcopy(payload)
    mutant["transition_count"] -= 1
    with pytest.raises(VerificationError, match="declared transition total"):
        verify_payload(mutant)


def test_antisymmetric_current_mutation_is_rejected(payload: dict) -> None:
    mutant = copy.deepcopy(payload)
    mutant["current_antisymmetric_part"][0][6] += 1
    with pytest.raises(VerificationError, match="antisymmetric current"):
        verify_payload(mutant)


def test_raw_pair_designation_mutation_is_rejected(payload: dict) -> None:
    mutant = copy.deepcopy(payload)
    mutant["designated_pair"]["count_forward"] -= 1
    with pytest.raises(VerificationError, match="serialized raw pair"):
        verify_payload(mutant)


def test_raw_cycle_designation_mutation_is_rejected(payload: dict) -> None:
    mutant = copy.deepcopy(payload)
    mutant["designated_cycle"]["buckets"] = [3, 5, 4]
    with pytest.raises(VerificationError, match="serialized raw cycle"):
        verify_payload(mutant)


def test_exact_row_normalized_maximizers_are_independently_recomputed() -> None:
    normalized_table = tuple(
        tuple(Fraction(value, sum(row)) for value in row) for row in EXPECTED_COUNTS
    )
    assert select_pair(normalized_table) == EXPECTED_NORMALIZED_PAIR
    assert select_cycle(normalized_table) == EXPECTED_NORMALIZED_CYCLE


def test_row_scaling_changes_raw_but_not_normalized_cycle_designation() -> None:
    mutant = [list(row) for row in EXPECTED_COUNTS]
    mutant[1] = [10_000 * value for value in mutant[1]]
    raw_cycle = select_cycle(mutant)
    assert raw_cycle is not None
    assert (raw_cycle.a, raw_cycle.b, raw_cycle.c) == (1, 2, 3)
    verify_normalized_designations(mutant)


def test_normalized_distribution_mutation_is_rejected() -> None:
    mutant = [list(row) for row in EXPECTED_COUNTS]
    mutant[2][3] -= 1
    mutant[2][5] += 1
    assert sum(mutant[2]) == sum(EXPECTED_COUNTS[2])
    with pytest.raises(VerificationError, match="normalized pair maximizer"):
        verify_normalized_designations(mutant)


def test_faulty_reversal_table_is_rejected() -> None:
    faulty_reverse = [list(row) for row in transpose(EXPECTED_COUNTS)]
    faulty_reverse[0][1] += 1
    with pytest.raises(VerificationError, match="not the transpose"):
        verify_transposition_semantics(EXPECTED_COUNTS, faulty_reverse)


def test_asymmetric_control_mutation_is_rejected() -> None:
    mutant = [list(row) for row in symmetric_control()]
    mutant[0][1] += 1
    with pytest.raises(VerificationError, match="not exactly symmetric"):
        verify_symmetric_control(mutant)
