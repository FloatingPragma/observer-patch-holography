"""Tests for the INS-03 export conformance validator.

Covered here: the shipped synthetic sample validates as
SCHEMA_CONFORMANT_SYNTHETIC and equals the builder's canonical text; each
required mutation guard fires its named error code (wrong effect entry,
non-summing context, wrong run-state literal, float outside
derived_for_display, wrong digest, missing field); the committed reference
constants match the literals of
Lean/EventAlgebra/LuedersPhaseInstrument.lean and
Lean/EventAlgebra/SourceBoundInstrumentInterface.lean, including the eight
outcome-0 entries 111/179, 111/179, 111/179, 315/716, 315/716, 315/716,
315/716, 1/2; and the synthetic/producer verdict distinction is exposed and
fail-closed.

What is not proved here: these tests exercise the validator against one
synthetic document and named mutations of it.  They certify no source
production, no provenance, and no custody, and they touch no register row.
"""

from __future__ import annotations

import copy
import json
import tempfile
from fractions import Fraction
from pathlib import Path

import ins03_export_validator as validator

SAMPLE_PATH = Path(__file__).parent / "sample_conforming_export.json"


def _load_sample() -> dict:
    return json.loads(SAMPLE_PATH.read_text(encoding="utf-8"))


def _validate_object(obj: dict, refresh_digest: bool = True) -> dict:
    """Write the object canonically to a scratch file and validate it."""
    if refresh_digest:
        obj["custody_digest_sha256"] = validator.compute_custody_digest(obj)
    with tempfile.TemporaryDirectory() as scratch:
        path = Path(scratch) / "export.json"
        path.write_text(validator.canonical_text(obj) + "\n", encoding="utf-8")
        return validator.validate_file(path)


def _codes(report: dict) -> set[str]:
    return {error["code"] for error in report["errors"]}


# ---------------------------------------------------------------------------
# The shipped sample
# ---------------------------------------------------------------------------


def test_sample_is_schema_conformant_synthetic() -> None:
    report = validator.validate_file(SAMPLE_PATH)
    assert report["errors"] == []
    assert report["verdict"] == "SCHEMA_CONFORMANT_SYNTHETIC"
    assert report["provenance_class"] == "synthetic"


def test_sample_file_equals_builder_canonical_text() -> None:
    assert SAMPLE_PATH.read_text(encoding="utf-8") == validator.sample_file_text()


def test_sample_declares_synthetic_provenance_class() -> None:
    sample = _load_sample()
    assert sample["provenance_class"] == "synthetic"
    assert sample["labels"]["evidential"] is False
    assert validator.SYNTHETIC_MARKER in sample["provenance"]["run_id"]


# ---------------------------------------------------------------------------
# Committed reference constants against the Lean literals
# ---------------------------------------------------------------------------


def test_committed_outcome0_entries_match_lean_table() -> None:
    """The eight outcome-0 entries of
    EventAlgebra.SourceBoundDeterminedData.publicTable_literal and
    EventAlgebra.luedersPhaseInstrument_run_table."""
    transcribed = (
        Fraction(111, 179), Fraction(111, 179), Fraction(111, 179),
        Fraction(315, 716), Fraction(315, 716), Fraction(315, 716),
        Fraction(315, 716), Fraction(1, 2),
    )
    assert validator.COMMITTED_OUTCOME0_FREQUENCIES == transcribed


def test_committed_run_state_matches_lean_literals() -> None:
    """EventAlgebra.committedRunState_eq_literal: diag(111/179, 68/179)."""
    assert validator.COMMITTED_RUN_STATE_DIAGONAL == (Fraction(111, 179),
                                                      Fraction(68, 179))


def test_committed_count_literals() -> None:
    """The fixture count literals: (111, 68) at mass 179 in the diagonal and
    record-conjugate contexts, (315, 401) at mass 716 in the rotated
    contexts (EventAlgebra.binaryFrequency_rotated_run), (179, 179) at mass
    358 in the phase context (EventAlgebra.modelFrequency_phase_zero)."""
    counts = validator.COMMITTED_COUNTS
    for context in ("web_diagonal", "web_conjugated_0", "web_conjugated_1"):
        assert counts[context] == (111, 68, 179)
    for context in ("web_conjugated_2", "web_conjugated_3",
                    "web_conjugated_4", "web_conjugated_5"):
        assert counts[context] == (315, 401, 716)
    assert counts["phase"] == (179, 179, 358)
    for context, (count0, _count1, mass) in counts.items():
        index = validator.COMMITTED_CONTEXTS.index(context)
        assert Fraction(count0, mass) == \
            validator.COMMITTED_OUTCOME0_FREQUENCIES[index]


def test_committed_offdiag_bound_matches_lean_literal() -> None:
    """EventAlgebra.prep_offdiag_normSq_le specialization: 7548/32041."""
    assert validator.COMMITTED_OFFDIAG_NORMSQ_BOUND == Fraction(7548, 32041)
    p, q = validator.COMMITTED_RUN_STATE_DIAGONAL
    assert p * q == Fraction(7548, 32041)


def test_committed_phase_effect_is_pauli_y_plus_projector() -> None:
    """EventAlgebra.sourcePhaseLift_entries: [[1/2, -i/2], [i/2, 1/2]]."""
    effect = validator.COMMITTED_EFFECT0["phase"]
    half = Fraction(1, 2)
    assert effect[0][0] == validator.c3_rational(half)
    assert effect[0][1] == validator.c3_imag(-half)
    assert effect[1][0] == validator.c3_imag(half)
    assert effect[1][1] == validator.c3_rational(half)


# ---------------------------------------------------------------------------
# Mutation guards
# ---------------------------------------------------------------------------


def test_mutation_wrong_effect_entry() -> None:
    obj = _load_sample()
    obj["outcome_maps"]["phase"]["0"]["declared_effect"][0][0] = \
        ["1/3", "0", "0", "0"]
    report = _validate_object(obj)
    assert report["verdict"] == "NONCONFORMANT"
    assert "EFFECT_DECLARED_MISMATCH" in _codes(report)


def test_mutation_non_summing_context() -> None:
    obj = _load_sample()
    obj["outcome_maps"]["web_diagonal"]["1"] = copy.deepcopy(
        obj["outcome_maps"]["web_diagonal"]["0"])
    obj["readback"]["web_diagonal"]["1"] = copy.deepcopy(
        obj["readback"]["web_diagonal"]["0"])
    report = _validate_object(obj)
    assert report["verdict"] == "NONCONFORMANT"
    assert "CONTEXT_SUM_NOT_IDENTITY" in _codes(report)


def test_mutation_wrong_run_state_literal() -> None:
    obj = _load_sample()
    prep = obj["preparation"]
    prep["rho_00"] = [110, 179]
    prep["rho_11"] = [69, 179]
    prep["positivity_certificate"]["diagonal_product"] = [7590, 32041]
    prep["preparation_content_sha256"] = \
        validator.compute_preparation_content_sha256(
            prep["carrier_id"], obj["provenance"]["simulator_commit"],
            prep["source_record_ids"], prep["rho_00"], prep["rho_01"],
            prep["rho_11"])
    report = _validate_object(obj)
    assert report["verdict"] == "NONCONFORMANT"
    assert "COMMITTED_RUN_STATE_MISMATCH" in _codes(report)
    assert "READBACK_TRACE_MISMATCH" in _codes(report)


def test_mutation_float_outside_display() -> None:
    obj = _load_sample()
    obj["provenance"]["import_graph_independence"]["edge_count"] = 0.5
    report = _validate_object(obj)
    assert report["verdict"] == "NONCONFORMANT"
    assert "FLOAT_OUTSIDE_DISPLAY" in _codes(report)


def test_float_inside_display_is_allowed() -> None:
    sample = _load_sample()
    display = sample["derived_for_display"]
    assert isinstance(display["phase_outcome0_frequency_float"], float)
    report = validator.validate_file(SAMPLE_PATH)
    assert "FLOAT_OUTSIDE_DISPLAY" not in _codes(report)


def test_mutation_wrong_digest() -> None:
    obj = _load_sample()
    obj["custody_digest_sha256"] = "0" * 64
    report = _validate_object(obj, refresh_digest=False)
    assert report["verdict"] == "NONCONFORMANT"
    assert "DIGEST_MISMATCH" in _codes(report)


def test_mutation_missing_field() -> None:
    obj = _load_sample()
    del obj["preparation"]
    report = _validate_object(obj)
    assert report["verdict"] == "NONCONFORMANT"
    assert "MISSING_FIELD" in _codes(report)


def test_mutation_missing_nested_field() -> None:
    obj = _load_sample()
    del obj["outcome_maps"]["phase"]["0"]["kraus"]
    report = _validate_object(obj)
    assert report["verdict"] == "NONCONFORMANT"
    assert "MISSING_FIELD" in _codes(report)


def test_noncanonical_serialization() -> None:
    obj = _load_sample()
    with tempfile.TemporaryDirectory() as scratch:
        path = Path(scratch) / "export.json"
        path.write_text(json.dumps(obj, sort_keys=True, indent=4) + "\n",
                        encoding="utf-8")
        report = validator.validate_file(path)
    assert report["verdict"] == "NONCONFORMANT"
    assert "NONCANONICAL_SERIALIZATION" in _codes(report)


def test_wrong_schema_id() -> None:
    obj = _load_sample()
    obj["schema"] = "oph.sim.ins03_phase_instrument_export.v0"
    report = _validate_object(obj)
    assert report["verdict"] == "NONCONFORMANT"
    assert "SCHEMA_ID_MISMATCH" in _codes(report)


# ---------------------------------------------------------------------------
# The synthetic/producer distinction
# ---------------------------------------------------------------------------


def test_producer_class_with_synthetic_markers_fails_closed() -> None:
    obj = _load_sample()
    obj["provenance_class"] = "producer"
    report = _validate_object(obj)
    assert report["verdict"] == "NONCONFORMANT"
    assert "SYNTHETIC_MARKER_IN_PRODUCER" in _codes(report)


def test_producer_classification_path_is_distinct_verdict() -> None:
    """A fabricated marker-free variant classifies as
    SCHEMA_CONFORMANT_PRODUCER: the verdict word differs from the synthetic
    one.  The inputs are fabricated test data; the validator certifies no
    provenance either way, per VALIDATOR_CONTRACT.md."""
    obj = _load_sample()
    obj["provenance_class"] = "producer"
    prep = obj["preparation"]
    prep["carrier_id"] = "carrier-under-test"
    prep["source_record_ids"] = ["record-under-test-0"]
    prep["operations"] = [{
        "operation": "unregisteredOp",
        "carrier_before": "carrier-under-test",
        "carrier_after": "carrier-under-test",
        "in_class": False,
        "description": "Fabricated test entry matching no constructor.",
    }]
    obj["provenance"] = {
        "producer_modules": [
            {"path": "oph_fpe/testdata/fabricated_module.py", "sha256": "a" * 64}
        ],
        "simulator_commit": "b" * 40,
        "repository_url": "https://github.com/FloatingPragma/oph-physics-sim",
        "rer_commit": "c" * 40,
        "run_id": "fabricated-test-run",
        "input_inventory": [],
        "runtime_read_log": [],
        "import_graph_independence": {
            "excluded_modules": list(validator.EXCLUDED_INPUT_FLOOR),
            "edge_count": 0,
            "dynamic_import_count": 0,
            "statement": "Fabricated test report.",
        },
    }
    prep["preparation_content_sha256"] = \
        validator.compute_preparation_content_sha256(
            prep["carrier_id"], obj["provenance"]["simulator_commit"],
            prep["source_record_ids"], prep["rho_00"], prep["rho_01"],
            prep["rho_11"])
    report = _validate_object(obj)
    assert report["errors"] == []
    assert report["verdict"] == "SCHEMA_CONFORMANT_PRODUCER"


def test_synthetic_evidential_conflict() -> None:
    obj = _load_sample()
    obj["labels"]["evidential"] = True
    report = _validate_object(obj)
    assert report["verdict"] == "NONCONFORMANT"
    assert "SYNTHETIC_EVIDENTIAL_CONFLICT" in _codes(report)


# ---------------------------------------------------------------------------
# Exact arithmetic guards
# ---------------------------------------------------------------------------


def test_q3_sign_is_exact() -> None:
    """sqrt(3) comparisons run on integers only: 7/4 > sqrt(3) > 433/250."""
    sqrt3 = validator.Q3(Fraction(0), Fraction(1))
    assert (validator.Q3(Fraction(7, 4)) - sqrt3).is_nonneg()
    assert not (validator.Q3(Fraction(433, 250)) - sqrt3).is_nonneg()
    assert validator.q3_le(sqrt3 * sqrt3, validator.Q3(Fraction(3)))
    assert validator.q3_le(validator.Q3(Fraction(3)), sqrt3 * sqrt3)


def test_positivity_gate_rejects_off_core_coordinate() -> None:
    """An off-diagonal coordinate outside |rho_01|^2 <= 7548/32041 fires the
    positivity error even with a consistent certificate structure."""
    obj = _load_sample()
    prep = obj["preparation"]
    prep["rho_01"] = ["1/2", "0", "1/2", "0"]
    prep["record_diagonal"] = False
    prep["record_diagonal_offdiag"] = ["1/2", "0", "1/2", "0"]
    prep["positivity_certificate"]["offdiag_norm_sq"] = ["1/2", "0"]
    prep["preparation_content_sha256"] = \
        validator.compute_preparation_content_sha256(
            prep["carrier_id"], obj["provenance"]["simulator_commit"],
            prep["source_record_ids"], prep["rho_00"], prep["rho_01"],
            prep["rho_11"])
    report = _validate_object(obj)
    assert report["verdict"] == "NONCONFORMANT"
    assert "PREP_POSITIVITY_VIOLATION" in _codes(report)


def test_mutation_renamed_committed_context_refused() -> None:
    obj = _load_sample()
    for section in ("outcome_maps", "readback"):
        obj[section]["renamed_context"] = obj[section].pop("phase")
    obj["contexts"] = [
        "renamed_context" if name == "phase" else name
        for name in obj["contexts"]
    ]
    report = _validate_object(obj)
    assert report["verdict"] == "NONCONFORMANT"
    assert "PARTIAL_COMMITTED_CONTEXT_SET" in _codes(report)
