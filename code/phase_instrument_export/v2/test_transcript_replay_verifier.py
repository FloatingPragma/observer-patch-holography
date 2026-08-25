"""Tests for the INS-03 v2 transcript replay verifier.

Covered here: the shipped synthetic sample transcript replays and verifies
as TRANSCRIPT_REPLAY_VERIFIED_UNAUTHENTICATED and derives the committed
static fixture byte-for-byte; the target-blind mutation controls each
produce their named fail-closed error (a mutated argument, a wrong
intermediate digest, a reordered step, a deleted step, a wrong final
export, and a producer-claim field); the committed reference constants
seen through v2 are the v1 validator's, checked by identity and by
independent literals; and the verdict vocabulary extends v1 exactly.

The mutation controls are target-blind: each mutation is chosen by
structural position (a step index, a field name), not by any endpoint
value, and each control asserts a named failure code rather than a value
comparison against any physical target.

What is not proved here: these tests exercise the verifier against one
synthetic transcript and named mutations of it.  They certify no source
production, no provenance, and no custody; the self-digest refresh test
demonstrates the opposite, that a recomputed self-digest authenticates
nothing.  No register row is touched; PR-03, PR-64, and PR-65 are open.
"""

from __future__ import annotations

import copy
import json
import tempfile
from fractions import Fraction
from pathlib import Path

import transcript_replay_verifier as trv

v1 = trv.v1

SAMPLE_PATH = Path(__file__).parent / "sample_synthetic_transcript.json"
V1_SAMPLE_PATH = Path(__file__).parent.parent / "sample_conforming_export.json"

# Structural positions in the shipped sample transcript: one load, one
# prepare, eight assembles, eight readouts, six binds, one emit.
STEP_LOAD = 0
STEP_PREPARE = 1
STEP_FIRST_ASSEMBLE = 2
STEP_FIRST_READOUT = 10
STEP_LAST_READOUT = 17
STEP_EMIT = 24


def _load_sample() -> dict:
    return json.loads(SAMPLE_PATH.read_text(encoding="utf-8"))


def _refresh_transcript_digest(transcript: dict) -> None:
    transcript["transcript_digest_sha256"] = \
        trv.compute_transcript_digest(transcript)


def _refresh_claimed_export_digests(transcript: dict) -> None:
    export = transcript["claimed_export"]
    export["custody_digest_sha256"] = v1.compute_custody_digest(export)
    transcript["claimed_export_digest_sha256"] = \
        trv.compute_claimed_export_digest(export)


def _verify_object(transcript: dict, refresh_digest: bool = True) -> dict:
    """Write the transcript canonically to a scratch file and verify it.

    The default digest refresh mirrors what any editor of the file can do:
    the self-digest is recomputable by whoever holds the bytes, so
    refreshing it before verification shows that it protects nothing."""
    if refresh_digest:
        _refresh_transcript_digest(transcript)
    with tempfile.TemporaryDirectory() as scratch:
        path = Path(scratch) / "transcript.json"
        path.write_text(v1.canonical_text(transcript) + "\n",
                        encoding="utf-8")
        return trv.verify_file(path)


def _codes(report: dict) -> set[str]:
    return {error["code"] for error in report["errors"]}


# ---------------------------------------------------------------------------
# The shipped sample
# ---------------------------------------------------------------------------


def test_sample_transcript_replays_and_verifies() -> None:
    report = trv.verify_file(SAMPLE_PATH)
    assert report["errors"] == []
    assert report["verdict"] == trv.REPLAY_VERDICT
    assert report["transcript_class"] == "synthetic"
    assert report["class_marker"] == trv.SYNTHETIC_CLASS_MARKER
    assert report["semantics_status"] == trv.SEMANTICS_STATUS_PROPOSED


def test_sample_file_equals_builder_canonical_text() -> None:
    assert SAMPLE_PATH.read_text(encoding="utf-8") \
        == trv.sample_transcript_file_text()


def test_sample_derives_committed_static_fixture() -> None:
    """The claimed export of the sample transcript is byte-identical to the
    committed v1 static fixture, both as the builder object and as the
    committed sample_conforming_export.json file bytes."""
    sample = _load_sample()
    fixture = v1.build_synthetic_sample()
    assert v1.canonical_text(sample["claimed_export"]) \
        == v1.canonical_text(fixture)
    assert v1.canonical_text(sample["claimed_export"]) + "\n" \
        == V1_SAMPLE_PATH.read_text(encoding="utf-8")


def test_sample_is_marked_synthetic_everywhere() -> None:
    sample = _load_sample()
    assert sample["transcript_class"] == "synthetic"
    assert sample["claimed_export"]["provenance_class"] == "synthetic"
    assert sample["claimed_export"]["labels"]["evidential"] is False
    assert v1.SYNTHETIC_MARKER in \
        sample["claimed_export"]["provenance"]["run_id"]


def test_sample_report_carries_boundary_and_v1_reading() -> None:
    report = trv.verify_file(SAMPLE_PATH)
    assert "not production" in report["boundary"]
    assert "not custody" in report["boundary"]
    assert report["v1_static_verdict"] == v1.STATIC_FIXTURE_VERDICT
    assert "not a replay gate" in report["v1_static_verdict_note"]


# ---------------------------------------------------------------------------
# Verdict vocabulary and committed reference constants
# ---------------------------------------------------------------------------


def test_verdict_vocabulary_extends_v1_exactly() -> None:
    assert trv.REPLAY_VERDICT \
        == "TRANSCRIPT_REPLAY_VERIFIED_UNAUTHENTICATED"
    assert trv.NONCONFORMANT_VERDICT == "NONCONFORMANT"
    assert trv.SYNTHETIC_CLASS_MARKER == "SYNTHETIC"
    assert v1.STATIC_FIXTURE_VERDICT == "STATIC_COMMITTED_FIXTURE_CONFORMANT"


def test_reference_constants_are_the_v1_validators() -> None:
    """Identity through the single-source import, plus independent literal
    checks against the committed values the v1 docstring cites."""
    assert trv.COMMITTED_CONTEXTS is v1.COMMITTED_CONTEXTS
    assert trv.COMMITTED_EFFECT0 is v1.COMMITTED_EFFECT0
    assert trv.COMMITTED_COUNTS is v1.COMMITTED_COUNTS
    assert trv.COMMITTED_RUN_STATE_DIAGONAL is v1.COMMITTED_RUN_STATE_DIAGONAL
    assert trv.COMMITTED_OFFDIAG_NORMSQ_BOUND \
        is v1.COMMITTED_OFFDIAG_NORMSQ_BOUND

    assert trv.COMMITTED_CONTEXTS == (
        "web_diagonal", "web_conjugated_0", "web_conjugated_1",
        "web_conjugated_2", "web_conjugated_3", "web_conjugated_4",
        "web_conjugated_5", "phase")
    assert trv.COMMITTED_COUNTS["web_diagonal"] == (111, 68, 179)
    assert trv.COMMITTED_COUNTS["web_conjugated_3"] == (315, 401, 716)
    assert trv.COMMITTED_COUNTS["phase"] == (179, 179, 358)
    assert trv.COMMITTED_RUN_STATE_DIAGONAL \
        == (Fraction(111, 179), Fraction(68, 179))
    assert trv.COMMITTED_OFFDIAG_NORMSQ_BOUND == Fraction(7548, 32041)


def test_sample_step_positions_match_structure() -> None:
    sample = _load_sample()
    steps = sample["steps"]
    assert len(steps) == 25
    assert steps[STEP_LOAD]["primitive"] == "load_committed_effect_table"
    assert steps[STEP_PREPARE]["primitive"] == "prepare_state"
    assert steps[STEP_FIRST_ASSEMBLE]["primitive"] == "assemble_lueders_kraus"
    assert steps[STEP_FIRST_READOUT]["primitive"] == "readout_table"
    assert steps[STEP_LAST_READOUT]["primitive"] == "readout_table"
    assert steps[STEP_LAST_READOUT]["args"]["context"] == "phase"
    assert steps[STEP_EMIT]["primitive"] == "emit_export"


# ---------------------------------------------------------------------------
# Target-blind mutation controls, each with its named failure
# ---------------------------------------------------------------------------


def test_mutated_argument_fails_with_step_digest_mismatch() -> None:
    """Mutating a readout count argument (keeping the count sum lawful)
    diverges the replayed state from the declared intermediate digest at
    that step; the named failure is STEP_DIGEST_MISMATCH there."""
    transcript = copy.deepcopy(_load_sample())
    args = transcript["steps"][STEP_FIRST_READOUT]["args"]
    args["count_0"] -= 1
    args["count_1"] += 1
    report = _verify_object(transcript)
    assert report["verdict"] == trv.NONCONFORMANT_VERDICT
    assert "STEP_DIGEST_MISMATCH" in _codes(report)
    assert any(error["code"] == "STEP_DIGEST_MISMATCH"
               and error["where"] == f"steps[{STEP_FIRST_READOUT}]"
               for error in report["errors"])
    assert "FINAL_EXPORT_MISMATCH" in _codes(report)


def test_wrong_intermediate_digest_fails_at_named_step() -> None:
    transcript = copy.deepcopy(_load_sample())
    transcript["steps"][5]["state_digest_sha256"] = "ab" * 32
    report = _verify_object(transcript)
    assert report["verdict"] == trv.NONCONFORMANT_VERDICT
    mismatches = [error for error in report["errors"]
                  if error["code"] == "STEP_DIGEST_MISMATCH"]
    assert [error["where"] for error in mismatches] == ["steps[5]"]


def test_reordered_step_fails_with_precondition() -> None:
    """Moving a readout before its assemble violates the step order; the
    named failure is STEP_PRECONDITION_UNMET at the early readout."""
    transcript = copy.deepcopy(_load_sample())
    steps = transcript["steps"]
    steps[STEP_FIRST_ASSEMBLE], steps[STEP_FIRST_READOUT] = \
        steps[STEP_FIRST_READOUT], steps[STEP_FIRST_ASSEMBLE]
    for position, step in enumerate(steps):
        step["index"] = position
    report = _verify_object(transcript)
    assert report["verdict"] == trv.NONCONFORMANT_VERDICT
    assert any(error["code"] == "STEP_PRECONDITION_UNMET"
               and error["where"] == f"steps[{STEP_FIRST_ASSEMBLE}]"
               for error in report["errors"])


def test_deleted_step_fails_with_emit_precondition() -> None:
    """Deleting the phase readout and renumbering leaves the emission
    without full context coverage; the named failure is
    EMIT_PRECONDITION_UNMET."""
    transcript = copy.deepcopy(_load_sample())
    del transcript["steps"][STEP_LAST_READOUT]
    for position, step in enumerate(transcript["steps"]):
        step["index"] = position
    report = _verify_object(transcript)
    assert report["verdict"] == trv.NONCONFORMANT_VERDICT
    assert "EMIT_PRECONDITION_UNMET" in _codes(report)


def test_wrong_final_export_fails_with_final_export_mismatch() -> None:
    """Editing the claimed export while leaving the steps intact (with
    every self-digest refreshed, as any editor can) is caught only by the
    replay comparison; the named failure is FINAL_EXPORT_MISMATCH."""
    transcript = copy.deepcopy(_load_sample())
    readback = transcript["claimed_export"]["readback"]["web_diagonal"]
    readback["0"]["count"] -= 1
    readback["1"]["count"] += 1
    _refresh_claimed_export_digests(transcript)
    report = _verify_object(transcript)
    assert report["verdict"] == trv.NONCONFORMANT_VERDICT
    assert "FINAL_EXPORT_MISMATCH" in _codes(report)
    assert "STEP_DIGEST_MISMATCH" not in _codes(report)
    assert "CLAIMED_EXPORT_DIGEST_MISMATCH" not in _codes(report)


def test_producer_class_fails_closed() -> None:
    transcript = copy.deepcopy(_load_sample())
    transcript["transcript_class"] = "producer"
    report = _verify_object(transcript)
    assert report["verdict"] == trv.NONCONFORMANT_VERDICT
    assert "PRODUCER_AUTHENTICATION_UNIMPLEMENTED" in _codes(report)
    assert report["class_marker"] is None


def test_authenticated_binding_hook_is_unimplemented_and_fails_closed() -> None:
    transcript = copy.deepcopy(_load_sample())
    transcript["authenticated_binding"] = {
        "signature": "00" * 32,
        "public_key": "00" * 16,
    }
    report = _verify_object(transcript)
    assert report["verdict"] == trv.NONCONFORMANT_VERDICT
    assert any(error["code"] == "PRODUCER_AUTHENTICATION_UNIMPLEMENTED"
               and error["where"] == "authenticated_binding"
               for error in report["errors"])


def test_evidential_label_claim_fails_closed() -> None:
    transcript = copy.deepcopy(_load_sample())
    transcript["claimed_export"]["labels"]["evidential"] = True
    _refresh_claimed_export_digests(transcript)
    report = _verify_object(transcript)
    assert report["verdict"] == trv.NONCONFORMANT_VERDICT
    assert "PRODUCER_AUTHENTICATION_UNIMPLEMENTED" in _codes(report)


def test_producer_provenance_class_in_export_fails_closed() -> None:
    transcript = copy.deepcopy(_load_sample())
    transcript["claimed_export"]["provenance_class"] = "producer"
    _refresh_claimed_export_digests(transcript)
    report = _verify_object(transcript)
    assert report["verdict"] == trv.NONCONFORMANT_VERDICT
    assert "PRODUCER_AUTHENTICATION_UNIMPLEMENTED" in _codes(report)


# ---------------------------------------------------------------------------
# Grammar and format failures
# ---------------------------------------------------------------------------


def test_unknown_primitive_named_failure() -> None:
    transcript = copy.deepcopy(_load_sample())
    transcript["steps"][STEP_LOAD]["primitive"] = "conjure_state"
    report = _verify_object(transcript)
    assert report["verdict"] == trv.NONCONFORMANT_VERDICT
    assert any(error["code"] == "UNKNOWN_PRIMITIVE"
               and error["where"] == f"steps[{STEP_LOAD}]"
               for error in report["errors"])


def test_non_string_primitive_fails_closed_without_exception() -> None:
    """A JSON-valid but unhashable primitive name must return a report.

    This is the malformed-content boundary promised by
    ``verify_transcript_text``; previously ``primitive = []`` escaped as a
    ``TypeError`` from the executor dictionary lookup.
    """
    transcript = copy.deepcopy(_load_sample())
    transcript["steps"][STEP_LOAD]["primitive"] = []
    report = _verify_object(transcript)
    assert report["verdict"] == trv.NONCONFORMANT_VERDICT
    assert any(
        error["code"] == "UNKNOWN_PRIMITIVE"
        and error["where"] == f"steps[{STEP_LOAD}]"
        and "must be a string" in error["detail"]
        for error in report["errors"]
    )


def test_argument_outside_grammar_named_failure() -> None:
    """A preparation whose diagonal does not sum to one is outside the
    grammar's domain; the named failure is ARGUMENT_OUTSIDE_GRAMMAR."""
    transcript = copy.deepcopy(_load_sample())
    transcript["steps"][STEP_PREPARE]["args"]["rho_00"] = [112, 179]
    report = _verify_object(transcript)
    assert report["verdict"] == trv.NONCONFORMANT_VERDICT
    assert any(error["code"] == "ARGUMENT_OUTSIDE_GRAMMAR"
               and error["where"] == f"steps[{STEP_PREPARE}]"
               for error in report["errors"])


def test_noncanonical_serialization_named_failure() -> None:
    transcript = _load_sample()
    with tempfile.TemporaryDirectory() as scratch:
        path = Path(scratch) / "transcript.json"
        path.write_text(json.dumps(transcript) + "\n", encoding="utf-8")
        report = trv.verify_file(path)
    assert report["verdict"] == trv.NONCONFORMANT_VERDICT
    assert "NONCANONICAL_SERIALIZATION" in _codes(report)


def test_frozen_semantics_claim_fails_closed() -> None:
    """No freeze registration exists, so a claimed frozen status is
    unverifiable and fails closed with SEMANTICS_STATUS_INVALID."""
    transcript = copy.deepcopy(_load_sample())
    transcript["semantics_status"] = "FROZEN"
    report = _verify_object(transcript)
    assert report["verdict"] == trv.NONCONFORMANT_VERDICT
    assert "SEMANTICS_STATUS_INVALID" in _codes(report)


def test_stale_transcript_digest_is_detected_when_not_refreshed() -> None:
    transcript = copy.deepcopy(_load_sample())
    transcript["steps"][STEP_FIRST_READOUT]["args"]["count_0"] -= 1
    transcript["steps"][STEP_FIRST_READOUT]["args"]["count_1"] += 1
    report = _verify_object(transcript, refresh_digest=False)
    assert report["verdict"] == trv.NONCONFORMANT_VERDICT
    assert "TRANSCRIPT_DIGEST_MISMATCH" in _codes(report)


def test_self_digest_refresh_authenticates_nothing() -> None:
    """The custody-is-not-self-digest ruling, executable: refreshing every
    self-digest after a mutation still fails replay, and a fully
    self-consistent transcript still cannot reach any producer verdict,
    because no producer verdict exists in this verifier."""
    transcript = copy.deepcopy(_load_sample())
    transcript["steps"][STEP_FIRST_READOUT]["args"]["count_0"] -= 1
    transcript["steps"][STEP_FIRST_READOUT]["args"]["count_1"] += 1
    _refresh_transcript_digest(transcript)
    report = _verify_object(transcript, refresh_digest=False)
    assert report["verdict"] == trv.NONCONFORMANT_VERDICT
    assert "TRANSCRIPT_DIGEST_MISMATCH" not in _codes(report)
    assert "STEP_DIGEST_MISMATCH" in _codes(report)


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def test_cli_exit_codes() -> None:
    assert trv.main([str(SAMPLE_PATH)]) == 0
    with tempfile.TemporaryDirectory() as scratch:
        missing = Path(scratch) / "missing.json"
        assert trv.main([str(missing)]) == 2
        garbled = Path(scratch) / "garbled.json"
        garbled.write_text("{not json", encoding="utf-8")
        assert trv.main([str(garbled)]) == 2
        transcript = copy.deepcopy(_load_sample())
        transcript["transcript_class"] = "producer"
        _refresh_transcript_digest(transcript)
        producer_path = Path(scratch) / "producer.json"
        producer_path.write_text(v1.canonical_text(transcript) + "\n",
                                 encoding="utf-8")
        assert trv.main([str(producer_path)]) == 1
