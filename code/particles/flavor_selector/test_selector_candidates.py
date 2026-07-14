"""Determinism, schema, and blindness tests for the flavor-orbit selector menu."""

import json
import math
import pathlib
import re
import sys

MODULE_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE_DIR))

import evaluate_candidates as ev  # noqa: E402

ARTIFACT = MODULE_DIR / "runtime" / "candidates_evaluated.json"
EXPECTED_IDS = [f"C-{i:02d}" for i in range(1, 13)]


def test_determinism_repeated_build():
    first = ev.serialize(ev.build_payload())
    second = ev.serialize(ev.build_payload())
    assert first == second


def test_determinism_artifact_matches_rebuild():
    assert ARTIFACT.exists(), "run evaluate_candidates.py first"
    on_disk = ARTIFACT.read_text(encoding="utf-8")
    assert on_disk == ev.serialize(ev.build_payload())


def test_no_timestamps_or_run_identifiers():
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    text = json.dumps(payload)
    for token in ("utc", "timestamp", "generated_at", "uuid"):
        assert token not in text.lower()


def test_schema_frozen_count_and_ids():
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert payload["frozen_candidate_count"] == 12
    ids = [c["id"] for c in payload["candidates"]]
    assert ids == EXPECTED_IDS


def test_schema_every_candidate_fixes_both_moduli():
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    for cand in payload["candidates"]:
        for key in ("id", "name", "closed_form", "moduli", "implied"):
            assert key in cand
        su = cand["moduli"]["sigma_u"]
        sd = cand["moduli"]["sigma_d"]
        assert isinstance(su, float) and su > 0.0
        assert isinstance(sd, float) and sd > 0.0
        assert cand["closed_form"]["sigma_u"]
        assert cand["closed_form"]["sigma_d"]


def test_schema_implied_ratio_block_consistent():
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    rho = payload["inputs"]["rho_ord_display_only"]
    for cand in payload["candidates"]:
        su = cand["moduli"]["sigma_u"]
        sd = cand["moduli"]["sigma_d"]
        imp = cand["implied"]
        assert math.isclose(imp["ln_yt_over_yu"], 2.0 * su, rel_tol=1e-12)
        assert math.isclose(imp["ln_yb_over_yd"], 2.0 * sd, rel_tol=1e-12)
        assert math.isclose(
            imp["ln_yc_over_yu"] + imp["ln_yt_over_yc"], 2.0 * su, rel_tol=1e-12
        )
        assert math.isclose(
            imp["ln_yc_over_yu"] / imp["ln_yt_over_yc"], rho, rel_tol=1e-12
        )
        assert math.isclose(
            imp["ln_yb_over_ys"] / imp["ln_ys_over_yd"], rho, rel_tol=1e-12
        )


def test_blindness_flags_and_status():
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    blind = payload["blindness"]
    assert blind["reads_measured_masses"] is False
    assert blind["reads_measured_ratios"] is False
    assert blind["reads_reference_artifacts"] is False
    assert payload["status"] == "evaluated_no_selection"


def test_blindness_evaluator_source_tokens():
    source = (MODULE_DIR / "evaluate_candidates.py").read_text(encoding="utf-8")
    lowered = source.lower()
    for token in ("pdg", "codata", "reference_values", "particle_reference"):
        assert token not in lowered
    # The evaluator reads no artifact from the runs surface.
    assert "runs/" not in source
    assert re.search(r"json\.load\b", source) is None
