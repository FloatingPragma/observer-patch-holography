from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import build_public_quantitative_section as builder
import public_surface_claims as claims


REPO_ROOT = Path(__file__).resolve().parents[1]


def _copy(relative: str, target_root: Path) -> None:
    source = REPO_ROOT / relative
    target = target_root / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def _fixture_root(tmp_path: Path) -> Path:
    root = tmp_path / "repo"
    root.mkdir(parents=True)
    _copy(str(claims.MANIFEST_RELATIVE), root)
    _copy(str(claims.REGISTRY_RELATIVE), root)

    manifest = json.loads(
        (root / claims.MANIFEST_RELATIVE).read_text(encoding="utf-8")
    )
    needed: set[str] = set()
    for row in manifest["rows"]:
        needed.add(row["producer"]["script"])
        needed.add(row["producer"]["artifact"])
        if "reference" in row:
            needed.add(row["reference"]["artifact"])
        for support in row.get("supporting_artifacts", []):
            needed.add(support["artifact"])
    for relative in sorted(needed):
        _copy(relative, root)

    for surface in manifest["surfaces"]:
        path = root / surface["path"]
        path.write_text(
            "# Fixture\n\n"
            f"{claims.BLOCK_START}\n"
            f"{claims.BLOCK_END}\n",
            encoding="utf-8",
        )
    assert builder.build(root, check=False) == []
    return root


def _edit_manifest(root: Path, mutate) -> None:
    path = root / claims.MANIFEST_RELATIVE
    manifest = json.loads(path.read_text(encoding="utf-8"))
    mutate(manifest)
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def _row(manifest: dict, row_id: str) -> dict:
    return next(row for row in manifest["rows"] if row["row_id"] == row_id)


def test_clean_generated_fixture_passes_and_is_deterministic(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    assert claims.check_repository(root) == []
    assert builder.build(root, check=True) == []

    first = (root / "README.md").read_bytes()
    assert builder.build(root, check=False) == []
    assert (root / "README.md").read_bytes() == first


def test_unknown_claim_and_class_drift_fail_closed(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    _edit_manifest(
        root,
        lambda manifest: _row(
            manifest, "bottom_quark_clebsch"
        ).update(
            {
                "claim_id": "OPH-GHOST-CLAIM",
                "claim_class": "emitted_artifact",
            }
        ),
    )
    issues = claims.check_repository(root)
    assert any("unknown registry claim ID" in issue for issue in issues)

    root = _fixture_root(tmp_path / "second")
    _edit_manifest(
        root,
        lambda manifest: _row(
            manifest, "bottom_quark_clebsch"
        ).update({"claim_class": "emitted_artifact"}),
    )
    issues = claims.check_repository(root)
    assert any("does not match registry class" in issue for issue in issues)


def test_missing_emitter_and_unresolved_artifact_pointer_fail_closed(
    tmp_path,
) -> None:
    root = _fixture_root(tmp_path)
    _edit_manifest(
        root,
        lambda manifest: _row(manifest, "bottom_quark_clebsch")["producer"].update(
            {"script": "code/particles/flavor/missing_emitter.py"}
        ),
    )
    issues = claims.check_repository(root)
    assert any("emitting script does not exist" in issue for issue in issues)

    root = _fixture_root(tmp_path / "second")
    _edit_manifest(
        root,
        lambda manifest: _row(manifest, "bottom_quark_clebsch")["producer"].update(
            {"value_pointer": "/predictions/unregistered_number"}
        ),
    )
    issues = claims.check_repository(root)
    assert any("unresolved producer value" in issue for issue in issues)


def test_nonpromoted_claim_cannot_be_presented_as_an_oph_result(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    _edit_manifest(
        root,
        lambda manifest: _row(
            manifest, "bottom_quark_clebsch"
        ).update({"role": "oph_result"}),
    )
    issues = claims.check_repository(root)
    assert any(
        "OPH result value requires physical_establishment or empirical_implementation"
        in issue
        for issue in issues
    )


def test_self_comparison_fails_without_target_anchored_role(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    artifact_path = (
        root / "code/particles/runs/flavor/down_type_register_clebsch_lane.json"
    )
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    artifact["predictions"]["mb_mb_gev"] = artifact["compare_only"]["references"][
        "mb_mb_gev"
    ]
    artifact_path.write_text(json.dumps(artifact), encoding="utf-8")

    issues, *_ = claims.validate_manifest(root)
    assert any("self-comparison agrees within" in issue for issue in issues)


def test_target_anchored_near_match_requires_explicit_nonprediction_label(
    tmp_path,
) -> None:
    root = _fixture_root(tmp_path)

    def remove_label(manifest: dict) -> None:
        row = _row(manifest, "top_target_anchored")
        row["status"]["en"] = "Back-solved from the measured pair."

    _edit_manifest(root, remove_label)
    issues = claims.check_repository(root)
    assert any(
        "target-anchored back-solves must say they are never a prediction" in issue
        for issue in issues
    )


def test_rejected_candidate_requires_explicit_disposition_label(tmp_path) -> None:
    root = _fixture_root(tmp_path)

    def remove_label(manifest: dict) -> None:
        row = _row(manifest, "bottom_quark_clebsch")
        row["status"]["en"] = "Conditional comparison output."

    _edit_manifest(root, remove_label)
    issues = claims.check_repository(root)
    assert any(
        "rejected candidates must state the rejection disposition" in issue
        for issue in issues
    )


def test_alpha_endpoint_definition_and_issue_cap_are_source_bound(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    calibration_path = (
        root / "code/P_derivation/runtime/measured_endpoint_calibration_current.json"
    )
    calibration = json.loads(calibration_path.read_text(encoding="utf-8"))
    calibration["calibrated_values"]["definition"] = "unregistered definition"
    calibration_path.write_text(json.dumps(calibration), encoding="utf-8")
    issues = claims.check_repository(root)
    assert any(
        "supporting artifact 0" in issue and "expected" in issue for issue in issues
    )

    root = _fixture_root(tmp_path / "second")
    bridge_path = root / "code/P_derivation/runtime/anchor_scheme_bridge_current.json"
    bridge = json.loads(bridge_path.read_text(encoding="utf-8"))
    bridge["verdict"]["source_only_reduction"] = "untracked blocker"
    bridge_path.write_text(json.dumps(bridge), encoding="utf-8")
    issues = claims.check_repository(root)
    assert any(
        "supporting artifact 1" in issue and "does not contain '#425'" in issue
        for issue in issues
    )


def test_unmanaged_numeric_oph_external_table_fails_closed(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    bad = root / "docs" / "bad_table.md"
    bad.parent.mkdir(parents=True)
    bad.write_text(
        "# Bad table\n\n"
        "| Quantity | OPH | PDG reference |\n"
        "| --- | ---: | ---: |\n"
        "| mass | 172.3523553288312 | 172.1 |\n",
        encoding="utf-8",
    )
    issues = claims.check_repository(root)
    assert any("unmanaged numeric OPH-versus-external" in issue for issue in issues)

    bad.write_text(
        "# Compact bad table\n\n"
        "|Quantity|OPH|NIST|\n"
        "|---|---:|---:|\n"
        "|constant|137.035999177|137.035999177|\n",
        encoding="utf-8",
    )
    issues = claims.check_repository(root)
    assert any("unmanaged numeric OPH-versus-external" in issue for issue in issues)


def test_prose_numerals_are_not_false_positive_claims(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    context = root / "docs" / "numeric_context.md"
    context.parent.mkdir(parents=True)
    context.write_text(
        "# Context\n\n"
        "The 2026-07-25 audit references issue #425, more than 800 theorem and "
        "lemma declarations, `code/D10/path_2.py`, and the representation "
        "$A_5\\times Z_6$. These are provenance or structural notation rather "
        "than an OPH-versus-external quantitative table.\n",
        encoding="utf-8",
    )
    assert claims.check_repository(root) == []


def test_manual_generated_block_edit_is_detected(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    path = root / "README.md"
    path.write_text(
        path.read_text(encoding="utf-8").replace(
            "`6.03 GeV`",
            "`4.18 GeV`",
            1,
        ),
        encoding="utf-8",
    )
    issues = claims.check_repository(root)
    assert any("generated quantitative claim block is stale" in issue for issue in issues)


def test_root_cli_rejects_mutated_fixture(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    path = root / "README_FR.md"
    path.write_text(
        path.read_text(encoding="utf-8").replace(
            "`140 MeV`",
            "`93,5 MeV`",
            1,
        ),
        encoding="utf-8",
    )
    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "tools/check_public_surface_claims.py"),
            "--root",
            str(root),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "generated quantitative claim block is stale" in result.stdout


def test_generated_block_has_required_boundaries(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    english = (root / "README.md").read_text(encoding="utf-8")
    french = (root / "README_FR.md").read_text(encoding="utf-8")
    for text in (english, french):
        assert "0 claims of class `physical_establishment`" in text or (
            "0 énoncé de classe `physical_establishment`" in text
        )
        assert "6.6742999959" not in text
        assert "299792458" not in text
    assert "target-anchored fit" in english
    assert "never a prediction" in english
    assert "structural" in english.lower()
    assert "non-discriminating" in english.lower()
    assert "1.6\\times10^{4}" in english
    assert "$G_{\\rm geom}/\\ell_\\star^2$" in english
    assert "unit-bookkeeping identity, not an SI prediction" in english
    assert "no $G_{\\rm SI}$ value or sigma distance" in english


def test_rejected_clebsch_rows_carry_rejected_candidate_role() -> None:
    manifest = claims.load_json(
        REPO_ROOT / "claims/public_surface_quantitative_claims.json"
    )
    clebsch_rows = [
        row for row in manifest["rows"] if "clebsch" in row["row_id"]
    ]
    assert len(clebsch_rows) == 4
    assert {row["role"] for row in clebsch_rows} == {"rejected_candidate"}
    assert {row["claim_id"] for row in clebsch_rows} == {
        "OPH-QUARK-REGISTER-CLEBSCH"
    }
    assert {
        row["producer"]["guards"]["/status"] for row in clebsch_rows
    } == {"CONDITIONAL_DECLARED_ROUTE_RETROSPECTIVELY_REJECTED"}
