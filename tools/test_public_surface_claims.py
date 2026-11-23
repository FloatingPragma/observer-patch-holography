from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import build_public_quantitative_section as builder
import public_surface_claims as claims


REPO_ROOT = Path(__file__).resolve().parents[1]
TABLE_README = "code/fixture_lane/README.md"
ALPHA_ROW = "alpha_gauge_width_diagnostic"
TABLE_HEADING = "## Diagnostic Values"
TABLE_TEXT = (
    "# Fixture lane\n\n"
    f"{TABLE_HEADING}\n\n"
    "| Observable | Branch value | Comparison coordinate |\n"
    "| --- | ---: | --- |\n"
    "| `alpha_inv` | `137.03599` (diagnostic) | `137.035999177` (CODATA 2022) |\n"
    "| `x_branch` | `1.0` (rejected candidate) | `1.1` (measured) |\n"
    "| `x_fit` | `2.0` (target-anchored fit; never a prediction) | `2.0` (measured) |\n"
    "| `y_label` | `n/a` | `n/a` |\n"
)
TABLE_DECLARATION = {
    "path": TABLE_README,
    "heading": TABLE_HEADING,
    "label_column": "Observable",
    "comparison_column": "Comparison coordinate",
    "rows": [
        {
            "label": "`alpha_inv`",
            "role": "diagnostic",
            "lane": {"kind": "manifest_row", "row_id": ALPHA_ROW},
        },
        {
            "label": "`x_branch`",
            "role": "rejected_candidate",
            "lane": {"kind": "registry_claim", "claim_id": "OPH-ALPHA-PIXEL"},
        },
        {
            "label": "`x_fit`",
            "role": "target_anchored_backsolve",
            "lane": {"kind": "registry_claim", "claim_id": "OPH-ALPHA-PIXEL"},
        },
        {"label": "`y_label`", "role": "conditional_branch"},
    ],
}


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
    for surface in manifest["comparison_table_surfaces"]:
        _copy(surface["path"], root)

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


def _table_fixture_root(tmp_path: Path) -> Path:
    """A fixture with one declared governed comparison table."""
    root = _fixture_root(tmp_path)
    table = root / TABLE_README
    table.parent.mkdir(parents=True, exist_ok=True)
    table.write_text(TABLE_TEXT, encoding="utf-8")
    _edit_manifest(
        root,
        lambda manifest: manifest["comparison_table_surfaces"].append(
            json.loads(json.dumps(TABLE_DECLARATION))
        ),
    )
    assert claims.check_repository(root) == []
    return root


def _edit_manifest(root: Path, mutate) -> None:
    path = root / claims.MANIFEST_RELATIVE
    manifest = json.loads(path.read_text(encoding="utf-8"))
    mutate(manifest)
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def _row(manifest: dict, row_id: str) -> dict:
    return next(row for row in manifest["rows"] if row["row_id"] == row_id)


def _surface(manifest: dict, path: str) -> dict:
    return next(
        surface
        for surface in manifest["comparison_table_surfaces"]
        if surface["path"] == path
    )


def _declared_row(manifest: dict, label: str) -> dict:
    return next(
        row
        for row in _surface(manifest, TABLE_README)["rows"]
        if row["label"] == label
    )


def _edit_surface_text(root: Path, relative: str, old: str, new: str) -> None:
    path = root / relative
    text = path.read_text(encoding="utf-8")
    assert text.count(old) == 1, old
    path.write_text(text.replace(old, new, 1), encoding="utf-8")


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
            manifest, ALPHA_ROW
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
            manifest, ALPHA_ROW
        ).update({"claim_class": "emitted_artifact"}),
    )
    issues = claims.check_repository(root)
    assert any("does not match registry class" in issue for issue in issues)


def test_readme_table_policy_drift_fails_closed(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    _edit_manifest(
        root,
        lambda manifest: manifest["scope_policy"].update(
            {"readme_table_condition": "always_render"}
        ),
    )
    issues = claims.check_repository(root)
    assert any("readme_table_condition must be" in issue for issue in issues)


def test_missing_emitter_and_unresolved_artifact_pointer_fail_closed(
    tmp_path,
) -> None:
    root = _fixture_root(tmp_path)
    _edit_manifest(
        root,
        lambda manifest: _row(manifest, ALPHA_ROW)["producer"].update(
            {"script": "code/P_derivation/missing_emitter.py"}
        ),
    )
    issues = claims.check_repository(root)
    assert any("emitting script does not exist" in issue for issue in issues)

    root = _fixture_root(tmp_path / "second")
    _edit_manifest(
        root,
        lambda manifest: _row(manifest, ALPHA_ROW)["producer"].update(
            {"value_pointer": "/modes/unregistered_number"}
        ),
    )
    issues = claims.check_repository(root)
    assert any("unresolved producer value" in issue for issue in issues)


def test_nonpromoted_claim_cannot_be_presented_as_an_oph_result(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    _edit_manifest(
        root,
        lambda manifest: _row(
            manifest, ALPHA_ROW
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
        root / "code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json"
    )
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    reference = json.loads(
        (root / "code/P_derivation/codata_2022_alpha_fixture.json").read_text(
            encoding="utf-8"
        )
    )["inverse_fine_structure_constant"]["value"]
    artifact["modes"]["thomson_structured_running_plus_gauge_width"][
        "fixed_point_point_estimate_display_only"
    ]["alpha_inv"] = reference
    artifact_path.write_text(json.dumps(artifact), encoding="utf-8")

    issues, *_ = claims.validate_manifest(root)
    assert any("self-comparison agrees within" in issue for issue in issues)


def test_target_anchored_near_match_requires_explicit_nonprediction_label(
    tmp_path,
) -> None:
    root = _fixture_root(tmp_path)

    def remove_label(manifest: dict) -> None:
        row = _row(manifest, ALPHA_ROW)
        row["role"] = "target_anchored_backsolve"
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
        row = _row(manifest, ALPHA_ROW)
        row["role"] = "rejected_candidate"
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


def test_audit_table_rows_are_not_misread_as_headers(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    audit = root / "docs" / "two_sided_audit.md"
    audit.parent.mkdir(parents=True)
    audit.write_text(
        "# Audit\n\n"
        "| Source statement | Finding |\n"
        "| --- | --- |\n"
        "| OPH predicts a measured value of 24 | The external comparison is "
        "diagnostic and supplies no physical claim. |\n",
        encoding="utf-8",
    )
    assert claims.check_repository(root) == []


def test_manual_generated_block_edit_is_detected(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    path = root / "README.md"
    marker = claims.BLOCK_END
    path.write_text(
        path.read_text(encoding="utf-8").replace(
            marker,
            "unregistered comparison result\n" + marker,
            1,
        ),
        encoding="utf-8",
    )
    issues = claims.check_repository(root)
    assert any("generated quantitative claim block is stale" in issue for issue in issues)


def test_root_cli_rejects_mutated_fixture(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    path = root / "README_FR.md"
    marker = claims.BLOCK_END
    path.write_text(
        path.read_text(encoding="utf-8").replace(
            marker,
            "résultat de comparaison non enregistré\n" + marker,
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


def test_zero_physical_establishment_suppresses_readme_tables(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    english = (root / "README.md").read_text(encoding="utf-8")
    french = (root / "README_FR.md").read_text(encoding="utf-8")
    for text in (english, french):
        assert "physical_establishment count is zero" in text
        assert "Quantitative Claim Status" not in text
        assert "Statut des énoncés quantitatifs" not in text
        assert "| Quantity | Registered branch value" not in text
        assert "| Quantité | Valeur de la branche" not in text
        assert "6.6742999959" not in text
        assert "299792458" not in text


def test_governed_readme_table_is_declared_row_for_row(tmp_path) -> None:
    root = _table_fixture_root(tmp_path)

    _edit_manifest(
        root,
        lambda manifest: _surface(manifest, TABLE_README)["rows"].pop(0),
    )
    issues = claims.check_repository(root)
    assert any("carries no declaration" in issue for issue in issues)

    root = _table_fixture_root(tmp_path / "second")
    _edit_manifest(
        root,
        lambda manifest: _declared_row(manifest, "`y_label`").update(
            {"label": "`y_label` (retired label)"}
        ),
    )
    issues = claims.check_repository(root)
    assert any(
        "declared row '`y_label` (retired label)' is absent from the rendered table"
        in issue
        for issue in issues
    )


def test_comparison_row_without_a_lane_must_drop_its_comparison(tmp_path) -> None:
    root = _table_fixture_root(tmp_path)

    def drop_lane(manifest: dict) -> None:
        _declared_row(manifest, "`alpha_inv`").pop("lane")

    _edit_manifest(root, drop_lane)
    issues = claims.check_repository(root)
    assert any(
        "must carry no numeral in its comparison cell" in issue for issue in issues
    )


@pytest.mark.parametrize(
    ("label", "lane", "diagnostic"),
    [
        (
            "`x_branch`",
            {"kind": "registry_claim", "claim_id": "OPH-GHOST-CLAIM"},
            "lane names unknown registry claim ID 'OPH-GHOST-CLAIM'",
        ),
        (
            "`alpha_inv`",
            {"kind": "manifest_row", "row_id": "ghost_manifest_row"},
            "lane names unknown manifest row ID 'ghost_manifest_row'",
        ),
        (
            "`x_fit`",
            {"kind": "prose_assertion", "row_id": "unsupported"},
            "lane kind must be one of",
        ),
    ],
)
def test_comparison_lane_must_resolve(tmp_path, label, lane, diagnostic) -> None:
    root = _table_fixture_root(tmp_path)
    _edit_manifest(
        root,
        lambda manifest: _declared_row(manifest, label).update({"lane": lane}),
    )
    issues = claims.check_repository(root)
    assert any(diagnostic in issue for issue in issues), issues


def test_comparison_role_wording_is_required(tmp_path) -> None:
    root = _table_fixture_root(tmp_path)
    _edit_surface_text(
        root,
        TABLE_README,
        "target-anchored fit; never a prediction",
        "target-anchored fit",
    )
    issues = claims.check_repository(root)
    assert any(
        "requires the row to state 'never a prediction'" in issue for issue in issues
    )

    root = _table_fixture_root(tmp_path / "second")
    _edit_manifest(
        root,
        lambda manifest: _declared_row(manifest, "`y_label`").update(
            {"role": "rejected_candidate"}
        ),
    )
    issues = claims.check_repository(root)
    assert any("requires the row to state 'rejected'" in issue for issue in issues)


def test_undeclared_readme_comparison_table_fails_closed(tmp_path) -> None:
    root = _fixture_root(tmp_path)
    undeclared = root / "code/lane/README.md"
    undeclared.parent.mkdir(parents=True, exist_ok=True)
    undeclared.write_text(
        "# Lane\n\n"
        "| Observable | Conditional value | Comparison coordinate |\n"
        "| --- | ---: | --- |\n"
        "| `m_x` | `1.234 GeV` | `1.200 GeV (measured)` |\n",
        encoding="utf-8",
    )
    issues = claims.check_repository(root)
    assert any(
        "code/lane/README.md:3: governed comparison table is not declared" in issue
        for issue in issues
    )


def test_declared_comparison_surface_must_exist(tmp_path) -> None:
    root = _table_fixture_root(tmp_path)
    (root / TABLE_README).unlink()
    issues = claims.check_repository(root)
    assert any(
        f"comparison table surface {TABLE_README}: declared surface does not "
        "exist" in issue
        for issue in issues
    )


def test_governed_table_heading_must_stay_in_place(tmp_path) -> None:
    root = _table_fixture_root(tmp_path)
    _edit_surface_text(root, TABLE_README, TABLE_HEADING, "## Values")
    issues = claims.check_repository(root)
    assert any(
        f"declared heading '{TABLE_HEADING}' is absent" in issue for issue in issues
    )


def test_separate_tables_in_one_readme_require_individual_declarations(tmp_path) -> None:
    root = _table_fixture_root(tmp_path)
    second = "## Held observations"
    surface = root / TABLE_README
    second_table = second + TABLE_TEXT.split(TABLE_HEADING, 1)[1]
    surface.write_text(TABLE_TEXT + "\n" + second_table, encoding="utf-8")
    assert any(
        "governed comparison table is not declared" in issue
        for issue in claims.check_repository(root)
    )

    def add_second(manifest):
        declaration = dict(_surface(manifest, TABLE_README), heading=second)
        manifest["comparison_table_surfaces"].append(declaration)

    _edit_manifest(root, add_second)
    assert claims.check_repository(root) == []
    _edit_manifest(root, add_second)
    assert any(
        "duplicate declaration for the same governed comparison table" in issue
        for issue in claims.check_repository(root)
    )


def test_colliding_table_selectors_cannot_cover_another_table(tmp_path) -> None:
    root = _table_fixture_root(tmp_path)
    def add_colliding_selector(manifest):
        duplicate = dict(_surface(manifest, TABLE_README), heading="# Fixture lane")
        manifest["comparison_table_surfaces"].append(duplicate)

    _edit_manifest(root, add_colliding_selector)
    assert any(
        "duplicate declaration for the same governed comparison table" in issue
        for issue in claims.check_repository(root)
    )
    (root / TABLE_README).write_text(
        TABLE_TEXT + "\n" + TABLE_HEADING + "\n", encoding="utf-8"
    )
    assert any("is ambiguous" in issue for issue in claims.check_repository(root))


def test_live_governed_rows_name_lanes_or_carry_no_comparison() -> None:
    """Every declared row resolves through a lane or shows no comparison value."""
    manifest = claims.load_json(
        REPO_ROOT / "claims/public_surface_quantitative_claims.json"
    )
    registry, _ = claims._registry_by_id(REPO_ROOT)
    manifest_rows = {row["row_id"]: row for row in manifest["rows"]}
    for surface in manifest["comparison_table_surfaces"]:
        text = (REPO_ROOT / surface["path"]).read_text(encoding="utf-8")
        offset = text.index(surface["heading"])
        table = next(t for t in claims.comparison_tables(text) if t["offset"] > offset)
        headers = table["headers"]
        comparison_index = headers.index(surface["comparison_column"])
        rendered = {cells[0]: cells for cells in table["rows"]}
        assert len(rendered) == len(surface["rows"])
        for row in surface["rows"]:
            cells = rendered[row["label"]]
            assert row["role"] in claims.ROLES
            lane = row.get("lane")
            if lane is None:
                assert not claims.NUMERIC_TOKEN.search(cells[comparison_index])
                continue
            assert lane["kind"] in claims.LANE_KINDS
            if lane["kind"] == "registry_claim":
                assert lane["claim_id"] in registry
            else:
                assert manifest_rows[lane["row_id"]]["role"] == row["role"]
