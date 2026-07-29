"""Fail-closed tests for the #607 frozen-prediction ladder register."""

import importlib.util
import json
import shutil
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
BUILDER = REPO_ROOT / "tools" / "build_fz_registry.py"

spec = importlib.util.spec_from_file_location("build_fz_registry", BUILDER)
fz_tool = importlib.util.module_from_spec(spec)
sys.modules["build_fz_registry"] = fz_tool
spec.loader.exec_module(fz_tool)


def live_register() -> dict:
    return json.loads(fz_tool.REGISTER_PATH.read_text(encoding="utf-8"))


def test_live_register_validates_and_surface_is_current():
    register = live_register()
    rows = fz_tool.validate(register)
    rendered = fz_tool.render(register, rows)
    committed = fz_tool.SURFACE_PATH.read_text(encoding="utf-8")
    assert rendered == committed


def test_ladder_excludes_the_retrospective_fz04_reservation():
    rows = fz_tool.validate(live_register())
    assert [row["id"] for row in rows] == [
        "FZ-01",
        "FZ-02",
        "FZ-03",
        "FZ-05",
        "FZ-06",
        "FZ-07",
        "FZ-08",
        "FZ-09",
        "FZ-10",
    ]
    result = live_register()["retrospective_results"][0]
    assert result["id"] == "RR-506-ALPHA-HVP"
    assert result["former_ladder_reservation"] == "FZ-04"


def test_fz06_is_attested_history_without_prediction_eligibility():
    register = live_register()
    fz06 = next(row for row in register["rows"] if row["id"] == "FZ-06")
    assert fz06["status"] == "superseded_void"
    assert fz06["frozen_utc"] is None
    assert fz06["comparison_protocol"].lower().startswith("none;")
    assert fz06["kill_band"].lower().startswith("none;")
    assert "alpha = 4" in fz06["content"]
    assert "void" in fz06["content"].lower()
    fz_tool.validate(register)

    register = live_register()
    fz06 = next(row for row in register["rows"] if row["id"] == "FZ-06")
    fz06["comparison_protocol"] = "score every loud event"
    with pytest.raises(SystemExit, match="must refuse comparison"):
        fz_tool.validate(register)

    register = live_register()
    fz06 = next(row for row in register["rows"] if row["id"] == "FZ-06")
    fz06["frozen_utc"] = "2026-07-17T07:18:00Z"
    with pytest.raises(SystemExit, match="cannot retain a freeze time"):
        fz_tool.validate(register)

    rendered = fz_tool.render(live_register(), fz_tool.validate(live_register()))
    active_table = rendered.split("## Superseded records outside the ladder", 1)[0]
    assert "| FZ-06 |" not in active_table
    assert "These identifiers preserve attested historical bytes" in rendered
    assert "| FZ-06 |" in rendered


def test_fz02_hash_is_bound_to_the_live_receipt():
    register = live_register()
    register["rows"][1]["content_sha256"] = "0" * 64
    with pytest.raises(SystemExit, match="does not equal the live"):
        fz_tool.validate(register)


def test_frozen_row_requires_attestation_fields():
    register = live_register()
    register["rows"][0]["attestation"] = None
    with pytest.raises(SystemExit, match="requires attestation"):
        fz_tool.validate(register)


def test_frozen_row_requires_parseable_nonfuture_utc():
    register = live_register()
    register["rows"][1]["frozen_utc"] = "not-a-time"
    with pytest.raises(SystemExit, match="ISO-8601 UTC"):
        fz_tool.validate(register)

    register = live_register()
    register["rows"][1]["frozen_utc"] = "2999-01-01T00:00:00Z"
    with pytest.raises(SystemExit, match="cannot be in the future"):
        fz_tool.validate(register)


def test_fz02_time_must_equal_corrected_custody_commit_time():
    register = live_register()
    register["external_custody_contracts"]["FZ-02"][
        "custody_commit_utc"
    ] = "2026-07-26T06:41:54Z"
    with pytest.raises(SystemExit, match="must equal the corrected"):
        fz_tool.validate(register)


def test_fz02_source_and_custody_commits_are_pinned():
    register = live_register()
    register["external_custody_contracts"]["FZ-02"]["source_commit"] = "0" * 40
    with pytest.raises(SystemExit, match="source commit must contain"):
        fz_tool.validate(register)

    register = live_register()
    register["external_custody_contracts"]["FZ-02"]["custody_commit"] = "0" * 40
    with pytest.raises(SystemExit, match="custody commit must equal"):
        fz_tool.validate(register)


def test_in_repo_custody_artifact_hashes_are_verified():
    register = live_register()
    register["external_custody_contracts"]["FZ-02"]["in_repo_artifact_sha256"][
        "Lean/Screen/A5AngularMultiplets.lean"
    ] = ("0" * 64)
    with pytest.raises(SystemExit, match="in-repo custody artifact hash mismatch"):
        fz_tool.validate(register)


def test_pending_row_requires_a_live_owning_issue():
    register = live_register()
    register["rows"][2]["owning_issue"] = 599
    with pytest.raises(SystemExit, match="not open in the snapshot"):
        fz_tool.validate(register)


def test_pending_row_requires_a_kill_band():
    register = live_register()
    register["rows"][3]["kill_band"] = ""
    with pytest.raises(SystemExit, match="kill_band must be nonempty"):
        fz_tool.validate(register)


def test_issue506_result_is_retrospective_and_bound_to_live_verdict():
    register = live_register()
    result = register["retrospective_results"][0]
    assert result["status"] == "retrospective_not_evaluable"
    assert result["former_ladder_reservation"] not in {
        row["id"] for row in register["rows"]
    }
    fz_tool.validate(register)

    register["retrospective_results"][0]["payload_sha256"] = "0" * 64
    with pytest.raises(SystemExit, match="retrospective payload hash"):
        fz_tool.validate(register)


def test_retrospective_reservation_cannot_be_inserted_into_ladder():
    register = live_register()
    register["rows"].insert(
        3,
        {
            "attestation": None,
            "comparison_protocol": "retrospective",
            "content": "not a prospective freeze",
            "content_sha256": None,
            "custody": None,
            "frozen_utc": None,
            "id": "FZ-04",
            "kill_band": "none",
            "milestone": "C1",
            "owning_issue": 506,
            "status": "registered_pending_freeze",
        },
    )
    with pytest.raises(SystemExit, match="must not appear"):
        fz_tool.validate(register)


def write_verdict(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload), encoding="utf-8")


def refresh_verdict_digest(payload: dict) -> None:
    canonical = {
        key: value for key, value in payload.items() if key != "verdict_sha256"
    }
    payload["verdict_sha256"] = "sha256:" + fz_tool.sha256_bytes(
        fz_tool.canonical_json_bytes(canonical)
    )


def test_issue506_self_reported_digest_is_not_trusted(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    payload = json.loads(fz_tool.FZ04_VERDICT_PATH.read_text(encoding="utf-8"))
    payload["class_matrix"]["tabulated_dispersive"]["claim_boundary"] += " altered"
    path = tmp_path / "verdict.json"
    write_verdict(path, payload)
    monkeypatch.setattr(fz_tool, "FZ04_VERDICT_PATH", path)
    with pytest.raises(SystemExit, match="self-digest"):
        fz_tool.validate(live_register())


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("prospective_freeze", True, "scope differs"),
        ("physical_alpha_prediction_emitted", True, "scope differs"),
    ],
)
def test_issue506_scope_mutations_fail_even_with_refreshed_self_digest(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    field: str,
    value: bool,
    message: str,
):
    payload = json.loads(fz_tool.FZ04_VERDICT_PATH.read_text(encoding="utf-8"))
    payload["scope"][field] = value
    refresh_verdict_digest(payload)
    path = tmp_path / "verdict.json"
    write_verdict(path, payload)
    monkeypatch.setattr(fz_tool, "FZ04_VERDICT_PATH", path)
    with pytest.raises(SystemExit, match=message):
        fz_tool.validate(live_register())


def test_issue506_claim_mutation_fails_even_with_refreshed_self_digest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    payload = json.loads(fz_tool.FZ04_VERDICT_PATH.read_text(encoding="utf-8"))
    payload["claim_boundary"] = "A stronger claim."
    refresh_verdict_digest(payload)
    path = tmp_path / "verdict.json"
    write_verdict(path, payload)
    monkeypatch.setattr(fz_tool, "FZ04_VERDICT_PATH", path)
    with pytest.raises(SystemExit, match="claim boundary"):
        fz_tool.validate(live_register())


def test_issue506_payload_mutation_fails_with_both_digests_refreshed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    canonical = json.loads(fz_tool.FZ04_VERDICT_PATH.read_text(encoding="utf-8"))
    payload = json.loads(json.dumps(canonical))
    payload["class_matrix"]["tabulated_dispersive"]["claim_boundary"] += " altered"
    refresh_verdict_digest(payload)
    path = tmp_path / "verdict.json"
    write_verdict(path, payload)

    register = live_register()
    register["retrospective_results"][0]["payload_sha256"] = payload[
        "verdict_sha256"
    ].removeprefix("sha256:")
    monkeypatch.setattr(fz_tool, "FZ04_VERDICT_PATH", path)
    monkeypatch.setattr(
        fz_tool, "rebuild_issue506_verdict", lambda: canonical
    )
    with pytest.raises(SystemExit, match="canonical producer replay"):
        fz_tool.validate(register)


def test_attestation_status_cannot_promote_calendar_pending_to_bitcoin():
    register = live_register()
    register["external_custody_contracts"]["FZ-02"][
        "attestation_state"
    ] = "bitcoin_attested"
    with pytest.raises(SystemExit, match="Bitcoin custody requires"):
        fz_tool.validate(register)


def test_external_custody_is_verified_when_present_or_explicitly_classified():
    register = live_register()
    fz_tool.validate(register)
    result = fz_tool.verify_external_custody(register)
    assert result["state"] in {"verified", "external_custody_not_present"}
    if result["state"] == "verified":
        assert result["contracts"]["FZ-01"] == {
            "verification": "verified",
            "attestation_state": "bitcoin_attested",
        }
        assert result["contracts"]["FZ-02"] == {
            "verification": "verified",
            "attestation_state": "calendar_pending",
        }
    else:
        assert {row["verification"] for row in result["contracts"].values()} == {
            "external_custody_not_present"
        }


def test_isolated_clone_reports_external_custody_not_present(tmp_path: Path):
    register = live_register()
    result = fz_tool.verify_external_custody(register, tmp_path / "isolated")
    assert result["state"] == "external_custody_not_present"
    assert {row["verification"] for row in result["contracts"].values()} == {
        "external_custody_not_present"
    }


def copy_external_custody(tmp_path: Path) -> Path:
    custody_root = tmp_path / "oph-meta"
    source = fz_tool.DEFAULT_CUSTODY_ROOT / "falsification"
    if not source.is_dir():
        pytest.skip("sibling oph-meta custody checkout is not present")
    shutil.copytree(source, custody_root / "falsification")
    return custody_root


def test_external_custody_rejects_artifact_tampering(tmp_path: Path):
    register = live_register()
    custody_root = copy_external_custody(tmp_path)
    target = (
        custody_root
        / "falsification"
        / "frozen_targets"
        / "fz02_2026-07-26"
        / "A5AngularMultiplets_pinned.lean"
    )
    target.write_bytes(target.read_bytes() + b"\n")
    with pytest.raises(SystemExit, match="custody artifact hash mismatch"):
        fz_tool.verify_external_custody(register, custody_root)


def test_external_custody_rejects_detached_timestamp_digest_tampering(
    tmp_path: Path,
):
    register = live_register()
    custody_root = copy_external_custody(tmp_path)
    stamp = (
        custody_root
        / "falsification"
        / "frozen_targets"
        / "fz02_2026-07-26"
        / "a5_angular_multiplet_reference.receipt.json.ots"
    )
    proof = bytearray(stamp.read_bytes())
    digest_offset = len(fz_tool.OTS_DETACHED_HEADER) + len(fz_tool.OTS_SHA256_TAG)
    proof[digest_offset] ^= 1
    stamp.write_bytes(proof)
    with pytest.raises(SystemExit, match="OpenTimestamps digest mismatch"):
        fz_tool.verify_external_custody(register, custody_root)
