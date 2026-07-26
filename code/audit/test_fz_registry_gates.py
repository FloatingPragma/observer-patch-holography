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


def test_ladder_is_contiguous_fz01_through_fz09():
    rows = fz_tool.validate(live_register())
    assert [row["id"] for row in rows] == [f"FZ-{i:02d}" for i in range(1, 10)]


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
