"""Fail-closed tests for the #607 frozen-prediction ladder register."""

import copy
import importlib.util
import json
import shutil
import subprocess
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


def live_owner_successors() -> dict:
    return json.loads(fz_tool.OWNER_SUCCESSOR_PATH.read_text(encoding="utf-8"))


def write_owner_successors(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload), encoding="utf-8")


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
        "FZ-11",
        "FZ-12",
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


def test_fz11_is_a_frozen_unarmed_branch_prediction():
    register = live_register()
    fz11 = next(row for row in register["rows"] if row["id"] == "FZ-11")
    assert fz11["status"] == "frozen_stamped_upgrade_pending"
    assert fz11["owning_issue"] == 655
    assert "B6/C4^2=32/315" in fz11["content"]
    assert "Physical comparison is unarmed" in fz11["comparison_protocol"]
    assert "intrinsic C4 is positive" in fz11["kill_band"]
    assert "anisotropic coefficient at angular rank one through five" in fz11[
        "kill_band"
    ]
    assert "real, reciprocal, finite-range cosine kinetic branch" in fz11["content"]
    assert "no independent isotropic k^4 or k^6 term" in fz11["content"]
    assert "three-dimensional quotient SO(3)/A5" in fz11["content"]
    assert "not one scalar parameter" in fz11["content"]
    assert "FAIL rejects the primitive twelve-port" in fz11["kill_band"]
    fz_tool.validate(register)

    fz11["kill_band"] = fz11["kill_band"].replace(
        "intrinsic C4 is positive", "intrinsic C4 is negative"
    )
    with pytest.raises(SystemExit, match="kill_band must equal"):
        fz_tool.validate(register)


@pytest.mark.parametrize(
    "field", ["content", "comparison_protocol", "kill_band"]
)
def test_fz11_registry_contract_rejects_arbitrary_overclaim_prefix(field: str):
    register = live_register()
    fz11 = next(row for row in register["rows"] if row["id"] == "FZ-11")
    fz11[field] = "UNCONDITIONAL OPH-WIDE PREDICTION. " + fz11[field]
    with pytest.raises(SystemExit, match=rf"FZ-11 {field} must equal"):
        fz_tool.validate(register)


def test_fz11_hash_and_custody_commits_are_bound():
    register = live_register()
    fz11 = next(row for row in register["rows"] if row["id"] == "FZ-11")
    fz11["content_sha256"] = "0" * 64
    with pytest.raises(SystemExit, match="content hash must equal"):
        fz_tool.validate(register)

    register = live_register()
    register["external_custody_contracts"]["FZ-11"]["source_commit"] = "0" * 40
    with pytest.raises(SystemExit, match="source commit must contain"):
        fz_tool.validate(register)

    register = live_register()
    register["external_custody_contracts"]["FZ-11"]["custody_commit"] = "0" * 40
    with pytest.raises(SystemExit, match="custody commit must contain"):
        fz_tool.validate(register)

    register = live_register()
    register["external_custody_contracts"]["FZ-11"][
        "lean_repair_source_commit"
    ] = "0" * 40
    with pytest.raises(SystemExit, match="Lean repair source commit"):
        fz_tool.validate(register)

    register = live_register()
    register["external_custody_contracts"]["FZ-11"][
        "lean_repair_custody_commit"
    ] = "0" * 40
    with pytest.raises(SystemExit, match="Lean repair custody commit"):
        fz_tool.validate(register)

    register = live_register()
    register["external_custody_contracts"]["FZ-11"][
        "decision_rule_custody_commit"
    ] = "0" * 40
    with pytest.raises(SystemExit, match="decision-rule custody commit"):
        fz_tool.validate(register)


def test_fz12_is_a_frozen_source_native_unarmed_branch_prediction():
    register = live_register()
    fz12 = next(row for row in register["rows"] if row["id"] == "FZ-12")
    assert fz12["status"] == "frozen_stamped_upgrade_pending"
    assert fz12["owning_issue"] == 666
    assert "exact D6 image" in fz12["content"]
    assert "B6/C4^2=-2/63" in fz12["content"]
    assert "without proving the physical field action" in fz12["content"]
    assert "Physical comparison is ineligible and unarmed" in fz12[
        "comparison_protocol"
    ]
    assert "same-action, same-sector source theorem" in fz12["kill_band"]
    assert "including every a >= a_min" in fz12["kill_band"]
    assert "forced and exclusive physical bridge theorem" in fz12["kill_band"]
    fz_tool.validate(register)


def test_fz12_hash_source_custody_time_and_ownership_are_bound():
    register = live_register()
    fz12 = next(row for row in register["rows"] if row["id"] == "FZ-12")
    fz12["content_sha256"] = "0" * 64
    with pytest.raises(SystemExit, match="content hash must equal"):
        fz_tool.validate(register)

    register = live_register()
    register["external_custody_contracts"]["FZ-12"]["source_commit"] = "0" * 40
    with pytest.raises(SystemExit, match="source commit must contain"):
        fz_tool.validate(register)

    register = live_register()
    register["external_custody_contracts"]["FZ-12"]["custody_commit"] = "0" * 40
    with pytest.raises(SystemExit, match="custody commit must contain"):
        fz_tool.validate(register)

    register = live_register()
    register["external_custody_contracts"]["FZ-12"]["frozen_utc"] = (
        "2026-08-02T11:52:28Z"
    )
    with pytest.raises(SystemExit, match="fixed freeze time"):
        fz_tool.validate(register)

    register = live_register()
    fz12 = next(row for row in register["rows"] if row["id"] == "FZ-12")
    fz12["owning_issue"] = 655
    with pytest.raises(SystemExit, match="ownership, status, or custody"):
        fz_tool.validate(register)


def test_fz12_source_commit_resolves_all_exact_historical_blobs():
    contract = live_register()["external_custody_contracts"]["FZ-12"]
    result = fz_tool.verify_fz12_source_history(contract)
    assert result["state"] == "verified"
    assert result["source_commit"] == fz_tool.FZ12_SOURCE_COMMIT
    assert set(result["blobs"]) == fz_tool.FZ12_IN_REPO_ARTIFACTS
    assert set(result["historical_parent_blobs"]) == {
        "code/a5_fingerprint/runtime/a5_multipole_fixed_point_receipt.json",
        fz_tool.FZ12_CARRIER_LEAN_REL,
        fz_tool.FZ12_MOMENT_LEAN_REL,
        fz_tool.FZ12_RAY_LEAN_REL,
    }
    assert result["historical_parent_blobs"][
        "code/a5_fingerprint/runtime/a5_multipole_fixed_point_receipt.json"
    ] == "d96b80c71a64f48bb7a2a7b2592bad5e122ce9c5f3ca3b4636af29693c426ecd"


def write_fz12_receipt(path: Path, receipt: dict) -> str:
    body = {key: value for key, value in receipt.items() if key != "receipt_sha256"}
    receipt["receipt_sha256"] = "sha256:" + fz_tool.sha256_bytes(
        fz_tool.canonical_json_bytes(body)
    )
    payload = fz_tool.canonical_json_bytes(receipt)
    path.write_bytes(payload)
    return fz_tool.sha256_bytes(payload)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("comparison_permitted", "target exposure or comparison boundary"),
        ("gate_promoted", "promotion gates drifted"),
        ("supersedes_fz11", "must remain distinct"),
        ("lower_bound_removed", "lower-bound gate"),
    ],
)
def test_fz12_receipt_boundary_mutations_fail_with_refreshed_digests(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
    message: str,
):
    receipt = json.loads(fz_tool.FZ12_RECEIPT_PATH.read_text(encoding="utf-8"))
    if mutation == "comparison_permitted":
        receipt["exposure_and_custody_boundary"]["comparison_permitted"] = True
    elif mutation == "gate_promoted":
        receipt["promotion_gates"]["gates"][0]["status"] = "ATTAINED"
    elif mutation == "supersedes_fz11":
        receipt["fz11_separation"]["supersedes_fz11"] = True
    elif mutation == "lower_bound_removed":
        receipt["physical_premises"]["positive_scale_lower_bound"] = (
            "mathematical a > 0"
        )
    else:
        raise AssertionError("unhandled mutation")
    receipt_path = tmp_path / "fz12.json"
    raw_digest = write_fz12_receipt(receipt_path, receipt)
    monkeypatch.setattr(fz_tool, "FZ12_RECEIPT_PATH", receipt_path)
    register = live_register()
    rows_by_id = {row["id"]: row for row in register["rows"]}
    rows_by_id["FZ-12"]["content_sha256"] = raw_digest
    contract = register["external_custody_contracts"]["FZ-12"]
    contract["in_repo_artifact_sha256"][fz_tool.FZ12_RECEIPT_REL] = raw_digest
    with pytest.raises(SystemExit, match=message):
        fz_tool.validate_fz12_prediction(rows_by_id, contract)


def test_fz11_source_commits_resolve_exact_historical_blobs_and_parentage():
    contract = live_register()["external_custody_contracts"]["FZ-11"]
    result = fz_tool.verify_fz11_source_history(contract)
    assert result["state"] == "verified"
    assert result["source_commit"] == fz_tool.FZ11_SOURCE_COMMIT
    assert result["repair_commit"] == fz_tool.FZ11_LEAN_REPAIR_SOURCE_COMMIT
    assert result["repair_is_direct_child"] is True
    assert result["original_blobs"][fz_tool.FZ11_RECEIPT_REL] == result[
        "repaired_blobs"
    ][fz_tool.FZ11_RECEIPT_REL]
    assert result["original_blobs"][fz_tool.FZ11_LEAN_REL] != result[
        "repaired_blobs"
    ][fz_tool.FZ11_LEAN_REL]


@pytest.mark.parametrize(
    ("mapping", "path", "message"),
    [
        (
            "original_in_repo_artifact_sha256",
            fz_tool.FZ11_RECEIPT_REL,
            "historical blob hash mismatch",
        ),
        (
            "original_in_repo_artifact_sha256",
            fz_tool.FZ11_LEAN_REL,
            "historical blob hash mismatch",
        ),
        (
            "in_repo_artifact_sha256",
            fz_tool.FZ11_LEAN_REL,
            "historical blob hash mismatch",
        ),
    ],
)
def test_fz11_source_history_rejects_historical_blob_hash_mutations(
    mapping: str, path: str, message: str
):
    contract = copy.deepcopy(
        live_register()["external_custody_contracts"]["FZ-11"]
    )
    contract[mapping][path] = "0" * 64
    with pytest.raises(SystemExit, match=message):
        fz_tool.verify_fz11_source_history(contract)


def test_fz11_source_repair_requires_a_direct_single_parent():
    with pytest.raises(SystemExit, match="direct single-parent child"):
        fz_tool.verify_direct_parent(
            fz_tool.ROOT,
            fz_tool.FZ11_LEAN_REPAIR_SOURCE_COMMIT,
            "0" * 40,
            "mutated source ancestry",
        )


def test_exported_tree_classifies_absent_source_history(tmp_path: Path):
    contract = live_register()["external_custody_contracts"]["FZ-11"]
    result = fz_tool.verify_fz11_source_history(contract, tmp_path)
    assert result["state"] == "git_history_not_present"


def test_present_but_invalid_git_metadata_fails_closed(tmp_path: Path):
    (tmp_path / ".git").write_text("not a gitdir", encoding="utf-8")
    with pytest.raises(SystemExit, match="checkout cannot be resolved"):
        fz_tool.git_checkout_root(tmp_path)


def test_fz02_frame_lock_is_ineligible_pending_issue_643():
    register = live_register()
    fz02 = next(row for row in register["rows"] if row["id"] == "FZ-02")
    assert "frame-lock clause is not established" in fz02["content"]
    assert "ineligible pending issue #643" in fz02["content"]
    assert "Status correction 2026-07-30" in fz02["content"]
    assert "retired from the scientific target" in fz02["content"]
    assert "no frame-lock verdict may be issued" in fz02["comparison_protocol"]
    assert "FZ02-R03a" in fz02["kill_band"]
    assert "FZ02-R03b (INELIGIBLE; #643)" in fz02["kill_band"]
    assert "FZ02-R03b correction 2026-07-30" in fz02["kill_band"]
    fz_tool.validate(register)

    register = live_register()
    fz02 = next(row for row in register["rows"] if row["id"] == "FZ-02")
    fz02["comparison_protocol"] = (
        "The registered frame-lock clause is READY for comparison."
    )
    with pytest.raises(SystemExit, match="unsupported frame-lock clause"):
        fz_tool.validate(register)


def test_fz05_requires_physical_attachment_and_retrospective_exposure():
    register = live_register()
    fz05 = next(row for row in register["rows"] if row["id"] == "FZ-05")
    assert all(issue in fz05["content"] for issue in ("#729", "#738", "#736"))
    assert "retrospective" in fz05["comparison_protocol"].lower()
    assert all(issue in fz05["kill_band"] for issue in ("#729", "#738", "#736"))

    fz05["content"] = "A positive finite N is a cosmological prediction."
    with pytest.raises(SystemExit, match="FZ-05 must keep"):
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
    with pytest.raises(SystemExit, match="does not match FZ-03 owning_issue"):
        fz_tool.validate(register)


def test_historical_owners_are_preserved_and_v2_successors_are_live():
    register = live_register()
    rows = register["rows"]
    open_issues = {
        row["number"]
        for row in json.loads(fz_tool.SNAPSHOT_PATH.read_text(encoding="utf-8"))[
            "rows"
        ]
    }
    mapping = fz_tool.validate_owner_successors(rows, open_issues)
    assert mapping == {
        "FZ-03": (736,),
        "FZ-05": (738,),
        "FZ-10": (736,),
        "FZ-11": (733, 736),
        "FZ-12": (733, 736),
    }
    assert {
        row["id"]: row["owning_issue"]
        for row in rows
        if row["id"] in mapping
    } == {
        "FZ-03": 508,
        "FZ-05": 639,
        "FZ-10": 546,
        "FZ-11": 655,
        "FZ-12": 666,
    }
    rendered = fz_tool.render(register, fz_tool.validate(register))
    assert "historical [#508]" in rendered
    assert "active V2 [#736]" in rendered
    assert "historical [#639]" in rendered
    assert "active V2 [#738]" in rendered


def test_owner_successor_map_fails_closed_when_a_historical_row_is_missing(
    tmp_path: Path,
):
    payload = live_owner_successors()
    payload["mappings"] = [
        mapping for mapping in payload["mappings"] if mapping["row_id"] != "FZ-03"
    ]
    path = tmp_path / "successors.json"
    write_owner_successors(path, payload)
    rows = live_register()["rows"]
    open_issues = {
        row["number"]
        for row in json.loads(fz_tool.SNAPSHOT_PATH.read_text(encoding="utf-8"))[
            "rows"
        ]
    }
    with pytest.raises(SystemExit, match=r"missing=\['FZ-03'\]"):
        fz_tool.validate_owner_successors(rows, open_issues, path)


def test_owner_successor_map_fails_closed_on_unknown_active_issue(tmp_path: Path):
    payload = live_owner_successors()
    payload["mappings"][0]["active_successor_issues"] = [999999]
    path = tmp_path / "successors.json"
    write_owner_successors(path, payload)
    rows = live_register()["rows"]
    open_issues = {
        row["number"]
        for row in json.loads(fz_tool.SNAPSHOT_PATH.read_text(encoding="utf-8"))[
            "rows"
        ]
    }
    with pytest.raises(SystemExit, match=r"snapshot: \[999999\]"):
        fz_tool.validate_owner_successors(rows, open_issues, path)


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


def lean_axiom_output(*, extra: str = "") -> str:
    reports = []
    for name in fz_tool.FZ11_AXIOM_DECLARATIONS:
        reports.append(
            f"'{name}' depends on axioms: [propext, Classical.choice, Quot.sound]"
        )
    return "\n".join(reports) + extra


def test_fz11_axiom_report_parser_accepts_exact_five_standard_reports():
    reports = fz_tool.parse_fz11_axiom_reports(lean_axiom_output())
    assert tuple(reports) == fz_tool.FZ11_AXIOM_DECLARATIONS
    assert all(
        set(axioms) == fz_tool.FZ11_ALLOWED_AXIOMS
        for axioms in reports.values()
    )


@pytest.mark.parametrize(
    ("output", "message"),
    [
        (
            "\n".join(lean_axiom_output().splitlines()[:-1]),
            "exactly five",
        ),
        (
            lean_axiom_output(
                extra="\n'extra.theorem' depends on axioms: [propext]"
            ),
            "exactly five",
        ),
        (lean_axiom_output(extra="\nsorryAx"), "reported sorryAx"),
        (
            lean_axiom_output().replace(
                "[propext, Classical.choice, Quot.sound]",
                "[propext, Classical.choice, Quot.sound, Classical.sorryAx]",
                1,
            ),
            "reported sorryAx",
        ),
        (
            lean_axiom_output().replace("Quot.sound]", "Classical.choice]", 1),
            "unexpected axioms",
        ),
    ],
)
def test_fz11_axiom_report_parser_rejects_mutations(output: str, message: str):
    with pytest.raises(SystemExit, match=message):
        fz_tool.parse_fz11_axiom_reports(output)


def test_explicit_fz11_lean_replay_fails_closed_without_lake(
    monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setattr(fz_tool.shutil, "which", lambda executable: None)
    with pytest.raises(SystemExit, match="never skips"):
        fz_tool.verify_fz11_lean_replay()


@pytest.mark.skipif(
    shutil.which("lake") is None,
    reason="the dedicated replay is mandatory in .github/workflows/lean-ci.yml",
)
def test_fz11_lean_replay_when_toolchain_is_available():
    result = fz_tool.verify_fz11_lean_replay()
    assert result["lean_exit_code"] == 0
    assert result["report_count"] == 5
    assert result["sorry_ax_present"] is False


def test_official_ots_cli_is_used_when_available(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    stamp = tmp_path / "stamp.ots"
    stamp.write_bytes(b"proof")
    executable = tmp_path / "ots"
    executable.write_text("stub", encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setattr(fz_tool, "_OTS_CLI_USABLE", None)
    monkeypatch.setattr(fz_tool, "_OTS_CLI_PATH", None)
    monkeypatch.setattr(fz_tool, "_OTS_OFFICIAL_PARSED_DIGESTS", set())
    monkeypatch.setattr(fz_tool.shutil, "which", lambda name: str(executable))
    monkeypatch.setattr(
        fz_tool.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args=args[0], returncode=0, stdout="parsed", stderr=""
        ),
    )
    assert fz_tool.verify_ots_with_official_cli(stamp) == "verified_by_ots_info"


def test_official_ots_cli_rejection_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    stamp = tmp_path / "stamp.ots"
    stamp.write_bytes(b"malformed")
    executable = tmp_path / "ots"
    executable.write_text("stub", encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setattr(fz_tool, "_OTS_CLI_USABLE", None)
    monkeypatch.setattr(fz_tool, "_OTS_CLI_PATH", None)
    monkeypatch.setattr(fz_tool, "_OTS_OFFICIAL_PARSED_DIGESTS", set())
    monkeypatch.setattr(fz_tool.shutil, "which", lambda name: str(executable))
    monkeypatch.setattr(
        fz_tool.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args=args[0], returncode=2, stdout="", stderr="invalid timestamp"
        ),
    )
    with pytest.raises(SystemExit, match="official OpenTimestamps parser rejected"):
        fz_tool.verify_ots_with_official_cli(stamp)


def test_structural_ots_prefix_cannot_produce_verified_custody(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    expected = "ab" * 32
    stamp = tmp_path / "fabricated.ots"
    stamp.write_bytes(
        fz_tool.OTS_DETACHED_HEADER
        + fz_tool.OTS_SHA256_TAG
        + bytes.fromhex(expected)
        + fz_tool.OTS_PENDING_ATTESTATION_TAG
    )
    # This is enough to fool the bounded prefix extractor, which is precisely
    # why the external-custody verdict requires an official parser.
    assert fz_tool.ots_binding(stamp)["file_sha256"] == expected
    monkeypatch.setattr(fz_tool, "_OTS_CLI_USABLE", None)
    monkeypatch.setattr(fz_tool, "_OTS_CLI_PATH", None)
    monkeypatch.setattr(fz_tool, "_OTS_OFFICIAL_PARSED_DIGESTS", set())
    monkeypatch.setattr(fz_tool.shutil, "which", lambda name: None)
    monkeypatch.setattr(fz_tool.os, "access", lambda *args: False)
    with pytest.raises(SystemExit, match="official OpenTimestamps parser is required"):
        fz_tool.verify_ots_binding(stamp, expected, "calendar_pending")


def test_fabricated_structural_stamp_cannot_verify_external_custody(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    register = live_register()
    custody_root = copy_external_custody(tmp_path)
    contract = register["external_custody_contracts"]["FZ-01"]
    filename, expected = next(iter(contract["artifact_sha256"].items()))
    stamp = custody_root / contract["custody_path"] / (filename + ".ots")
    stamp.write_bytes(
        fz_tool.OTS_DETACHED_HEADER
        + fz_tool.OTS_SHA256_TAG
        + bytes.fromhex(expected)
        + fz_tool.OTS_BITCOIN_ATTESTATION_TAG
    )
    assert fz_tool.ots_binding(stamp)["file_sha256"] == expected
    monkeypatch.setattr(
        fz_tool,
        "verify_ots_with_official_cli",
        lambda path: (
            "official_cli_not_available"
            if path == stamp
            else "verified_by_ots_info"
        ),
    )
    with pytest.raises(SystemExit, match="official OpenTimestamps parser is required"):
        fz_tool.verify_external_custody(register, custody_root)


def test_external_custody_is_verified_when_present_or_explicitly_classified():
    register = live_register()
    fz_tool.validate(register)
    result = fz_tool.verify_external_custody(register)
    assert result["state"] in {"verified", "external_custody_not_present"}
    if result["state"] == "verified":
        assert result["fz11_custody_git_history"]["state"] in {
            "verified",
            "git_history_not_present",
        }
        assert result["fz12_custody_git_history"]["state"] in {
            "verified",
            "git_history_not_present",
        }
        assert result["contracts"]["FZ-01"] == {
            "verification": "verified",
            "attestation_state": "bitcoin_attested",
        }
        assert result["contracts"]["FZ-02"] == {
            "verification": "verified",
            "attestation_state": "calendar_pending",
        }
        assert result["contracts"]["FZ-11"] == {
            "verification": "verified",
            "attestation_state": "calendar_pending",
        }
        assert result["contracts"]["FZ-12"] == {
            "verification": "verified",
            "attestation_state": "calendar_pending",
        }
    else:
        assert result["fz11_custody_git_history"]["state"] == (
            "external_custody_not_present"
        )
        assert result["fz12_custody_git_history"]["state"] == (
            "external_custody_not_present"
        )
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


@pytest.mark.parametrize(
    "field",
    ["prediction_bytes_changed", "source_receipt_bytes_changed", "comparison_data_read"],
)
def test_external_custody_rejects_fz12_clarification_boundary_drift(
    tmp_path: Path, field: str
):
    register = live_register()
    custody_root = copy_external_custody(tmp_path)
    contract = register["external_custody_contracts"]["FZ-12"]
    manifest = (
        custody_root
        / contract["custody_path"]
        / contract["decision_rule_manifest"]
    )
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    payload[field] = True
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    contract["decision_rule_manifest_sha256"] = fz_tool.sha256_file(manifest)
    with pytest.raises(SystemExit, match="clarification manifest bindings drifted"):
        fz_tool.verify_external_custody(register, custody_root)


def test_external_custody_rejects_weakened_fz12_null_quantifier(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    register = live_register()
    custody_root = copy_external_custody(tmp_path)
    contract = register["external_custody_contracts"]["FZ-12"]
    directory = custody_root / contract["custody_path"]
    rule = directory / fz_tool.FZ12_DECISION_RULE_JSON
    payload = json.loads(rule.read_text(encoding="utf-8"))
    payload["null_rule"] = "Every null is inconclusive."
    rule.write_text(json.dumps(payload), encoding="utf-8")
    rule_hash = fz_tool.sha256_file(rule)
    contract["decision_rule_artifact_sha256"][rule.name] = rule_hash
    manifest = directory / contract["decision_rule_manifest"]
    manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
    manifest_payload["append_only_artifacts"][rule.name] = rule_hash
    manifest.write_text(json.dumps(manifest_payload), encoding="utf-8")
    contract["decision_rule_manifest_sha256"] = fz_tool.sha256_file(manifest)
    monkeypatch.setattr(fz_tool, "verify_ots_binding", lambda *args: {})
    with pytest.raises(SystemExit, match="clarified decision-rule bindings drifted"):
        fz_tool.verify_external_custody(register, custody_root)


def test_coordinated_fz11_root_custody_history_is_strict_when_available():
    contract = live_register()["external_custody_contracts"]["FZ-11"]
    directory = fz_tool.DEFAULT_CUSTODY_ROOT / contract["custody_path"]
    if not directory.is_dir():
        pytest.skip("coordinated oph-meta custody checkout is not present")
    result = fz_tool.verify_fz11_custody_history(
        contract, fz_tool.DEFAULT_CUSTODY_ROOT, directory
    )
    if fz_tool.git_checkout_root(fz_tool.DEFAULT_CUSTODY_ROOT) is None:
        assert result["state"] == "git_history_not_present"
    else:
        assert result["state"] == "verified"
        assert result["repair_is_direct_child"] is True
        assert result["original_blob_count"] == 3
        assert result["repair_blob_count"] > result["original_blob_count"]
        assert result["decision_rule_blob_count"] > result["repair_blob_count"]
        assert result["original_proof_count"] == 3
        assert result["repair_proof_count"] == 6
        assert result["decision_rule_proof_count"] == 9
        assert result["decision_rule_is_direct_child"] is True
        with pytest.raises(SystemExit, match="direct single-parent child"):
            fz_tool.verify_direct_parent(
                fz_tool.DEFAULT_CUSTODY_ROOT,
                fz_tool.FZ11_LEAN_REPAIR_CUSTODY_COMMIT,
                "0" * 40,
                "mutated root custody ancestry",
            )
        with pytest.raises(SystemExit, match="direct single-parent child"):
            fz_tool.verify_direct_parent(
                fz_tool.DEFAULT_CUSTODY_ROOT,
                fz_tool.FZ11_DECISION_RULE_CUSTODY_COMMIT,
                "0" * 40,
                "mutated decision-rule ancestry",
            )


def test_coordinated_fz12_root_custody_history_is_strict_when_available():
    contract = live_register()["external_custody_contracts"]["FZ-12"]
    directory = fz_tool.DEFAULT_CUSTODY_ROOT / contract["custody_path"]
    if not directory.is_dir():
        pytest.skip("coordinated oph-meta custody checkout is not present")
    result = fz_tool.verify_fz12_custody_history(
        contract, fz_tool.DEFAULT_CUSTODY_ROOT
    )
    if fz_tool.git_checkout_root(fz_tool.DEFAULT_CUSTODY_ROOT) is None:
        assert result["state"] == "git_history_not_present"
    else:
        assert result["state"] == "verified"
        assert result["custody_commit"] == fz_tool.FZ12_CUSTODY_COMMIT
        assert result["decision_rule_commit"] == (
            fz_tool.FZ12_DECISION_RULE_CUSTODY_COMMIT
        )
        assert result["original_blob_count"] == 3
        assert result["decision_rule_blob_count"] == 6
        assert result["original_proof_count"] == 3
        assert result["decision_rule_proof_count"] == 6
        assert result["decision_rule_is_direct_child"] is True
        with pytest.raises(SystemExit, match="direct single-parent child"):
            fz_tool.verify_direct_parent(
                fz_tool.DEFAULT_CUSTODY_ROOT,
                fz_tool.FZ12_DECISION_RULE_CUSTODY_COMMIT,
                "0" * 40,
                "mutated FZ-12 decision-rule ancestry",
            )


def test_historical_custody_does_not_pin_future_live_ots_bytes(tmp_path: Path):
    contract = live_register()["external_custody_contracts"]["FZ-11"]
    if fz_tool.git_checkout_root(fz_tool.DEFAULT_CUSTODY_ROOT) is None:
        pytest.skip("coordinated oph-meta Git history is not present")
    # Historical proofs come from each commit. No live custody proof is read
    # through this deliberately empty directory, so a future in-place OTS
    # upgrade cannot invalidate the original history chain.
    result = fz_tool.verify_fz11_custody_history(
        contract, fz_tool.DEFAULT_CUSTODY_ROOT, tmp_path
    )
    assert result["state"] == "verified"
    assert result["decision_rule_proof_count"] == 9


def test_copied_fz11_custody_has_no_false_git_history_claim(tmp_path: Path):
    register = live_register()
    custody_root = copy_external_custody(tmp_path)
    result = fz_tool.verify_external_custody(register, custody_root)
    assert result["state"] == "verified"
    assert result["fz11_custody_git_history"]["state"] == (
        "git_history_not_present"
    )
    assert result["fz12_custody_git_history"]["state"] == (
        "git_history_not_present"
    )


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


def test_external_custody_rejects_scientific_erratum_tampering(tmp_path: Path):
    register = live_register()
    custody_root = copy_external_custody(tmp_path)
    erratum = (
        custody_root
        / "falsification"
        / "frozen_targets"
        / "fz02_2026-07-26"
        / "SCIENTIFIC_ERRATUM_2026-07-29.md"
    )
    erratum.write_bytes(erratum.read_bytes() + b"\n")
    with pytest.raises(SystemExit, match="scientific erratum hash mismatch"):
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


def test_external_custody_rejects_fz11_snapshot_tampering(tmp_path: Path):
    register = live_register()
    custody_root = copy_external_custody(tmp_path)
    snapshot = (
        custody_root
        / "falsification"
        / "frozen_targets"
        / "fz11_2026-07-31"
        / "spin_six_primitive_port_prediction_frozen_2026-07-31.json"
    )
    snapshot.write_bytes(snapshot.read_bytes() + b"\n")
    with pytest.raises(SystemExit, match="custody artifact hash mismatch"):
        fz_tool.verify_external_custody(register, custody_root)


def test_external_custody_rejects_fz12_snapshot_boundary_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    register = live_register()
    custody_root = copy_external_custody(tmp_path)
    contract = register["external_custody_contracts"]["FZ-12"]
    directory = custody_root / contract["custody_path"]
    snapshot = directory / contract["prediction_file"]
    payload = json.loads(snapshot.read_text(encoding="utf-8"))
    payload["new_target_or_comparison_data_read"] = True
    snapshot.write_text(json.dumps(payload), encoding="utf-8")
    digest = fz_tool.sha256_file(snapshot)
    contract["artifact_sha256"][snapshot.name] = digest
    manifest = directory / contract["registration_manifest"]
    manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
    manifest_payload["artifacts"][snapshot.name] = digest
    manifest.write_text(json.dumps(manifest_payload), encoding="utf-8")
    contract["registration_manifest_sha256"] = fz_tool.sha256_file(manifest)
    monkeypatch.setattr(fz_tool, "verify_ots_binding", lambda *args: {})
    with pytest.raises(SystemExit, match="snapshot content drifted"):
        fz_tool.verify_external_custody(register, custody_root)


def test_external_custody_rejects_fz12_manifest_binding_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    register = live_register()
    custody_root = copy_external_custody(tmp_path)
    contract = register["external_custody_contracts"]["FZ-12"]
    manifest = (
        custody_root / contract["custody_path"] / contract["registration_manifest"]
    )
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    payload["source_commit"] = "0" * 40
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    contract["registration_manifest_sha256"] = fz_tool.sha256_file(manifest)
    monkeypatch.setattr(fz_tool, "verify_ots_binding", lambda *args: {})
    with pytest.raises(SystemExit, match="registration manifest bindings drifted"):
        fz_tool.verify_external_custody(register, custody_root)


def test_external_custody_rejects_fz11_manifest_binding_drift(tmp_path: Path):
    register = live_register()
    custody_root = copy_external_custody(tmp_path)
    manifest = (
        custody_root
        / "falsification"
        / "frozen_targets"
        / "fz11_2026-07-31"
        / "registration_manifest_2026-07-31.json"
    )
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    payload["source_commit"] = "0" * 40
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    register["external_custody_contracts"]["FZ-11"][
        "registration_manifest_sha256"
    ] = fz_tool.sha256_file(manifest)
    with pytest.raises(SystemExit, match="registration manifest bindings drifted"):
        fz_tool.verify_external_custody(register, custody_root)


def test_external_custody_rejects_fz11_lean_repair_tampering(tmp_path: Path):
    register = live_register()
    custody_root = copy_external_custody(tmp_path)
    repaired = (
        custody_root
        / "falsification"
        / "frozen_targets"
        / "fz11_2026-07-31"
        / "A5PrimitivePortPrediction_repaired.lean"
    )
    repaired.write_bytes(repaired.read_bytes() + b"\n")
    with pytest.raises(SystemExit, match="Lean repair artifact hash mismatch"):
        fz_tool.verify_external_custody(register, custody_root)


def test_external_custody_rejects_fz11_lean_repair_manifest_drift(
    tmp_path: Path,
):
    register = live_register()
    custody_root = copy_external_custody(tmp_path)
    manifest = (
        custody_root
        / "falsification"
        / "frozen_targets"
        / "fz11_2026-07-31"
        / "lean_repair_manifest_2026-07-31.json"
    )
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    payload["prediction_bytes_changed"] = True
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    register["external_custody_contracts"]["FZ-11"][
        "lean_repair_manifest_sha256"
    ] = fz_tool.sha256_file(manifest)
    with pytest.raises(SystemExit, match="Lean repair manifest bindings drifted"):
        fz_tool.verify_external_custody(register, custody_root)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("prediction_bytes_changed", True),
        ("comparison_data_read", True),
    ],
)
def test_external_custody_rejects_fz11_decision_manifest_boundary_drift(
    tmp_path: Path, field: str, value: bool
):
    register = live_register()
    custody_root = copy_external_custody(tmp_path)
    manifest = (
        custody_root
        / "falsification"
        / "frozen_targets"
        / "fz11_2026-07-31"
        / "decision_rule_erratum_manifest_2026-07-31.json"
    )
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    payload[field] = value
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    register["external_custody_contracts"]["FZ-11"][
        "decision_rule_manifest_sha256"
    ] = fz_tool.sha256_file(manifest)
    with pytest.raises(
        SystemExit, match="decision-rule correction manifest bindings drifted"
    ):
        fz_tool.verify_external_custody(register, custody_root)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("prediction_bytes_changed", True),
        ("comparison_data_read", True),
    ],
)
def test_external_custody_rejects_fz11_decision_rule_boundary_drift(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    field: str,
    value: bool,
):
    register = live_register()
    custody_root = copy_external_custody(tmp_path)
    decision = (
        custody_root
        / "falsification"
        / "frozen_targets"
        / "fz11_2026-07-31"
        / "fz11_decision_rule_v2_2026-07-31.json"
    )
    payload = json.loads(decision.read_text(encoding="utf-8"))
    payload[field] = value
    decision.write_text(json.dumps(payload), encoding="utf-8")
    digest = fz_tool.sha256_file(decision)
    register["external_custody_contracts"]["FZ-11"][
        "decision_rule_artifact_sha256"
    ][decision.name] = digest
    manifest = (
        decision.parent / "decision_rule_erratum_manifest_2026-07-31.json"
    )
    manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
    manifest_payload["append_only_artifacts"][decision.name] = digest
    manifest.write_text(json.dumps(manifest_payload), encoding="utf-8")
    register["external_custody_contracts"]["FZ-11"][
        "decision_rule_manifest_sha256"
    ] = fz_tool.sha256_file(manifest)
    monkeypatch.setattr(fz_tool, "verify_ots_binding", lambda *args: {})
    with pytest.raises(SystemExit, match="corrected decision-rule bindings drifted"):
        fz_tool.verify_external_custody(register, custody_root)


def test_external_custody_rejects_fz11_decision_artifact_tampering(
    tmp_path: Path,
):
    register = live_register()
    custody_root = copy_external_custody(tmp_path)
    erratum = (
        custody_root
        / "falsification"
        / "frozen_targets"
        / "fz11_2026-07-31"
        / "FZ11_DECISION_RULE_ERRATUM_2026-07-31.md"
    )
    erratum.write_bytes(erratum.read_bytes() + b"\n")
    with pytest.raises(SystemExit, match="decision-rule artifact hash mismatch"):
        fz_tool.verify_external_custody(register, custody_root)
