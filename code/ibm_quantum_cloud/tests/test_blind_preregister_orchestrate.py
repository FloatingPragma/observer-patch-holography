from __future__ import annotations

import importlib.util
import json
import stat
import sys
from pathlib import Path

import pytest


pytest.importorskip("qiskit")

PROGRAMS = Path(__file__).resolve().parents[1] / "programs"
sys.path.insert(0, str(PROGRAMS))


def _load(name: str):
    path = PROGRAMS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


blind = _load("blind_preregister")
orchestrate = _load("blind_preregister_orchestrate")


def _catalog_builder(call_count: list[int]):
    def build():
        call_count.append(1)
        return blind.build_blinded_preregistration(
            test_seed="orchestrator-test",
            _test_scope={
                "cayley_models": ("z5",),
                "cayley_protocols": blind.CAYLEY_PROTOCOLS,
                "mh_spectra": ("z3_cyclic_control",),
                "calibration_codes": (0, 15),
            },
        )

    return build


def _analysis_builder(public, reveal):
    document = dict(reveal["sealed_payload"]["analysis_document"])
    document.pop("analysis_lock_sha256")
    document["hardened_orchestration_test"] = True
    document["catalog_precommitment_sha256"] = public[
        "catalog_precommitment_sha256"
    ]
    document["analysis_lock_sha256"] = blind.sha256_json(document)
    return document


def _analysis_validator(document):
    body = dict(document)
    claimed = body.pop("analysis_lock_sha256")
    assert blind.sha256_json(body) == claimed
    assert body["catalog_precommitment_sha256"]


@pytest.fixture(scope="module")
def orchestrated_bundle():
    calls: list[int] = []
    bundle = orchestrate.build_and_bind_production_bundle(
        catalog_builder=_catalog_builder(calls),
        analysis_builder=_analysis_builder,
        analysis_validator=_analysis_validator,
        _allow_test_mode=True,
    )
    assert calls == [1]
    return bundle


def test_two_phase_orchestration_builds_catalog_once_and_verifies_logical_hashes(
    orchestrated_bundle,
) -> None:
    public, reveal, analysis = orchestrated_bundle
    assert reveal["sealed_payload"]["analysis_document"] == analysis
    assert public["analysis_lock_sha256"] == analysis["analysis_lock_sha256"]
    assert blind.verify_bundle(public, reveal, analysis, rebuild_circuits=True)
    assert reveal["sealed_payload"]["mode"] == "deterministic_test_only"


def test_production_entry_rejects_deterministic_catalog() -> None:
    calls: list[int] = []
    with pytest.raises(orchestrate.SealOrchestrationError, match="production-random"):
        orchestrate.build_and_bind_production_bundle(
            catalog_builder=_catalog_builder(calls),
            analysis_builder=_analysis_builder,
            analysis_validator=_analysis_validator,
        )
    assert calls == [1]


def test_analysis_builder_cannot_switch_catalog() -> None:
    calls: list[int] = []

    def wrong_catalog(public, reveal):
        document = _analysis_builder(public, reveal)
        document["catalog_precommitment_sha256"] = "0" * 64
        body = dict(document)
        body.pop("analysis_lock_sha256")
        document["analysis_lock_sha256"] = blind.sha256_json(body)
        return document

    with pytest.raises(orchestrate.SealOrchestrationError, match="different catalog"):
        orchestrate.build_and_bind_production_bundle(
            catalog_builder=_catalog_builder(calls),
            analysis_builder=wrong_catalog,
            analysis_validator=_analysis_validator,
            _allow_test_mode=True,
        )


def test_transaction_writes_reveal_then_analysis_then_public_with_safe_modes(
    orchestrated_bundle,
    tmp_path: Path,
    monkeypatch,
) -> None:
    public, reveal, analysis = orchestrated_bundle
    checkout = tmp_path / "checkout"
    checkout.mkdir()
    public_path = checkout / "public.json"
    analysis_path = checkout / "analysis.json"
    reveal_path = tmp_path / "private" / "reveal.json"
    link_order: list[str] = []
    real_link = orchestrate.os.link

    def recording_link(source, destination, **kwargs):
        link_order.append(Path(destination).name)
        return real_link(source, destination, **kwargs)

    monkeypatch.setattr(orchestrate.os, "link", recording_link)
    result = orchestrate.write_sealed_bundle_atomically(
        public_path=public_path,
        reveal_path=reveal_path,
        analysis_path=analysis_path,
        checkout_root=checkout,
        public_manifest=public,
        reveal=reveal,
        analysis_lock=analysis,
    )
    assert result == (public_path, reveal_path, analysis_path)
    assert link_order == ["reveal.json", "analysis.json", "public.json"]
    assert stat.S_IMODE(reveal_path.stat().st_mode) == 0o600
    assert stat.S_IMODE(public_path.stat().st_mode) == 0o644
    assert stat.S_IMODE(analysis_path.stat().st_mode) == 0o644
    assert json.loads(public_path.read_text()) == public
    assert json.loads(reveal_path.read_text()) == reveal
    assert json.loads(analysis_path.read_text()) == analysis


def test_transaction_rolls_back_all_targets_on_commit_failure(
    orchestrated_bundle,
    tmp_path: Path,
    monkeypatch,
) -> None:
    public, reveal, analysis = orchestrated_bundle
    checkout = tmp_path / "checkout"
    checkout.mkdir()
    public_path = checkout / "public.json"
    analysis_path = checkout / "analysis.json"
    reveal_path = tmp_path / "private" / "reveal.json"
    real_link = orchestrate.os.link
    count = 0

    def failing_link(source, destination, **kwargs):
        nonlocal count
        count += 1
        if count == 2:
            raise OSError("synthetic commit failure")
        return real_link(source, destination, **kwargs)

    monkeypatch.setattr(orchestrate.os, "link", failing_link)
    with pytest.raises(OSError, match="synthetic"):
        orchestrate.write_sealed_bundle_atomically(
            public_path=public_path,
            reveal_path=reveal_path,
            analysis_path=analysis_path,
            checkout_root=checkout,
            public_manifest=public,
            reveal=reveal,
            analysis_lock=analysis,
        )
    assert not public_path.exists()
    assert not reveal_path.exists()
    assert not analysis_path.exists()
    assert not list(tmp_path.rglob("*.tmp"))


def test_transaction_rejects_overwrite_and_in_checkout_reveal(
    orchestrated_bundle,
    tmp_path: Path,
) -> None:
    public, reveal, analysis = orchestrated_bundle
    checkout = tmp_path / "checkout"
    checkout.mkdir()
    public_path = checkout / "public.json"
    public_path.write_text("owner data")
    with pytest.raises(orchestrate.SealOrchestrationError, match="outside"):
        orchestrate.write_sealed_bundle_atomically(
            public_path=public_path,
            reveal_path=checkout / "reveal.json",
            analysis_path=checkout / "analysis.json",
            checkout_root=checkout,
            public_manifest=public,
            reveal=reveal,
            analysis_lock=analysis,
        )
    with pytest.raises(FileExistsError, match="overwrite"):
        orchestrate.write_sealed_bundle_atomically(
            public_path=public_path,
            reveal_path=tmp_path / "external-reveal.json",
            analysis_path=checkout / "analysis.json",
            checkout_root=checkout,
            public_manifest=public,
            reveal=reveal,
            analysis_lock=analysis,
        )
    assert public_path.read_text() == "owner data"


def test_stdout_summary_contains_only_public_values(orchestrated_bundle, capsys) -> None:
    public, reveal, _ = orchestrated_bundle
    summary = orchestrate.public_seal_summary(public)
    print(json.dumps(summary, sort_keys=True))
    output = capsys.readouterr().out
    assert reveal["secret_hex"] not in output
    assert reveal["opaque_id_key_hex"] not in output
    assert "sealed_payload" not in output
    assert set(summary) == {
        "analysis_lock_sha256",
        "catalog_precommitment_sha256",
        "circuit_count",
        "manifest_sha256",
        "total_planned_shots",
    }


def test_orchestrator_has_no_credential_or_submission_imports() -> None:
    source = (PROGRAMS / "blind_preregister_orchestrate.py").read_text(encoding="utf-8")
    assert "qiskit_ibm_runtime" not in source
    assert "record_gated_cayley_runtime" not in source
    assert ".api-key" not in source
