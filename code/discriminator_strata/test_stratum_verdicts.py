from __future__ import annotations

import json

import pytest

import stratum_verdicts as sv


def test_all_four_verdicts_build_with_expected_statuses() -> None:
    verdicts = sv.build_all()
    assert sorted(verdicts) == [642, 643, 645, 646]
    for issue, verdict in verdicts.items():
        assert verdict["status"] == sv.VERDICTS[issue]
        boundary = verdict["comparison_boundary"]
        assert boundary["public_measurement_read"] is False
        assert boundary["comparison_permitted"] is False
        assert len(verdict["parent_pins"]) == 2
        assert verdict["reopen_condition"]


def test_dyonic_lattice_certificate_is_exact() -> None:
    lattice = sv.dyonic_line_lattice_certificate()
    assert lattice["lattice_isotropic"] is True
    assert lattice["lattice_maximal_by_exhaustion"] is True
    assert lattice["allowed_classes"] == 6
    assert lattice["forbidden_classes"] == 30
    assert lattice["lattice"] == [[0, m] for m in range(6)]


def test_dyonic_pairing_rejects_a_tilted_lattice() -> None:
    def pairing(e1: int, m1: int, e2: int, m2: int) -> int:
        return (e1 * m2 - e2 * m1) % 6

    tilted = [(m % 2, m) for m in range(6)]
    isotropic = all(
        pairing(*left, *right) == 0 for left in tilted for right in tilted
    )
    assert isotropic is False


def test_registry_slot_extraction_counts() -> None:
    registry = sv._load_registry()
    assert len(sv._slot_candidates(registry, "z6_charge_line_congruences")) == 5
    assert len(sv._slot_candidates(registry, "a5_angular_rules")) == 4
    assert len(sv._slot_candidates(registry, "wz_scale_free_response")) == 2
    assert len(sv._slot_candidates(registry, "observer_overlap_cross_spectra")) == 3


def test_committed_verdicts_are_byte_exact() -> None:
    for issue, verdict in sv.build_all().items():
        committed = (sv.RUNTIME / f"stratum_verdict_{issue}.json").read_bytes()
        assert committed == sv.canonical_json_bytes(verdict)


def test_registry_drift_fails_closed(tmp_path, monkeypatch: pytest.MonkeyPatch) -> None:
    tampered = json.loads(sv.REGISTRY_PATH.read_text(encoding="utf-8"))
    tampered["scoring_boundary"]["comparison_access_permitted"] = True
    path = tmp_path / "candidate_registry.json"
    path.write_text(json.dumps(tampered), encoding="utf-8")
    monkeypatch.setattr(sv, "REGISTRY_PATH", path)
    with pytest.raises(sv.VerdictError, match="scoring boundary drift"):
        sv.build_all()


def test_wz_promotion_drift_fails_closed(tmp_path, monkeypatch: pytest.MonkeyPatch) -> None:
    tampered = json.loads(sv.WZ_MANIFEST_PATH.read_text(encoding="utf-8"))
    tampered["promotion_allowed"] = True
    path = tmp_path / "INTEGRATION_MANIFEST.json"
    path.write_text(json.dumps(tampered), encoding="utf-8")
    monkeypatch.setattr(sv, "WZ_MANIFEST_PATH", path)
    registry = sv._load_registry()
    with pytest.raises(sv.VerdictError, match="promotion drift"):
        sv.build_646(registry)
